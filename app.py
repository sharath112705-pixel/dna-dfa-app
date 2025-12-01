import streamlit as st
import pandas as pd
import re
import graphviz
from collections import defaultdict
from itertools import count

st.set_page_config(page_title="DNA DFA Pattern Detector", layout="wide")

# ----------------------------------
# TITLE & CUSTOM CSS
# ----------------------------------
st.markdown("""
    <h1 style="text-align:center;margin-bottom:0;">🧬 DNA Pattern Detection using DFA</h1>
    <p style="text-align:center;margin-top:5px;">Regex → NFA → DFA → Match in DNA Sequence</p>
    <style>
        .match-highlight { background-color: yellow; font-weight: bold; }
        .info-box {
            background: #1c1c1c;
            padding: 12px;
            border-radius: 10px;
            margin-top: 10px;
            color: white;
        }
    </style>
""", unsafe_allow_html=True)


# ----------------------------------
# NFA & DFA CLASSES
# ----------------------------------
class State:
    _id = count()

    def __init__(self, is_accept=False):
        self.id = next(self._id)
        self.transitions = defaultdict(set)
        self.is_accept = is_accept


class NFA:
    def __init__(self, start, accept):
        self.start = start
        self.accept = accept


def regex_to_nfa(pattern):
    prev_state = State()
    start = prev_state

    for char in pattern:
        new_state = State()
        prev_state.transitions[char].add(new_state)
        prev_state = new_state

    prev_state.is_accept = True
    return NFA(start, prev_state)


class DFA:
    def __init__(self):
        self.states = {}
        self.start_state = None
        self.accept_states = set()
        self.transitions = defaultdict(dict)

    def add_transition(self, src, symbol, dest):
        self.transitions[src][symbol] = dest

    def get_next_state(self, state, symbol):
        return self.transitions.get(state, {}).get(symbol)


def nfa_to_dfa(nfa):
    dfa = DFA()
    initial = frozenset([nfa.start])
    dfa.start_state = initial
    dfa.states[initial] = True

    unmarked = [initial]
    alphabet = {"A", "T", "G", "C"}

    while unmarked:
        current = unmarked.pop()
        for symbol in alphabet:
            next_set = set()
            for state in current:
                next_set.update(state.transitions.get(symbol, []))

            next_state = frozenset(next_set)
            if not next_state:
                continue

            if next_state not in dfa.states:
                dfa.states[next_state] = True
                unmarked.append(next_state)

            dfa.add_transition(current, symbol, next_state)

    # Identify accepting states
    for s in dfa.states:
        if any(st.is_accept for st in s):
            dfa.accept_states.add(s)

    return dfa


# ----------------------------------
# SEARCH MATCHES (supports overlapping)
# ----------------------------------
def search_pattern(dfa, seq):
    matches = []
    n = len(seq)

    for start in range(n):
        current = dfa.start_state
        for end in range(start, n):
            char = seq[end]
            current = dfa.get_next_state(current, char)
            if not current:
                break
            if current in dfa.accept_states:
                matches.append((start, end + 1))

    return matches


# ----------------------------------
# GRAPHVIZ DFA VISUAL
# ----------------------------------
def draw_dfa(dfa):
    dot = graphviz.Digraph()
    dot.attr(rankdir='LR')

    for s in dfa.states:
        dot.node(str(id(s)),
                 shape="doublecircle" if s in dfa.accept_states else "circle")

    edges_added = set()
    for src, trans in dfa.transitions.items():
        for char, dest in trans.items():
            edge_key = (id(src), char, id(dest))
            if edge_key not in edges_added:
                dot.edge(str(id(src)), str(id(dest)), label=char)
                edges_added.add(edge_key)

    return dot


# ----------------------------------
# UI INPUTS
# ----------------------------------
seq = st.text_area("Enter DNA Sequence (A,T,G,C only):").upper()
pattern = st.text_input("Enter Pattern (simple regex): ").upper()

if st.button("🔍 Search Now"):
    if not re.fullmatch("[ATGC]+", seq):
        st.error("❌ Invalid DNA sequence! Only A,T,G,C allowed.")
        st.stop()

    if not re.fullmatch("[ATGC]+", pattern):
        st.error("❌ Pattern must only contain A,T,G,C (no special regex yet)")
        st.stop()

    # Processing pipeline
    nfa = regex_to_nfa(pattern)
    dfa = nfa_to_dfa(nfa)
    matches = search_pattern(dfa, seq)

    # Tabs
    tab1, tab2, tab3 = st.tabs(["📍 Highlight", "🚀 DFA Graph", "📋 Match List"])

    # TAB 1: Highlight in Sequence
    with tab1:
        display_seq = seq
        for start, end in sorted(matches, reverse=True):
            display_seq = (display_seq[:start] +
                           f"<span class='match-highlight'>{display_seq[start:end]}</span>" +
                           display_seq[end:])
        st.markdown(display_seq, unsafe_allow_html=True)

    # TAB 2: Graphviz DFA
    with tab2:
        st.graphviz_chart(draw_dfa(dfa))

    # TAB 3: List Matches
    with tab3:
        if matches:
            df = pd.DataFrame(matches, columns=["Start", "End"])
            df["Matched Pattern"] = [seq[s:e] for s, e in matches]
            st.dataframe(df)
        else:
            st.warning("No matches found ❌")


st.info("📌 Example: Sequence = ATGCATGC, Pattern = ATG")
