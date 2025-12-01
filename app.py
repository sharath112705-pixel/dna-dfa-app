import streamlit as st
import pandas as pd
import re
import graphviz
from collections import defaultdict
from itertools import count

st.set_page_config(page_title="DNA Pattern Detector", layout="wide")

# CSS
st.markdown("""
<style>
.match-highlight {
    background-color: yellow;
    font-weight: bold;
}
.info-box {
    background: #1c1c1c;
    padding: 12px;
    border-radius: 10px;
    margin-top: 10px;
    color: white;
}
</style>
""", unsafe_allow_html=True)

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
    prev = State()
    start = prev

    for char in pattern:
        new = State()
        prev.transitions[char].add(new)
        prev = new

    prev.is_accept = True
    return NFA(start, prev)


class DFA:
    def __init__(self):
        self.states = {}
        self.start_state = None
        self.accept_states = set()
        self.transitions = defaultdict(dict)

    def add_transition(self, src, symbol, dest):
        self.transitions[src][symbol] = dest

    def next_state(self, state, symbol):
        return self.transitions.get(state, {}).get(symbol)


def nfa_to_dfa(nfa):
    alphabet = {"A", "T", "C", "G"}
    dfa = DFA()

    start = frozenset([nfa.start])
    dfa.start_state = start
    dfa.states[start] = True

    unprocessed = [start]

    while unprocessed:
        current = unprocessed.pop()
        for sym in alphabet:
            nxt = set()
            for s in current:
                nxt.update(s.transitions.get(sym, []))

            nxt = frozenset(nxt)
            if not nxt:
                continue

            if nxt not in dfa.states:
                dfa.states[nxt] = True
                unprocessed.append(nxt)

            dfa.add_transition(current, sym, nxt)

    for s in dfa.states:
        if any(st.is_accept for st in s):
            dfa.accept_states.add(s)

    return dfa


def search_overlapping(dfa, seq):
    matches = []
    n = len(seq)

    for i in range(n):
        curr = dfa.start_state
        for j in range(i, n):
            curr = dfa.next_state(curr, seq[j])
            if not curr:
                break
            if curr in dfa.accept_states:
                matches.append({
                    "start": i,
                    "end": j + 1,
                    "sequence": seq[i:j + 1]
                })
    return matches


def draw_dfa(dfa):
    dot = graphviz.Digraph()
    dot.attr(rankdir="LR")

    for s in dfa.states:
        dot.node(str(id(s)), shape="doublecircle" if s in dfa.accept_states else "circle")

    added = set()
    for src, trans in dfa.transitions.items():
        for sym, dst in trans.items():
            key = (id(src), sym, id(dst))
            if key not in added:
                added.add(key)
                dot.edge(str(id(src)), str(id(dst)), label=sym)

    return dot


def highlight(seq, matches):
    seq_list = list(seq)
    for m in matches:
        for i in range(m["start"], m["end"]):
            seq_list[i] = f"<span class='match-highlight'>{seq_list[i]}</span>"
    return "".join(seq_list)


def show():
    st.title("🧬 DNA Pattern Detection using Automata")

    seq = st.text_area("Enter DNA Sequence (A,T,G,C only):").upper()
    pattern = st.text_input("Enter Pattern:").upper()

    if st.button("🔍 Search"):

        if not re.fullmatch("[ATGC]+", seq):
            st.error("❌ Only A,T,G,C allowed in DNA")
            st.stop()

        if not re.fullmatch("[ATGC]+", pattern):
            st.error("❌ Pattern invalid. Only A,T,G,C supported.")
            st.stop()

        nfa = regex_to_nfa(pattern)
        dfa = nfa_to_dfa(nfa)
        matches = search_overlapping(dfa, seq)

        tab1, tab2, tab3 = st.tabs(["Sequence", "DFA", "Matches"])

        with tab1:
            st.markdown(highlight(seq, matches), unsafe_allow_html=True)

        with tab2:
            c1, c2, c3 = st.columns(3)
            c1.metric("DFA States", len(dfa.states))
            c2.metric("Transitions", sum(len(v) for v in dfa.transitions.values()))
            c3.metric("Accept States", len(dfa.accept_states))
            st.graphviz_chart(draw_dfa(dfa))

        with tab3:
            if not matches:
                st.warning("No match found")
            else:
                df = pd.DataFrame(matches)
                st.dataframe(df)


if __name__ == "__main__":
    show()
