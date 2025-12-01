# app.py - Final Complete Streamlit App (v5 Stable)

import streamlit as st
import graphviz
from collections import defaultdict

# ============== PAGE CONFIG ==============
st.set_page_config(
    page_title="DNA Pattern Matcher",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# ============== CUSTOM CSS ==============
st.markdown("""
<style>
    .main-header {
        text-align: center;
        background: linear-gradient(90deg, #00d9ff, #00ff88);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        font-size: 2.5rem;
        font-weight: bold;
    }
    .match-highlight {
        background-color: #00ff88;
        color: black;
        padding: 2px 4px;
        border-radius: 3px;
        font-weight: bold;
    }
    .sequence-box {
        background-color: #1a1a2e;
        padding: 15px;
        border-radius: 10px;
        font-family: monospace;
        word-wrap: break-word;
        line-height: 1.8;
    }
    .info-box {
        background-color: #16213e;
        padding: 15px;
        border-radius: 10px;
        border-left: 4px solid #00d9ff;
    }
</style>
""", unsafe_allow_html=True)


# ============== SAMPLE DNA DATA ==============
SAMPLE_SEQUENCES = {
    "Simple Test": "ATGATCGATCGTAGATGCTAGCTGATCGATGTAAATAGCTGATCG",
    "CAG Repeat Test": "ATGAAAGCAGCAGCAGCAGCAGTAAG",
    "BRCA1 Fragment": "ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTT"
}

# =============================================
# AUTOMATA STRUCTURES
# =============================================

class State:
    def __init__(self, id, accept=False):
        self.id = id
        self.accept = accept
        self.transitions = defaultdict(list)
        self.eps = []


class NFA:
    def __init__(self):
        self.states = []
        self.start = None
        self.accept = []
        self.alpha = ['A', 'T', 'G', 'C']

    def add(self, accept=False):
        s = State(len(self.states), accept)
        self.states.append(s)
        if accept:
            self.accept.append(s)
        return s


class DFA:
    def __init__(self):
        self.start = 0
        self.num = 0
        self.accept = []
        self.trans = []
        self.alpha = ['A', 'T', 'G', 'C']


# ========== REGEX → NFA ==========

def pattern_to_nfa(pattern):
    nfa = NFA()
    s0 = nfa.add()
    nfa.start = s0
    cur = s0

    for i, c in enumerate(pattern):
        final = (i == len(pattern)-1)
        nxt = nfa.add(final)
        cur.transitions[c].append(nxt)
        cur = nxt

    return nfa


def alt_nfa(n1, n2):
    n = NFA()
    s = n.add()
    n.start = s

    s.eps.append(n1.start)
    s.eps.append(n2.start)
    n.states = [s] + n1.states + n2.states
    n.accept = n1.accept + n2.accept

    # Fix IDs
    for i, st in enumerate(n.states):
        st.id = i

    return n


def regex_to_nfa(regex):
    parts = regex.split("|")
    nfa = pattern_to_nfa(parts[0])
    for p in parts[1:]:
        nfa = alt_nfa(nfa, pattern_to_nfa(p))
    return nfa


# ========== NFA → DFA ==========

def eps_close(states):
    clk = set(states)
    stack = list(states)
    while stack:
        s = stack.pop()
        for ns in s.eps:
            if ns not in clk:
                clk.add(ns)
                stack.append(ns)
    return frozenset(clk)


def step(states, sym):
    r = set()
    for s in states:
        if sym in s.transitions:
            r.update(s.transitions[sym])
    return r


def nfa_to_dfa(nfa):
    d = DFA()
    M = {}

    def acc(ss):
        return any(s.accept for s in ss)

    start = eps_close([nfa.start])
    M[start] = 0
    if acc(start): d.accept.append(0)

    Q = [start]
    count = 1

    while Q:
        x = Q.pop(0)
        i = M[x]
        for c in d.alpha:
            y = step(x, c)
            if not y: continue
            y = eps_close(y)

            if y not in M:
                M[y] = count
                if acc(y): d.accept.append(count)
                Q.append(y)
                count += 1

            d.trans.append({"from": i, "symbol": c, "to": M[y]})

    d.num = count
    return d


# ========== KMP DFA BUILDING (Optimized Searching) ==========

def kmp(pattern):
    m = len(pattern)
    A = ['A', 'T', 'G', 'C']

    if m == 0:
        return {}, 0

    T = {c: [0]*m for c in A}
    T[pattern[0]][0] = 1
    x = 0

    for j in range(1, m):
        for c in A: T[c][j] = T[c][x]
        T[pattern[j]][j] = j+1
        x = T[pattern[j]][x]

    return T, m


def kmp_find(table, m, text):
    r = []
    s = 0
    for i, c in enumerate(text):
        if c in table:
            s = table[c][s]
        if s == m:
            r.append((i-m+1, i))
            s = 0
    return r


# ========== HIGHLIGHT MATCHES UI ==========

def show_highlights(text, hits):
    if not hits:
        return f'<div class="sequence-box">{text}</div>'
    out = []
    last = -1
    for a, b in hits:
        if a > last+1: out.append(text[last+1:a])
        out.append(f'<span class="match-highlight">{text[a:b+1]}</span>')
        last = b
    if last < len(text)-1:
        out.append(text[last+1:])
    return f'<div class="sequence-box">{"".join(out)}</div>'


# ========== INFO DB ==========

DB = {
    "ATG": ("Start Codon", "Begins translation"),
    "TAA": ("Stop Codon", "Ends protein"),
    "TAG": ("Stop Codon", "Ends protein"),
    "TGA": ("Stop Codon", "Ends protein"),
    "CAG": ("CAG Repeat", "Linked to Huntington’s")
}


# =============================================
# STREAMLIT APP
# =============================================

def main():

    # Sidebar
    with st.sidebar:
        st.header("⚙ Settings")
        choice = st.selectbox("Sample DNA", list(SAMPLE_SEQUENCES.keys()))
        dna_default = SAMPLE_SEQUENCES[choice]

        mode = st.radio("Pattern Type", ["Custom", "ATG", "TAA|TAG|TGA", "CAG"])
        patt = st.text_input("Pattern:", value="ATG" if mode == "Custom" else mode)

    # Title
    st.markdown('<h1 class="main-header">🧬 DNA Pattern Matcher</h1>', unsafe_allow_html=True)
    st.markdown("---")

    dna = st.text_area("Enter DNA Sequence", value=dna_default, height=150).upper()

    if st.button("🔍 Search & Automaton"):
        dna = "".join(c for c in dna if c in "ATGC")
        patt = patt.upper()

        if "|" in patt:
            parts = patt.split("|")
            hits = []
            for p in parts:
                t, m = kmp(p)
                hits += kmp_find(t, m, dna)
        else:
            t, m = kmp(patt)
            hits = kmp_find(t, m, dna)

        st.success(f"Matches Found: {len(hits)}")

        # Tabs
        T1, T2, T3 = st.tabs(["📊 Results", "🔄 Automaton", "📖 Info"])

        with T1:
            st.markdown(show_highlights(dna, hits), unsafe_allow_html=True)

        with T2:
            st.subheader("Deterministic Finite Automaton")

            if "|" in patt:
                nfa = regex_to_nfa(patt)
                dfa = nfa_to_dfa(nfa)
            else:
                nfa = pattern_to_nfa(patt)
                dfa = nfa_to_dfa(nfa)

            dot = graphviz.Digraph()
            dot.attr(rankdir="LR")

            # nodes
            for i in range(dfa.num):
                shape = "doublecircle" if i in dfa.accept else "circle"
                dot.node(f"q{i}", shape=shape)

            # edges
            dot.node("start", shape="none")
            dot.edge("start", "q0")
            for tr in dfa.trans:
                dot.edge(f"q{tr['from']}", f"q{tr['to']}", label=tr["symbol"])

            st.graphviz_chart(dot, use_container_width=True)

        with T3:
            st.subheader("Pattern Information")
            if patt in DB:
                st.info(f"**{DB[patt][0]}** — {DB[patt][1]}")
            else:
                st.warning("Custom Pattern: No biological metadata")



if __name__ == "__main__":
    main()
