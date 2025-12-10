# app.py - Upgraded DNA Pattern Matcher with fixes for IndexError and sidebar field updates
import streamlit as st
import graphviz
from collections import defaultdict, deque
import pandas as pd
import time
import io
import base64

# Try to import reportlab for PDF export
try:
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfgen import canvas
    REPORTLAB_AVAILABLE = True
except Exception:
    REPORTLAB_AVAILABLE = False

# ============== PAGE CONFIG & CSS (logo + watermark) ==============
st.set_page_config(page_title="DNA Pattern Matcher", page_icon="🧬", layout="wide")

# ============== SAMPLE DATA ==============
SAMPLE_SEQUENCES = {
    'Simple Test': 'ATGATCGATCGTAGATGCTAGCTGATCGATGTAAATAGCTGATCG',
    'BRCA1 Fragment': 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTG',
    'CAG Repeat Test': 'ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCG'
}

# ================= SESSION STATE FOR ANIMATION =================
if 'anim_index' not in st.session_state:
    st.session_state.anim_index = 0
if 'anim_steps' not in st.session_state:
    st.session_state.anim_steps = []
if 'anim_playing' not in st.session_state:
    st.session_state.anim_playing = False
if 'anim_last_time' not in st.session_state:   # ✅ NEW
    st.session_state.anim_last_time = 0        # ✅ NEW

# ================= KMP =================
def build_kmp_automaton(pattern):
    m = len(pattern)
    alphabet = ['A', 'T', 'G', 'C']
    if m == 0:
        return {c: [0] for c in alphabet}, 0
    dfa_table = {c: [0] * m for c in alphabet}
    dfa_table[pattern[0]][0] = 1
    x = 0
    for j in range(1, m):
        for c in alphabet:
            dfa_table[c][j] = dfa_table[c][x]
        dfa_table[pattern[j]][j] = j + 1
        x = dfa_table[pattern[j]][x]
    return dfa_table, m

def search_kmp(dfa_table, m, text):
    matches = []
    state = 0
    for i, ch in enumerate(text):
        if ch not in dfa_table:
            state = 0
            continue
        state = dfa_table[ch][state]
        if state == m:
            matches.append((i - m + 1, i))
            state = 0
    return matches

# ================= DFA VISUAL =================
def create_dfa_graph(dfa_struct, highlight_state=None):
    dot = graphviz.Digraph(engine='dot')
    dot.attr(rankdir='LR')
    for i in range(dfa_struct['num']):
        fill = '#0ea5a4' if i in dfa_struct['accepts'] else '#1e293b'
        if highlight_state is not None and i == highlight_state:
            dot.node(f"q{i}", f"q{i}", style='filled', fillcolor='#f97316')
        else:
            dot.node(f"q{i}", f"q{i}", style='filled', fillcolor=fill)
    dot.node('start', shape='none')
    dot.edge('start', 'q0')
    for t in dfa_struct['trans']:
        dot.edge(f"q{t['from']}", f"q{t['to']}", label=t['symbol'])
    return dot

# ================= MAIN =================
st.title("🧬 DNA Pattern Matcher")

dna = st.text_area("Enter DNA:", SAMPLE_SEQUENCES['Simple Test'])
pattern = st.text_input("Enter Pattern:", "ATG")

if st.button("Search"):
    table, m = build_kmp_automaton(pattern)
    results = search_kmp(table, m, dna)

    st.success(f"Matches Found: {len(results)}")

    if results:
        start, end = results[0]
        subseq = dna[start:end+1]

        steps = []
        state = 0
        for ch in subseq:
            state = table[ch][state]
            steps.append({"char": ch, "state": state})

        st.session_state.anim_steps = steps
        st.session_state.anim_index = 0

# ================= ✅ FIXED TRAVERSAL ANIMATION =================
st.subheader("🎬 Traversal Animation")

if st.session_state.anim_steps:
    steps = st.session_state.anim_steps

    col1, col2, col3, col4 = st.columns(4)
    delay = col1.slider("Delay", 0.1, 1.5, 0.5)

    if col2.button("◀ Back"):
        st.session_state.anim_index = max(0, st.session_state.anim_index - 1)
        st.session_state.anim_playing = False

    if col3.button("Next ▶"):
        st.session_state.anim_index = min(len(steps)-1, st.session_state.anim_index + 1)
        st.session_state.anim_playing = False

    if col4.button("⏯ Play / Pause"):
        st.session_state.anim_playing = not st.session_state.anim_playing

    # ✅ ✅ ✅ REAL STEP-BY-STEP ANIMATION
    if st.session_state.anim_playing:
        now = time.time()
        if now - st.session_state.anim_last_time >= delay:
            st.session_state.anim_last_time = now
            if st.session_state.anim_index < len(steps) - 1:
                st.session_state.anim_index += 1
                st.rerun()
            else:
                st.session_state.anim_playing = False

    idx = st.session_state.anim_index
    st.info(f"Current State: q{steps[idx]['state']} | Character: {steps[idx]['char']}")
