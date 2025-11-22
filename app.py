import streamlit as st
from dfa_graph.py import plot_dna_dfa
import re

st.set_page_config(page_title="DNA Pattern DFA", layout="wide")

st.title("🧬 DNA Mutation Pattern Detector using DFA")
st.markdown("### Pattern: `ATG` (Start Codon)")

dna_input = st.text_input("Enter DNA Sequence (A,T,G,C):").upper()

if dna_input:
    if not re.fullmatch("[ATGC]+", dna_input):
        st.error("Invalid! Only A, T, G, C are allowed.")
    else:
        st.subheader("🔹 DFA Visualization")
        graph = plot_dna_dfa()
        st.graphviz_chart(graph)

        st.subheader("🔹 Pattern Matching Results")
        pattern = "ATG"
        matches = [(m.start(), m.group()) for m in re.finditer(pattern, dna_input)]

        if matches:
            st.success(f"Found {len(matches)} match(es) of '{pattern}'")

            for i, (pos, seq) in enumerate(matches, 1):
                st.write(f"Match {i}: Position {pos}, Sequence: {seq}")

            highlighted = dna_input.replace(pattern, f":green[{pattern}]")
            st.markdown(f"### Result with Highlight:\n{highlighted}")
        else:
            st.warning("No mutation pattern found!")
