import streamlit as st
import re
from dfa_graph import plot_dna_dfa

st.set_page_config(page_title="DNA Pattern DFA", layout="wide")

st.title("🧬 DNA Pattern Detection using DFA")
st.markdown("### Detecting Start Codon Pattern: **ATG**")

dna_input = st.text_input("Enter DNA Sequence (A, T, G, C only):").upper()

if dna_input:
    # Validation
    if not re.fullmatch("[ATGC]+", dna_input):
        st.error("❌ Invalid sequence! Only A, T, G, C characters allowed.")
    else:
        st.subheader("🔹 DFA State Diagram")
        graph = plot_dna_dfa()
        st.graphviz_chart(graph)

        st.subheader("🔹 Pattern Search Results")

        pattern = "ATG"
        matches = [(m.start(), m.group()) for m in re.finditer(pattern, dna_input)]

        if matches:
            st.success(f"🎯 Pattern Found! {len(matches)} match(es) of '{pattern}'")

            for idx, (pos, seq) in enumerate(matches, start=1):
                st.write(f"Match {idx}: Position {pos} → {seq}")

            # Highlight matches
            highlighted_seq = dna_input.replace(pattern, f":green[{pattern}]")
            st.markdown(f"### ✨ Highlighted Sequence:\n{highlighted_seq}")
        else:
            st.warning("⚠ No ATG pattern detected in sequence.")
