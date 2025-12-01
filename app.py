# app.py - Upgraded DNA Pattern Matcher with animation, NFA visual, CSV/PDF download, colors, performance
# Replace your existing file with this complete script.

import streamlit as st
import graphviz
from collections import defaultdict, deque
import pandas as pd
import time
import io
import base64
import sys

# Try to import reportlab for PDF export
try:
    from reportlab.lib.pagesizes import letter
    from reportlab.pdfgen import canvas
    REPORTLAB_AVAILABLE = True
except Exception:
    REPORTLAB_AVAILABLE = False

# ============== PAGE CONFIG & CSS (logo + watermark) ==============
st.set_page_config(page_title="DNA Pattern Matcher (Upgraded)", page_icon="🧬", layout="wide")

st.markdown("""
<style>
/* Header logo as inline SVG */
.header {
  display:flex;
  align-items:center;
  gap:12px;
  margin-bottom: 8px;
}
.logo-svg { width:56px; height:56px; filter: drop-shadow(0 2px 6px rgba(0,0,0,0.3)); }
.watermark {
  position: fixed;
  right: 10px;
  bottom: 10px;
  opacity: 0.06;
  font-size: 56px;
  transform: rotate(-20deg);
  z-index: 0;
  pointer-events:none;
}
.main-header {
  font-size: 28px;
  font-weight: 700;
  margin: 0;
}
.match-highlight {
    padding:2px 4px;
    border-radius:4px;
    font-weight:700;
    color:#000;
}
.sequence-box {
    background-color: #0f1724;
    padding: 12px;
    border-radius: 8px;
    font-family: monospace;
    color: #e6eef8;
    white-space: pre-wrap;
    word-break: break-word;
}
.info-box {
    background-color: #0b1220;
    padding:12px;
    border-radius:8px;
    border-left:4px solid #10b981;
    color:#dbeafe;
}
.small-muted { color: #94a3b8; font-size: 0.9rem; }
.controls { display:flex; gap:8px; align-items:center; }
</style>
<div class="header">
  <!-- inline SVG logo -->
  <svg class="logo-svg" viewBox="0 0 64 64" xmlns="http://www.w3.org/2000/svg">
    <defs>
      <linearGradient id="g1" x1="0" x2="1">
        <stop offset="0" stop-color="#00d9ff"/>
        <stop offset="1" stop-color="#00ff88"/>
      </linearGradient>
    </defs>
    <rect width="64" height="64" rx="12" fill="#071029"/>
    <g transform="translate(8,8)">
      <path d="M6 0 C10 6, 18 6, 22 0" stroke="url(#g1)" stroke-width="2.6" fill="none" />
      <circle cx="6" cy="24" r="6" fill="url(#g1)"/>
      <circle cx="22" cy="24" r="6" fill="#0ea5a4"/>
    </g>
  </svg>
  <div>
    <div class="main-header">🧬 DNA Pattern Matcher — Upgraded</div>
    <div class="small-muted">Automata visualization • fast multi-pattern search • traversal animation</div>
  </div>
</div>
<div class="watermark">DNA</div>
""", unsafe_allow_html=True)

# ============== SAMPLE DATA ==============
SAMPLE_SEQUENCES = {
    'Simple Test': 'ATGATCGATCGTAGATGCTAGCTGATCGATGTAAATAGCTGATCG',
    'BRCA1 Fragment': 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTG',
    'CAG Repeat Test': 'ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCG'
}

# ============== AUTOMATA / SEARCH IMPLEMENTATIONS ==============

# --- Aho-Corasick implementation (multi-pattern) for performance ---
class AhoNode:
    __slots__ = ('next', 'fail', 'out')
    def __init__(self):
        self.next = dict()
        self.fail = None
        self.out = []

class Aho:
    def __init__(self, patterns):
        self.root = AhoNode()
        self._build_trie(patterns)
        self._build_fail()

    def _build_trie(self, patterns):
        for pid, p in enumerate(patterns):
            node = self.root
            for ch in p:
                node = node.next.setdefault(ch, AhoNode())
            node.out.append(pid)

    def _build_fail(self):
        q = deque()
        for ch, node in self.root.next.items():
            node.fail = self.root
            q.append(node)
        while q:
            r = q.popleft()
            for ch, s in r.next.items():
                q.append(s)
                state = r.fail
                while state is not None and ch not in state.next:
                    state = state.fail
                s.fail = state.next[ch] if state and ch in state.next else self.root
                s.out += s.fail.out

    def search(self, text):
        node = self.root
        results = []  # list of (end_index, pattern_id)
        for i, ch in enumerate(text):
            while node is not None and ch not in node.next:
                node = node.fail
            node = node.next[ch] if node and ch in node.next else self.root
            for pid in node.out:
                results.append((i, pid))
        return results

# --- KMP for single pattern ---
def build_kmp_automaton(pattern):
    m = len(pattern)
    alphabet = ['A', 'T', 'G', 'C']
    if m == 0:
        return {c: [0] * 0 for c in alphabet}, 0
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

# ============== NFA/DFA constructions (for visualization) ==============

class StateV:
    def __init__(self, idx, accept=False):
        self.id = idx
        self.accept = accept
        self.trans = defaultdict(list)
        self.eps = []

class NFAVis:
    def __init__(self):
        self.states = []
        self.start = None

    def add(self, accept=False):
        s = StateV(len(self.states), accept)
        self.states.append(s)
        return s

def pattern_to_nfa_vis(pattern):
    n = NFAVis()
    if not pattern:
        s = n.add(True)
        n.start = s
        return n
    cur = n.add()
    n.start = cur
    for i, ch in enumerate(pattern):
        final = (i == len(pattern)-1)
        nxt = n.add(final)
        cur.trans[ch].append(nxt)
        cur = nxt
    return n

def regex_to_nfa_vis(regex):
    if '|' in regex:
        parts = regex.split('|')
        # build combined manually by making new start with eps to each
        combined = NFAVis()
        start = combined.add()
        combined.start = start
        offset = 1
        for part in parts:
            sub = pattern_to_nfa_vis(part)
            # shift ids for sub
            base = len(combined.states)
            # append states
            for s in sub.states:
                new = StateV(len(combined.states), s.accept)
                combined.states.append(new)
            # copy transitions with adjusted ids
            for s in sub.states:
                sid = base + s.id
                new_s = combined.states[sid]
                for ch, dests in s.trans.items():
                    for d in dests:
                        new_s.trans[ch].append(combined.states[base + d.id])
                for e in s.eps:
                    new_s.eps.append(combined.states[base + e.id])
            # connect start eps to the sub start
            combined.states[0].eps.append(combined.states[base + sub.start.id])
        return combined
    return pattern_to_nfa_vis(regex)

def nfa_vis_to_dfa(nfa_vis):
    # simple powerset conversion for visualization transitions (states are frozensets of NFA states)
    alphabet = ['A', 'T', 'G', 'C']
    # epsilon-closure
    def eps_close(states):
        stack = list(states)
        out = set(states)
        while stack:
            s = stack.pop()
            for e in s.eps:
                if e not in out:
                    out.add(e)
                    stack.append(e)
        return frozenset(out)
    start = eps_close([nfa_vis.start])
    mapping = {start: 0}
    q = [start]
    trans = []
    accepts = set()
    idx = 1
    if any(s.accept for s in start): accepts.add(0)
    while q:
        cur = q.pop(0)
        cur_id = mapping[cur]
        for ch in alphabet:
            nxtset = set()
            for s in cur:
                nxtset.update(s.trans.get(ch, []))
            if not nxtset: continue
            nxt = eps_close(nxtset)
            if nxt not in mapping:
                mapping[nxt] = idx
                if any(s.accept for s in nxt): accepts.add(idx)
                idx += 1
                q.append(nxt)
            trans.append({'from': cur_id, 'to': mapping[nxt], 'symbol': ch})
    dfa = {'num': idx, 'trans': trans, 'accepts': accepts}
    return dfa

# ============== UI helpers ==============

# Pattern colors (text highlights). Option 3: text colored, graphs unified.
PATTERN_COLORS = ['#A7F3D0', '#FDE68A', '#C7B3FF', '#FBCFE8', '#93C5FD', '#FDBA74']

def colored_highlight_text(dna, matches_by_pattern):
    """
    matches_by_pattern: list of lists of (start,end) per pattern
    produce a string with colored spans; patterns earlier in list get colors[0], etc.
    Overlaps: priority to earlier patterns in the list.
    """
    n = len(dna)
    owner = [None] * n
    for pid, matches in enumerate(matches_by_pattern):
        for (a,b) in matches:
            for i in range(a, b+1):
                if owner[i] is None:
                    owner[i] = pid
    parts = []
    i = 0
    while i < n:
        if owner[i] is None:
            parts.append(dna[i])
            i += 1
        else:
            pid = owner[i]
            color = PATTERN_COLORS[pid % len(PATTERN_COLORS)]
            j = i
            while j < n and owner[j] == pid:
                j += 1
            # span from i to j-1
            segment = dna[i:j]
            parts.append(f"<span class='match-highlight' style='background:{color};'>{segment}</span>")
            i = j
    return "<div class='sequence-box'>" + "".join(parts) + "</div>"

# CSV generator
def df_to_csv_bytes(df):
    return df.to_csv(index=False).encode('utf-8')

# PDF report using reportlab
def create_pdf_report(title, dna, patterns, matches_by_pattern):
    """
    Patterns list and matches_by_pattern: list of lists of (start,end)
    """
    if not REPORTLAB_AVAILABLE:
        return None
    buffer = io.BytesIO()
    c = canvas.Canvas(buffer, pagesize=letter)
    width, height = letter
    margin = 40
    y = height - margin

    c.setFont("Helvetica-Bold", 16)
    c.drawString(margin, y, title); y -= 20
    c.setFont("Helvetica", 10)
    c.drawString(margin, y, f"Sequence length: {len(dna)} bp"); y -= 16
    c.drawString(margin, y, f"Patterns: {', '.join(patterns)}"); y -= 20

    # Add first 600 chars of DNA with highlights as plain text (reportlab doesn't do HTML easily)
    c.setFont("Courier", 8)
    display = dna[:800]
    # naive: put markers where matches are
    y -= 10
    c.drawString(margin, y, "Sequence (first 800 bp):"); y -= 14
    lines = [display[i:i+100] for i in range(0, len(display), 100)]
    for ln in lines:
        c.drawString(margin, y, ln)
        y -= 12
        if y < margin + 50:
            c.showPage()
            y = height - margin
    # Add simple match summary
    c.showPage()
    y = height - margin
    c.setFont("Helvetica-Bold", 12)
    c.drawString(margin, y, "Match Summary"); y -= 18
    c.setFont("Helvetica", 10)
    for pid, matches in enumerate(matches_by_pattern):
        c.drawString(margin, y, f"Pattern '{patterns[pid]}': {len(matches)} matches")
        y -= 14
        if y < margin + 40:
            c.showPage()
            y = height - margin
    c.save()
    buffer.seek(0)
    return buffer.read()

# Graphviz rendering helpers
def create_kmp_dfa_visualization(pattern):
    if not pattern:
        return None
    dfa_table, _ = build_kmp_automaton(pattern)
    m = len(pattern)
    alphabet = ['A','T','G','C']
    dot = graphviz.Digraph(comment='KMP DFA', engine='dot')
    dot.attr(rankdir='LR', bgcolor='#0f1724', fontcolor='white')
    dot.attr('node', style='filled', fontcolor='white')
    for i in range(m+1):
        if i == m:
            dot.node(f"q{i}", f"q{i}", shape='doublecircle', fillcolor='#065f46')
        elif i == 0:
            dot.node(f"q{i}", f"q{i}", shape='circle', fillcolor='#1e40af')
        else:
            dot.node(f"q{i}", f"q{i}", shape='circle', fillcolor='#1e293b')
    dot.node('', shape='none', width='0')
    dot.edge('', 'q0', color='#ef4444')
    edge_map = defaultdict(list)
    for state in range(m):
        for char in alphabet:
            ns = dfa_table[char][state]
            edge_map[(state, ns)].append(char)
    for (frm, to), chars in edge_map.items():
        lbl = ','.join(sorted(chars))
        dot.edge(f"q{frm}", f"q{to}", label=lbl)
    return dot

def create_nfa_graph(nfa_vis):
    dot = graphviz.Digraph(engine='dot')
    dot.attr(rankdir='LR')
    dot.attr('node', style='filled', fontcolor='white')
    # nodes
    for s in nfa_vis.states:
        if s.accept:
            dot.node(f"s{s.id}", f"s{s.id}", shape='doublecircle', fillcolor='#065f46')
        else:
            dot.node(f"s{s.id}", f"s{s.id}", shape='circle', fillcolor='#1e293b')
    dot.node('start', shape='none')
    dot.edge('start', f"s{nfa_vis.start.id}")
    # transitions
    for s in nfa_vis.states:
        for ch, dests in s.trans.items():
            for d in dests:
                dot.edge(f"s{s.id}", f"s{d.id}", label=ch)
        for e in s.eps:
            dot.edge(f"s{s.id}", f"s{e.id}", label='ε', style='dashed')
    return dot

def create_dfa_graph(dfa_struct, highlight_state=None):
    # dfa_struct: dict with 'num','trans' (list), 'accepts' (set)
    dot = graphviz.Digraph(engine='dot')
    dot.attr(rankdir='LR')
    for i in range(dfa_struct['num']):
        fill = '#0ea5a4' if i in dfa_struct['accepts'] else '#1e293b'
        if highlight_state is not None and i == highlight_state:
            dot.node(f"q{i}", f"q{i}", style='filled', fillcolor='#f97316')  # highlight color
        else:
            dot.node(f"q{i}", f"q{i}", style='filled', fillcolor=fill)
    dot.node('start', shape='none')
    dot.edge('start', 'q0')
    for t in dfa_struct['trans']:
        dot.edge(f"q{t['from']}", f"q{t['to']}", label=t['symbol'])
    return dot

# ============== App state setup ==============
if 'autoplay' not in st.session_state:
    st.session_state.autoplay = False
if 'anim_index' not in st.session_state:
    st.session_state.anim_index = 0
if 'anim_steps' not in st.session_state:
    st.session_state.anim_steps = []
if 'anim_playing' not in st.session_state:
    st.session_state.anim_playing = False
if 'last_search_id' not in st.session_state:
    st.session_state.last_search_id = None

# ============== MAIN APP UI ==============
def main():
    st.sidebar.header("⚙️ Settings")
    sample_choice = st.sidebar.selectbox("Choose sample:", list(SAMPLE_SEQUENCES.keys()))
    quick_pattern = st.sidebar.radio("Quick patterns:", ['Custom', 'ATG', 'TAA|TAG|TGA', 'TATA', 'GAATTC', 'CAG', 'AATAAA'])
    st.sidebar.markdown("---")
    st.sidebar.info("TOC: NFA/DFA, KMP, Aho-Corasick (fast multi-pattern search)")

    st.markdown("### Sequence & Pattern Input")
    col1, col2 = st.columns([3,1])
    with col1:
        dna_input = st.text_area("Enter DNA sequence (A,T,G,C only):", value=SAMPLE_SEQUENCES[sample_choice], height=170, key='dna_text')
    with col2:
        if quick_pattern == 'Custom':
            pattern_input = st.text_input("Enter pattern:", value='ATG', key='pattern_field')
        else:
            pattern_input = st.text_input("Enter pattern:", value=quick_pattern, key='pattern_field')
        st.caption("Example patterns: ATG, TATA, TAA|TAG|TGA")

    # Preprocess
    dna_clean = ''.join([c for c in dna_input.upper() if c in 'ATGC'])
    pattern_clean = pattern_input.upper().strip()
    st.write("")  # spacing

    # Search options
    st.markdown("### Search Options")
    col_a, col_b, col_c = st.columns([1,1,1])
    with col_a:
        use_aho = st.checkbox("Use Aho–Corasick for multi-pattern (fast)", value=True)
    with col_b:
        enable_animation = st.checkbox("Enable traversal animation (manual + autoplay)", value=True)
    with col_c:
        max_matches = st.number_input("Max matches to show per pattern", min_value=10, max_value=10000, value=500)

    # --- Run Search ---
    run = st.button("🔍 Build Automaton & Search")

    if run:
        if not dna_clean:
            st.error("Invalid or empty DNA sequence after filtering. Use characters A,T,G,C.")
            st.stop()
        if not pattern_clean:
            st.error("Please enter a pattern.")
            st.stop()

        # handle alternation
        if '|' in pattern_clean:
            patterns = [p for p in (pattern_clean.split('|')) if p]
        else:
            patterns = [pattern_clean]

        # decide search method
        matches_by_pattern = []
        # Use Aho for multiple patterns with reasonable lengths
        if use_aho and len(patterns) > 1:
            aho = Aho(patterns)
            results = aho.search(dna_clean)
            # results are (end_index, pid) -> convert to (start,end)
            pid_to_matches = defaultdict(list)
            for end_idx, pid in results:
                plen = len(patterns[pid])
                start_idx = end_idx - plen + 1
                if len(pid_to_matches[pid]) < max_matches:
                    pid_to_matches[pid].append((start_idx, end_idx))
            matches_by_pattern = [pid_to_matches[i] for i in range(len(patterns))]
        else:
            # KMP per pattern
            for p in patterns:
                table, m = build_kmp_automaton(p)
                found = search_kmp(table, m, dna_clean)
                if max_matches:
                    found = found[:max_matches]
                matches_by_pattern.append(found)

        total_matches = sum(len(x) for x in matches_by_pattern)
        st.success(f"Search complete — {total_matches} matches found across {len(patterns)} pattern(s).")

        # prepare dataframe for download & display
        rows = []
        for pid, matches in enumerate(matches_by_pattern):
            for (a,b) in matches:
                rows.append({'Pattern': patterns[pid], 'Start': a, 'End': b, 'Sequence': dna_clean[a:b+1], 'PatternID': pid})
        df = pd.DataFrame(rows)

        # store search data in session for animation
        st.session_state.search_data = {
            'dna': dna_clean,
            'patterns': patterns,
            'matches_by_pattern': matches_by_pattern,
            'df': df
        }
        # reset animation indices
        st.session_state.anim_index = 0
        st.session_state.anim_steps = []
        st.session_state.anim_playing = False
        st.session_state.last_search_id = time.time()

    # If session has search data, show results
    if st.session_state.get('search_data') is not None:
        data = st.session_state['search_data']
        dna = data['dna']
        patterns = data['patterns']
        matches_by_pattern = data['matches_by_pattern']
        df = data['df']

        # Tabs: Results, DFA/NFA, Details, Animation/Traversal, Download
        tab_results, tab_auto, tab_nfa, tab_details, tab_anim, tab_download = st.tabs(["📊 Results", "🔄 DFA (KMP)", "🧩 NFA Visual", "📋 Match Details", "🎬 Traversal Animation", "📥 Export"])

        with tab_results:
            st.subheader("Highlighted Matches (text only, colored per pattern)")
            html = colored_highlight_text(dna, matches_by_pattern)
            st.markdown(html, unsafe_allow_html=True)
            st.markdown("---")
            st.info(f"Patterns: {', '.join(patterns)}  —  Total matches: {len(df)}")
            st.caption("Pattern highlight colors are only for text output; automaton graphs use unified coloring for clarity.")

        with tab_auto:
            st.subheader("Deterministic Automaton (constructed from selected pattern)")
            # For DFA visualization, if alternation present we'll convert full NFA->DFA (powerset)
            # Here show KMP DFA for single pattern or NFA->DFA for alternation
            if len(patterns) == 1:
                dot = create_kmp_dfa_visualization(patterns[0])
                st.caption("KMP DFA (All transitions shown).")
                if dot:
                    st.graphviz_chart(dot, use_container_width=True)
            else:
                # build NFA visual and compute DFA struct for visualization (unified)
                nfa_vis = regex_to_nfa_vis("|".join(patterns))
                dfa_struct = nfa_vis_to_dfa(nfa_vis)
                dot = create_dfa_graph(dfa_struct)
                st.caption("DFA (from alternation NFA via subset construction); unified colors.")
                st.graphviz_chart(dot, use_container_width=True)

        with tab_nfa:
            st.subheader("NFA Visualization (structure with ε transitions)")
            nfa_vis = regex_to_nfa_vis("|".join(patterns))
            dotn = create_nfa_graph(nfa_vis)
            st.caption("NFA graph (states from left to right). Epsilon transitions are dashed.")
            st.graphviz_chart(dotn, use_container_width=True)

        with tab_details:
            st.subheader("Match Table")
            if df.empty:
                st.warning("No matches to show.")
            else:
                # show table with selection
                # allow user to select a match to step through
                st.dataframe(df[['Pattern', 'Start', 'End', 'Sequence']])
                st.markdown("---")
                st.caption("Tip: Use the Traversal tab to step through DFA transitions for a particular match.")
        with tab_anim:
            st.subheader("Traversal Animation — manual steps + autoplay")
            st.markdown("Select a pattern (by index) and an occurrence to animate traversal across the corresponding KMP DFA.")
            colp, colm, colc = st.columns([1,1,2])
            with colp:
                pid = st.number_input("Pattern index (0-based)", min_value=0, max_value=max(0, len(patterns)-1), value=0)
            with colm:
                occ_idx = st.number_input("Occurrence index (0-based)", min_value=0, max_value=max(0, max(len(matches_by_pattern[pid])-1,0)), value=0)
            with colc:
                step_delay = st.slider("Auto-play delay (seconds)", 0.05, 1.5, 0.35)

            # Build KMP dfa_table for the chosen pattern
            chosen_pattern = patterns[pid]
            dfa_table, m = build_kmp_automaton(chosen_pattern)
            # Prepare sequence window for that occurrence
            matches_list = matches_by_pattern[pid]
            if not matches_list:
                st.warning("No occurrences for selected pattern.")
            else:
                start_idx, end_idx = matches_list[occ_idx]
                # We'll animate the window from start_idx - context to end_idx + context (bounded)
                context = 6
                win_start = max(0, start_idx - context)
                win_end = min(len(dna)-1, end_idx + context)
                subseq = dna[win_start: win_end+1]
                st.write(f"Animating pattern '{chosen_pattern}' occurrence at [{start_idx} - {end_idx}] — showing window [{win_start}:{win_end}]")
                # Build steps: each step is (pos_in_window, dfa_state)
                steps = []
                state = 0
                for i, ch in enumerate(subseq):
                    # compute next state
                    if ch not in dfa_table:
                        ns = 0
                    else:
                        ns = dfa_table[ch][state]
                    steps.append({'pos': win_start + i, 'state': ns, 'char': ch, 'idx': i})
                    state = ns

                # place steps into session
                st.session_state.anim_steps = steps
                # reset index if different search
                if st.session_state.last_search_id != st.session_state.get('last_search_id'):
                    st.session_state.anim_index = 0

                # Controls
                c1, c2, c3, c4 = st.columns([1,1,1,1])
                if c1.button("⏮ Reset"):
                    st.session_state.anim_index = 0
                    st.session_state.anim_playing = False
                if c2.button("◀ Back"):
                    st.session_state.anim_index = max(0, st.session_state.anim_index - 1)
                    st.session_state.anim_playing = False
                if c3.button("Next ▶"):
                    st.session_state.anim_index = min(len(steps)-1, st.session_state.anim_index + 1)
                    st.session_state.anim_playing = False
                if c4.button("⏯ Play/Pause"):
                    st.session_state.anim_playing = not st.session_state.anim_playing

                # Autoplay loop (will update iteratively in current request)
                if st.session_state.anim_playing:
                    # advance index with small sleep until finished or paused
                    if st.session_state.anim_index < len(steps)-1:
                        st.session_state.anim_index += 1
                        time.sleep(step_delay)
                        # rerun to show next frame
                        st.experimental_rerun()
                    else:
                        st.session_state.anim_playing = False

                # Display subsequence with current position highlighted and pattern color highlight
                cur_idx = st.session_state.anim_index
                # Build display string: show subseq with current char boxed and pattern region highlighted
                disp_parts = []
                for i_global, ch in enumerate(subseq, start=win_start):
                    rel = i_global - win_start
                    style = ""
                    if i_global >= start_idx and i_global <= end_idx:
                        style = f"background:{PATTERN_COLORS[pid % len(PATTERN_COLORS)]};"
                    if rel == steps[cur_idx]['idx']:
                        # current char marker
                        style += "border:2px solid #ffffff;"
                    disp_parts.append(f"<span style='padding:2px 3px; margin:1px; {style}'>{ch}</span>")
                st.markdown("<div class='sequence-box'>" + "".join(disp_parts) + "</div>", unsafe_allow_html=True)

                # Show DFA visualization with current state highlighted (KMP DFA)
                # Build a simple dfa_struct for visualization: number of states m+1, transitions from dfa_table
                dfa_struct = {'num': m+1, 'trans': [], 'accepts': {m}}
                for s in range(m+1):
                    for ch in ['A','T','G','C']:
                        # for state m (accept), KMP table not defined; if s==m, transition defined by fallback: use table[ch][0] if exists
                        if s < m:
                            ns = dfa_table[ch][s]
                        else:
                            ns = dfa_table[ch][0] if m>0 else 0
                        dfa_struct['trans'].append({'from': s, 'to': ns, 'symbol': ch})
                # highlight current state
                cur_state = steps[cur_idx]['state'] if steps else 0
                dotdfa = create_dfa_graph(dfa_struct, highlight_state=cur_state)
                st.graphviz_chart(dotdfa, use_container_width=True)

        with tab_download:
            st.subheader("Export Results")
            if df.empty:
                st.warning("No results to export.")
            else:
                csv_bytes = df_to_csv_bytes(df[['Pattern','Start','End','Sequence']])
                st.download_button("Download CSV", data=csv_bytes, file_name="dna_matches.csv", mime="text/csv")

                # PDF export
                if REPORTLAB_AVAILABLE:
                    pdf_bytes = create_pdf_report("DNA Pattern Matches Report", dna, patterns, matches_by_pattern)
                    if pdf_bytes:
                        st.download_button("Download PDF Report", data=pdf_bytes, file_name="dna_matches_report.pdf", mime="application/pdf")
                else:
                    st.info("PDF export requires reportlab. Install via `pip install reportlab` to enable PDF download.")

    else:
        st.caption("Click 'Build Automaton & Search' to run detection. Use the sidebar to choose sample/patterns.")

if __name__ == "__main__":
    main()
