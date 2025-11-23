# app.py - Main Streamlit Application (COMPLETE & FIXED v4)

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

# Sample DNA sequences
SAMPLE_SEQUENCES = {
    'Simple Test': 'ATGATCGATCGTAGATGCTAGCTGATCGATGTAAATAGCTGATCG',
    'BRCA1 Fragment': 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTG',
    'CAG Repeat Test': 'ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCG',
    'E. coli lacZ': 'ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAATAAA'
}

# ============== AUTOMATA CLASSES ==============

class State:
    def __init__(self, state_id, is_accept=False):
        self.id = state_id
        self.is_accept = is_accept
        self.transitions = defaultdict(list)
        self.epsilon_transitions = []


class NFA:
    def __init__(self):
        self.states = []
        self.start_state = None
        self.accept_states = []
        self.alphabet = ['A', 'T', 'G', 'C']

    def add_state(self, is_accept=False):
        state = State(len(self.states), is_accept)
        self.states.append(state)
        if is_accept:
            self.accept_states.append(state)
        return state


class DFA:
    def __init__(self):
        self.num_states = 0
        self.start_state = 0
        self.accept_states = []
        self.transitions = []
        self.alphabet = ['A', 'T', 'G', 'C']


# ============== NFA CONSTRUCTION ==============

def pattern_to_nfa(pattern):
    nfa = NFA()
    if not pattern:
        state = nfa.add_state(True)
        nfa.start_state = state
        return nfa
    
    prev_state = nfa.add_state()
    nfa.start_state = prev_state
    
    for i, char in enumerate(pattern):
        is_final = (i == len(pattern) - 1)
        new_state = nfa.add_state(is_final)
        prev_state.transitions[char].append(new_state)
        prev_state = new_state
    
    return nfa


def alternate_nfa(nfa1, nfa2):
    combined = NFA()
    new_start = combined.add_state()
    combined.start_state = new_start
    
    new_start.epsilon_transitions.append(nfa1.start_state)
    new_start.epsilon_transitions.append(nfa2.start_state)
    
    combined.states = [new_start] + nfa1.states + nfa2.states
    combined.accept_states = nfa1.accept_states + nfa2.accept_states
    
    for i, state in enumerate(combined.states):
        state.id = i
    
    return combined


def regex_to_nfa(regex):
    if '|' in regex:
        parts = regex.split('|')
        nfa = pattern_to_nfa(parts[0])
        for part in parts[1:]:
            nfa2 = pattern_to_nfa(part)
            nfa = alternate_nfa(nfa, nfa2)
        return nfa
    return pattern_to_nfa(regex)


# ============== NFA TO DFA CONVERSION ==============

def epsilon_closure(states):
    closure = set(states)
    stack = list(states)
    
    while stack:
        state = stack.pop()
        for eps_state in state.epsilon_transitions:
            if eps_state not in closure:
                closure.add(eps_state)
                stack.append(eps_state)
    
    return frozenset(closure)


def move(states, symbol):
    result = set()
    for state in states:
        if symbol in state.transitions:
            result.update(state.transitions[symbol])
    return result


def nfa_to_dfa(nfa):
    dfa = DFA()
    alphabet = ['A', 'T', 'G', 'C']
    
    state_map = {}
    state_counter = 0
    
    def is_accept_set(state_set):
        return any(s.is_accept for s in state_set)
    
    start_set = epsilon_closure([nfa.start_state])
    state_map[start_set] = state_counter
    state_counter += 1
    dfa.start_state = 0
    
    if is_accept_set(start_set):
        dfa.accept_states.append(0)
    
    unmarked = [start_set]
    
    while unmarked:
        current_set = unmarked.pop(0)
        current_id = state_map[current_set]
        
        for symbol in alphabet:
            move_set = move(current_set, symbol)
            if not move_set:
                continue
            
            closure_set = epsilon_closure(move_set)
            
            if closure_set not in state_map:
                state_map[closure_set] = state_counter
                state_counter += 1
                unmarked.append(closure_set)
                
                if is_accept_set(closure_set):
                    dfa.accept_states.append(state_map[closure_set])
            
            dfa.transitions.append({
                'from': current_id,
                'symbol': symbol,
                'to': state_map[closure_set]
            })
    
    dfa.num_states = state_counter
    return dfa


# ============== PATTERN MATCHING ==============

def build_kmp_automaton(pattern):
    m = len(pattern)
    alphabet = ['A', 'T', 'G', 'C']
    
    if m == 0:
        return {c: [] for c in alphabet}, 0
    
    dfa_table = {c: [0] * m for c in alphabet}
    dfa_table[pattern[0]][0] = 1
    x = 0
    
    for j in range(1, m):
        for c in alphabet:
            dfa_table[c][j] = dfa_table[c][x]
        dfa_table[pattern[j]][j] = j + 1
        x = dfa_table[pattern[j]][x]
    
    return dfa_table, m


def search_pattern(dfa_table, pattern_length, text):
    matches = []
    if pattern_length == 0:
        return matches
    
    state = 0
    for i, char in enumerate(text):
        if char not in dfa_table:
            state = 0
            continue
        
        if state < pattern_length:
            state = dfa_table[char][state]
        else:
            state = dfa_table[char][0]
        
        if state == pattern_length:
            start_pos = i - pattern_length + 1
            end_pos = i
            matches.append({
                'start': start_pos,
                'end': end_pos,
                'sequence': text[start_pos:i + 1]
            })
            state = dfa_table[char][0] if pattern_length > 0 else 0
    
    return matches


def search_multiple_patterns(patterns, text):
    all_matches = []
    seen = set()
    
    for pattern in patterns:
        if not pattern:
            continue
        dfa_table, length = build_kmp_automaton(pattern)
        matches = search_pattern(dfa_table, length, text)
        
        for match in matches:
            key = (match['start'], match['end'], match['sequence'])
            if key not in seen:
                seen.add(key)
                all_matches.append(match)
    
    all_matches.sort(key=lambda x: (x['start'], x['end']))
    return all_matches


# ============== VISUALIZATION (COMPLETELY REWRITTEN) ==============

def create_kmp_dfa_visualization(pattern):
    """Create complete KMP DFA with ALL transitions"""
    if not pattern:
        return None
    
    m = len(pattern)
    alphabet = ['A', 'T', 'G', 'C']
    
    # Build KMP DFA table
    dfa_table, _ = build_kmp_automaton(pattern)
    
    dot = graphviz.Digraph(comment='KMP DFA', engine='dot')
    dot.attr(rankdir='LR', bgcolor='#0d1117', splines='curved')
    dot.attr('node', style='filled', fontcolor='white', fontname='Arial Bold', fontsize='14')
    dot.attr('edge', fontname='Arial', fontsize='11')
    
    # Create all states (0 to m, where m is accept state)
    for i in range(m + 1):
        if i == m:
            dot.node(f'q{i}', f'q{i}', shape='doublecircle',
                     fillcolor='#22c55e', color='#16a34a', penwidth='3')
        elif i == 0:
            dot.node(f'q{i}', f'q{i}', shape='circle',
                     fillcolor='#3b82f6', color='#ef4444', penwidth='3')
        else:
            dot.node(f'q{i}', f'q{i}', shape='circle',
                     fillcolor='#3b82f6', color='#1e40af', penwidth='2')
    
    # Start arrow
    dot.node('', '', shape='none', width='0')
    dot.edge('', 'q0', color='#ef4444', penwidth='2')
    
    # Collect ALL transitions
    edge_map = defaultdict(list)
    
    for state in range(m):
        for char in alphabet:
            next_state = dfa_table[char][state]
            edge_map[(state, next_state)].append(char)
    
    # Draw all transitions with proper styling
    for (from_state, to_state), chars in edge_map.items():
        label = ','.join(sorted(chars))
        
        # Forward transition (progress in pattern)
        if to_state == from_state + 1:
            dot.edge(f'q{from_state}', f'q{to_state}', label=f' {label} ',
                    color='#4ade80', fontcolor='#4ade80', penwidth='2.5')
        
        # Backward transition (failure/mismatch - go back)
        elif to_state < from_state:
            dot.edge(f'q{from_state}', f'q{to_state}', label=f' {label} ',
                    color='#f59e0b', fontcolor='#fbbf24', style='dashed', penwidth='1.5',
                    constraint='false')
        
        # Self loop
        elif to_state == from_state:
            dot.edge(f'q{from_state}', f'q{to_state}', label=f' {label} ',
                    color='#94a3b8', fontcolor='#cbd5e1', dir='back', tailport='n', headport='n')
        
        # Jump forward (rare)
        else:
            dot.edge(f'q{from_state}', f'q{to_state}', label=f' {label} ',
                    color='#8b5cf6', fontcolor='#a78bfa', penwidth='1.5')
    
    return dot


def create_automaton_graph(dfa, pattern):
    """Create DFA visualization for regex patterns with alternation"""
    dot = graphviz.Digraph(comment='DFA', engine='dot')
    dot.attr(rankdir='LR', bgcolor='#1a1a2e', splines='curved')
    dot.attr('node', style='filled', fontcolor='white', fontname='Arial', fontsize='13')
    dot.attr('edge', fontcolor='#e74c3c', fontname='Arial', fontsize='11')
    
    # Create states
    for i in range(dfa.num_states):
        if i in dfa.accept_states:
            dot.node(f'q{i}', f'q{i}', shape='doublecircle',
                     fillcolor='#2ecc71', color='#27ae60', penwidth='2')
        elif i == dfa.start_state:
            dot.node(f'q{i}', f'q{i}', shape='circle',
                     fillcolor='#3498db', color='#e74c3c', penwidth='3')
        else:
            dot.node(f'q{i}', f'q{i}', shape='circle',
                     fillcolor='#3498db', color='#2c3e50', penwidth='2')
    
    # Start arrow
    dot.node('', '', shape='none', width='0')
    dot.edge('', f'q{dfa.start_state}', color='#e74c3c', penwidth='2')
    
    # Group transitions by edge
    edge_labels = defaultdict(list)
    for trans in dfa.transitions:
        key = (trans['from'], trans['to'])
        edge_labels[key].append(trans['symbol'])
    
    # Draw edges
    for (from_state, to_state), symbols in edge_labels.items():
        label = ','.join(sorted(symbols))
        if from_state == to_state:
            dot.edge(f'q{from_state}', f'q{to_state}', label=label,
                     color='#95a5a6', fontcolor='#ecf0f1')
        else:
            dot.edge(f'q{from_state}', f'q{to_state}', label=label,
                     color='#7f8c8d', fontcolor='#e74c3c', penwidth='1.5')
    
    return dot


# ============== HELPER FUNCTIONS ==============

def highlight_matches(dna, matches):
    if not matches:
        return f'<div class="sequence-box">{dna}</div>'
    
    html_parts = []
    last_end = -1
    
    sorted_matches = sorted(matches, key=lambda x: x['start'])
    
    for match in sorted_matches:
        start = match['start']
        end = match['end']
        
        if start > last_end + 1:
            html_parts.append(dna[last_end + 1:start])
        elif start <= last_end:
            continue
        
        html_parts.append(f'<span class="match-highlight">{dna[start:end + 1]}</span>')
        last_end = end
    
    if last_end < len(dna) - 1:
        html_parts.append(dna[last_end + 1:])
    
    result = ''.join(html_parts)
    return f'<div class="sequence-box">{result}</div>'


def get_pattern_info(pattern):
    pattern_db = {
        'ATG': {
            'name': 'Start Codon (Methionine)',
            'description': 'Signals the start of protein translation. Every protein-coding sequence begins with ATG.',
            'significance': 'Essential for gene expression. Mutations here prevent protein synthesis.',
            'category': 'Codon'
        },
        'TAA': {
            'name': 'Stop Codon (Ochre)',
            'description': 'One of three stop codons that signal the termination of protein translation.',
            'significance': 'Mutations creating premature stop codons cause truncated proteins.',
            'category': 'Codon'
        },
        'TAG': {
            'name': 'Stop Codon (Amber)',
            'description': 'Stop codon that terminates protein synthesis.',
            'significance': 'Used in genetic engineering for incorporating non-natural amino acids.',
            'category': 'Codon'
        },
        'TGA': {
            'name': 'Stop Codon (Opal)',
            'description': 'Stop codon, also known as Opal. Can code for Selenocysteine in special contexts.',
            'significance': 'Sometimes read through to incorporate rare amino acids.',
            'category': 'Codon'
        },
        'TATA': {
            'name': 'TATA Box',
            'description': 'DNA sequence found in promoter region approximately 25-30 base pairs upstream of transcription start.',
            'significance': 'Critical for transcription initiation. TATA-binding protein binds here.',
            'category': 'Promoter Element'
        },
        'GAATTC': {
            'name': 'EcoRI Restriction Site',
            'description': 'Recognition sequence for EcoRI restriction endonuclease from E. coli.',
            'significance': 'Widely used in molecular cloning and genetic engineering.',
            'category': 'Restriction Site'
        },
        'CAG': {
            'name': 'CAG Trinucleotide (Glutamine)',
            'description': 'Codes for the amino acid Glutamine. Abnormal expansion causes several diseases.',
            'significance': 'More than 36 CAG repeats in HTT gene causes Huntingtons disease.',
            'category': 'Trinucleotide Repeat'
        },
        'AATAAA': {
            'name': 'Polyadenylation Signal',
            'description': 'Consensus sequence that signals the addition of poly-A tail to mRNA.',
            'significance': 'Essential for mRNA stability and export from nucleus.',
            'category': 'Regulatory Element'
        }
    }
    
    return pattern_db.get(pattern.upper(), {
        'name': 'Custom Pattern',
        'description': 'User-defined search pattern.',
        'significance': 'Searching for exact sequence matches.',
        'category': 'Custom'
    })


# ============== MAIN APP (COMPLETELY FIXED) ==============

def main():
    # Sidebar with IMMEDIATE updates
    with st.sidebar:
        st.header("⚙️ Settings")
        
        st.subheader("📁 Load Sample Sequence")
        sample_choice = st.selectbox(
            "Choose sample:", 
            list(SAMPLE_SEQUENCES.keys()),
            key='sample_selector_unique'
        )
        
        st.subheader("🔬 Quick Patterns")
        pattern_choice = st.radio(
            "Select pattern:",
            ['Custom', 'ATG', 'TAA|TAG|TGA', 'TATA', 'GAATTC', 'CAG', 'AATAAA'],
            key='pattern_selector_unique'
        )
        
        st.markdown('---')
        st.subheader("📚 About")
        st.info("**TOC Concepts:**\n- Regular Expressions\n- NFA Construction\n- NFA to DFA\n- KMP Automata")
    
    # Header
    st.markdown('<h1 class="main-header">🧬 DNA Pattern Matcher</h1>', unsafe_allow_html=True)
    st.markdown('<p style="text-align: center; color: #888;">Using Finite Automata | Theory of Computation Project</p>', unsafe_allow_html=True)
    st.markdown('---')
    
    # Main content with values FROM sidebar
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("🧬 DNA Sequence Input")
        dna_input = st.text_area(
            "Enter DNA sequence (A, T, G, C only):",
            value=SAMPLE_SEQUENCES[sample_choice],  # Directly use sidebar value
            height=150,
            key="dna_textarea"
        )
    
    with col2:
        st.subheader("🔍 Pattern")
        # Use pattern from sidebar radio button
        if pattern_choice == 'Custom':
            pattern_input = st.text_input(
                "Enter pattern:",
                value='ATG',
                key="pattern_field"
            )
        else:
            pattern_input = st.text_input(
                "Enter pattern:",
                value=pattern_choice,
                key="pattern_field"
            )
        st.caption("**Examples:** ATG, TATA, TAA|TAG|TGA")
    
    # Search button
    search_clicked = st.button("🔍 Build Automaton & Search", type="primary", use_container_width=True)
    
    if search_clicked:
        if not dna_input or not pattern_input:
            st.error("❌ Please enter both DNA sequence and pattern.")
            return
        
        dna_clean = ''.join(c for c in dna_input.upper() if c in 'ATGC')
        pattern_clean = ''.join(c for c in pattern_input.upper() if c in 'ATGC|')
        
        if not dna_clean:
            st.error("❌ Invalid DNA sequence.")
            return
        
        if not pattern_clean or pattern_clean == '|':
            st.error("❌ Invalid pattern.")
            return
        
        if '|' in pattern_clean:
            patterns = [p for p in pattern_clean.split('|') if p]
            if not patterns:
                st.error("❌ Invalid pattern.")
                return
            pattern_clean = '|'.join(patterns)
        
        st.markdown('---')
        
        # Tabs
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Results", "🔄 Automaton", "📋 Details", "📖 Info"])
        
        with tab1:
            if '|' in pattern_clean:
                patterns = pattern_clean.split('|')
                matches = search_multiple_patterns(patterns, dna_clean)
            else:
                dfa_table, pattern_length = build_kmp_automaton(pattern_clean)
                matches = search_pattern(dfa_table, pattern_length, dna_clean)
            
            col_stat1, col_stat2, col_stat3 = st.columns(3)
            with col_stat1:
                st.metric("Sequence Length", f"{len(dna_clean)} bp")
            with col_stat2:
                st.metric("Pattern", pattern_clean)
            with col_stat3:
                st.metric("Matches Found", len(matches))
            
            st.subheader("Highlighted Matches")
            st.markdown(highlight_matches(dna_clean, matches), unsafe_allow_html=True)
        
        with tab2:
            st.subheader("Finite Automaton (DFA)")
            
            if '|' in pattern_clean:
                nfa = regex_to_nfa(pattern_clean)
                dfa = nfa_to_dfa(nfa)
                st.caption("DFA for multiple patterns (alternation)")
                dot = create_automaton_graph(dfa, pattern_clean)
            else:
                st.caption("KMP-style DFA")
                st.info("🟢 **Green** = Forward | 🟠 **Orange dashed** = Failure/back | ⚪ **Gray** = Self-loop")
                dot = create_kmp_dfa_visualization(pattern_clean)
                nfa = pattern_to_nfa(pattern_clean)
                dfa = nfa_to_dfa(nfa)
            
            if dot:
                st.graphviz_chart(dot, use_container_width=True)
            
            col_dfa1, col_dfa2, col_dfa3 = st.columns(3)
            with col_dfa1:
                st.metric("States", dfa.num_states)
            with col_dfa2:
                st.metric("Transitions", len(dfa.transitions))
            with col_dfa3:
                st.metric("Accept States", len(dfa.accept_states))
        
        with tab3:
            st.subheader("Match Details")
            
            if matches:
                import pandas as pd
                df = pd.DataFrame(matches)
                df.index = df.index + 1
                df.columns = ['Start', 'End', 'Sequence']
                st.dataframe(df, use_container_width=True)
                
                csv = df.to_csv(index=True)
                st.download_button(
                    "📥 Download CSV",
                    data=csv,
                    file_name="dna_matches.csv",
                    mime="text/csv"
                )
            else:
                st.warning("No matches found.")
        
        with tab4:
            display_pattern = pattern_clean.split('|')[0] if '|' in pattern_clean else pattern_clean
            info = get_pattern_info(display_pattern)
            
            st.subheader(f"📖 {info['name']}")
            st.markdown(f"""
<div class="info-box">
<p><strong>Category:</strong> {info['category']}</p>
<p><strong>Description:</strong> {info['description']}</p>
<p><strong>Significance:</strong> {info['significance']}</p>
</div>
""", unsafe_allow_html=True)
            
            if display_pattern == 'CAG':
                import re
                st.markdown('---')
                st.subheader("🔬 CAG Repeat Analysis")
                
                cag_repeats = re.findall(r'((?:CAG)+)', dna_clean)
                if cag_repeats:
                    longest = max(cag_repeats, key=len)
                    repeat_count = len(longest) // 3
                    
                    st.metric("Longest Consecutive Repeat", f"{repeat_count} CAGs")
                    
                    if repeat_count > 36:
                        st.error(f"⚠️ {repeat_count} repeats - Huntingtons disease range!")
                    elif repeat_count > 26:
                        st.warning(f"⚡ {repeat_count} repeats - Intermediate range")
                    else:
                        st.success(f"✅ {repeat_count} repeats - Normal range")
    
    # Footer
    st.markdown('---')
    st.markdown(
        '<p style="text-align: center; color: #666; font-size: 0.8rem;">'
        '🧬 DNA Pattern Matcher | TOC Mini Project'
        '</p>',
        unsafe_allow_html=True
    )


if __name__ == "__main__":
    main()
