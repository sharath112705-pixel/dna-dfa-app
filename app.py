# app.py - Main Streamlit Application

import streamlit as st
import graphviz
from collections import defaultdict
import json

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
    .stTabs [data-baseweb="tab-list"] {
        gap: 24px;
    }
    .stTabs [data-baseweb="tab"] {
        background-color: #1a1a2e;
        border-radius: 5px;
        padding: 10px 20px;
    }
</style>
""", unsafe_allow_html=True)

# ============== AUTOMATA CLASSES ==============

class State:
    """Represents a state in NFA/DFA"""
    def __init__(self, id, is_accept=False):
        self.id = id
        self.is_accept = is_accept
        self.transitions = defaultdict(list)  # symbol -> [states]
        self.epsilon_transitions = []  # ε-transitions

class NFA:
    """Non-deterministic Finite Automaton"""
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
    """Deterministic Finite Automaton"""
    def __init__(self):
        self.num_states = 0
        self.start_state = 0
        self.accept_states = []
        self.transitions = []  # list of {from, symbol, to}
        self.alphabet = ['A', 'T', 'G', 'C']

# ============== NFA CONSTRUCTION ==============

def pattern_to_nfa(pattern):
    """Convert simple pattern to NFA (Thompson's construction for concatenation)"""
    nfa = NFA()
    
    # Create start state
    prev_state = nfa.add_state()
    nfa.start_state = prev_state
    
    # Create chain of states for each character
    for i, char in enumerate(pattern):
        is_final = (i == len(pattern) - 1)
        new_state = nfa.add_state(is_final)
        prev_state.transitions[char].append(new_state)
        prev_state = new_state
    
    return nfa

def concatenate_nfa(nfa1, nfa2):
    """Concatenate two NFAs"""
    for accept_state in nfa1.accept_states:
        accept_state.is_accept = False
        accept_state.epsilon_transitions.append(nfa2.start_state)
    
    nfa1.accept_states = nfa2.accept_states
    nfa1.states.extend(nfa2.states)
    
    # Renumber states
    for i, state in enumerate(nfa1.states):
        state.id = i
    
    return nfa1

def alternate_nfa(nfa1, nfa2):
    """Create union of two NFAs (alternation)"""
    combined = NFA()
    new_start = combined.add_state()
    combined.start_state = new_start
    
    new_start.epsilon_transitions.append(nfa1.start_state)
    new_start.epsilon_transitions.append(nfa2.start_state)
    
    combined.states = [new_start] + nfa1.states + nfa2.states
    combined.accept_states = nfa1.accept_states + nfa2.accept_states
    
    # Renumber states
    for i, state in enumerate(combined.states):
        state.id = i
    
    return combined

def regex_to_nfa(regex):
    """Convert simple regex (with | for alternation) to NFA"""
    if '|' in regex:
        parts = regex.split('|')
        nfa = pattern_to_nfa(parts[0])
        for part in parts[1:]:
            nfa2 = pattern_to_nfa(part)
            nfa = alternate_nfa(nfa, nfa2)
        return nfa
    else:
        return pattern_to_nfa(regex)

# ============== NFA TO DFA CONVERSION ==============

def epsilon_closure(states):
    """Compute epsilon closure of a set of states"""
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
    """Compute move(states, symbol)"""
    result = set()
    for state in states:
        if symbol in state.transitions:
            result.update(state.transitions[symbol])
    return result

def nfa_to_dfa(nfa):
    """Convert NFA to DFA using subset construction"""
    dfa = DFA()
    alphabet = ['A', 'T', 'G', 'C']
    
    # Map from frozenset of NFA states to DFA state ID
    state_map = {}
    state_counter = 0
    
    def is_accept_set(state_set):
        return any(s.is_accept for s in state_set)
    
    # Start state is epsilon-closure of NFA start
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
    """Build KMP-style DFA for efficient pattern matching"""
    m = len(pattern)
    alphabet = ['A', 'T', 'G', 'C']
    dfa_table = {c: [0] * m for c in alphabet}
    
    if m == 0:
        return dfa_table, m
    
    dfa_table[pattern[0]][0] = 1
    x = 0  # restart state
    
    for j in range(1, m):
        for c in alphabet:
            dfa_table[c][j] = dfa_table[c][x]
        dfa_table[pattern[j]][j] = j + 1
        x = dfa_table[pattern[j]][x]
    
    return dfa_table, m

def search_pattern(dfa_table, pattern_length, text):
    """Search for pattern in text using DFA"""
    matches = []
    state = 0
    
    for i, char in enumerate(text):
        if char not in dfa_table:
            state = 0
            continue
        
        state = dfa_table[char][state]
        
        if state == pattern_length:
            matches.append({
                'start': i - pattern_length + 1,
                'end': i + 1,
                'sequence': text[i - pattern_length + 1:i + 1]
            })
            # Allow overlapping matches
            if pattern_length > 1:
                state = dfa_table[char][state - 1] if state > 0 else 0
            else:
                state = 0
    
    return matches

def search_multiple_patterns(patterns, text):
    """Search for multiple patterns (for regex with |)"""
    all_matches = []
    for pattern in patterns:
        dfa_table, length = build_kmp_automaton(pattern)
        matches = search_pattern(dfa_table, length, text)
        all_matches.extend(matches)
    
    # Sort by position and remove duplicates
    all_matches.sort(key=lambda x: (x['start'], x['end']))
    return all_matches

# ============== VISUALIZATION ==============

def create_automaton_graph(dfa, pattern):
    """Create Graphviz visualization of DFA"""
    dot = graphviz.Digraph(comment='DFA')
    dot.attr(rankdir='LR', bgcolor='#1a1a2e')
    dot.attr('node', style='filled', fontcolor='white', fontname='Arial')
    dot.attr('edge', fontcolor='#e74c3c', fontname='Arial', fontsize='12')
    
    # Add states
    for i in range(dfa.num_states):
        if i in dfa.accept_states:
            dot.node(f'q{i}', f'q{i}', shape='doublecircle', 
                    fillcolor='#2ecc71', color='#27ae60')
        elif i == dfa.start_state:
            dot.node(f'q{i}', f'q{i}', shape='circle', 
                    fillcolor='#3498db', color='#e74c3c', penwidth='3')
        else:
            dot.node(f'q{i}', f'q{i}', shape='circle', 
                    fillcolor='#3498db', color='#2c3e50')
    
    # Add start arrow
    dot.node('start', '', shape='none', width='0', height='0')
    dot.edge('start', f'q{dfa.start_state}', color='#e74c3c', penwidth='2')
    
    # Group transitions by (from, to) to combine labels
    edge_labels = defaultdict(list)
    for trans in dfa.transitions:
        key = (trans['from'], trans['to'])
        edge_labels[key].append(trans['symbol'])
    
    # Add transitions
    for (from_state, to_state), symbols in edge_labels.items():
        label = ','.join(sorted(symbols))
        if from_state == to_state:
            dot.edge(f'q{from_state}', f'q{to_state}', label=label, 
                    color='#666', fontcolor='#e74c3c')
        else:
            dot.edge(f'q{from_state}', f'q{to_state}', label=label, 
                    color='#666', fontcolor='#e74c3c')
    
    return dot

def create_simple_dfa_graph(pattern):
    """Create a simple linear DFA visualization for the pattern"""
    dot = graphviz.Digraph(comment='Pattern DFA')
    dot.attr(rankdir='LR', bgcolor='#0d1117')
    dot.attr('node', style='filled', fontcolor='white', fontname='Arial Bold', fontsize='14')
    dot.attr('edge', fontcolor='#ff6b6b', fontname='Arial Bold', fontsize='14', color='#4a5568')
    
    # Create states
    for i in range(len(pattern) + 1):
        if i == len(pattern):  # Accept state
            dot.node(f'q{i}', f'q{i}', shape='doublecircle', 
                    fillcolor='#22c55e', color='#16a34a', penwidth='2')
        elif i == 0:  # Start state
            dot.node(f'q{i}', f'q{i}', shape='circle', 
                    fillcolor='#3b82f6', color='#ef4444', penwidth='3')
        else:
            dot.node(f'q{i}', f'q{i}', shape='circle', 
                    fillcolor='#3b82f6', color='#1e40af', penwidth='2')
    
    # Start arrow
    dot.node('start', '', shape='none', width='0', height='0')
    dot.edge('start', 'q0', color='#ef4444', penwidth='2')
    
    # Main transitions
    for i, char in enumerate(pattern):
        dot.edge(f'q{i}', f'q{i+1}', label=f' {char} ', penwidth='2')
    
    return dot

# ============== HELPER FUNCTIONS ==============

def highlight_matches(dna, matches):
    """Create HTML with highlighted matches"""
    if not matches:
        return f'<div class="sequence-box">{dna}</div>'
    
    html_parts = []
    last_end = 0
    
    # Sort matches by start position
    sorted_matches = sorted(matches, key=lambda x: x['start'])
    
    for match in sorted_matches:
        # Add non-matching part
        if match['start'] > last_end:
            html_parts.append(dna[last_end:match['start']])
        
        # Add matching part with highlight
        html_parts.append(f'<span class="match-highlight">{dna[match["start"]:match["end"]]}</span>')
        last_end = match['end']
    
    # Add remaining part
    if last_end < len(dna):
        html_parts.append(dna[last_end:])
    
    # Word wrap
    result = ''.join(html_parts)
    return f'<div class="sequence-box">{result}</div>'

def get_pattern_info(pattern):
    """Get biological information about the pattern"""
    pattern_db = {
        'ATG': {
            'name': 'Start Codon (Methionine)',
            'description': 'Signals the start of protein translation. Every protein-coding sequence begins with ATG, which codes for the amino acid Methionine.',
            'significance': 'Essential for gene expression. Mutations here prevent protein synthesis.',
            'category': 'Codon'
        },
        'TAA': {
            'name': 'Stop Codon (Ochre)',
            'description': 'One of three stop codons that signal the termination of protein translation.',
            'significance': 'Mutations creating premature stop codons cause truncated, non-functional proteins.',
            'category': 'Codon'
        },
        'TAG': {
            'name': 'Stop Codon (Amber)',
            'description': 'Stop codon that terminates protein synthesis. Named "Amber" after discoverer Harris Bernstein.',
            'significance': 'Used in genetic engineering for incorporating non-natural amino acids.',
            'category': 'Codon'
        },
        'TGA': {
            'name': 'Stop Codon (Opal)',
            'description': 'Stop codon, also known as "Opal" or "Umber". Can code for Selenocysteine in special contexts.',
            'significance': 'Sometimes "read through" to incorporate rare amino acids.',
            'category': 'Codon'
        },
        'TATA': {
            'name': 'TATA Box',
            'description': 'A DNA sequence (5\'-TATAAA-3\') found in the promoter region approximately 25-30 base pairs upstream of the transcription start site.',
            'significance': 'Critical for transcription initiation. TATA-binding protein (TBP) binds here.',
            'category': 'Promoter Element'
        },
        'GAATTC': {
            'name': 'EcoRI Restriction Site',
            'description': 'Recognition sequence for EcoRI restriction endonuclease from E. coli. Cuts between G and A.',
            'significance': 'Widely used in molecular cloning and genetic engineering.',
            'category': 'Restriction Site'
        },
        'GGATCC': {
            'name': 'BamHI Restriction Site',
            'description': 'Recognition sequence for BamHI restriction enzyme from Bacillus amyloliquefaciens.',
            'significance': 'Common restriction site used in cloning vectors.',
            'category': 'Restriction Site'
        },
        'CAG': {
            'name': 'CAG Trinucleotide (Glutamine)',
            'description': 'Codes for the amino acid Glutamine. Abnormal expansion of CAG repeats causes several diseases.',
            'significance': 'More than 36 CAG repeats in HTT gene causes Huntington\'s disease.',
            'category': 'Trinucleotide Repeat'
        },
        'CGG': {
            'name': 'CGG Trinucleotide',
            'description': 'Codes for Arginine. Expansion in FMR1 gene causes Fragile X syndrome.',
            'significance': 'More than 200 repeats causes intellectual disability.',
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

# Sample DNA sequences
SAMPLE_SEQUENCES = {
    'Simple Test': 'ATGATCGATCGTAGATGCTAGCTGATCGATGTAAATAGCTGATCG',
    'BRCA1 Fragment': 'ATGGATTTATCTGCTCTTCGCGTTGAAGAAGTACAAAATGTCATTAATGCTATGCAGAAAATCTTAGAGTGTCCCATCTGTCTGGAGTTGATCAAGGAACCTGTCTCCACAAAGTGTGACCACATATTTTGCAAATTTTGCATGCTGAAACTTCTCAACCAGAAGAAAGGGCCTTCACAGTGTCCTTTATGTAAGAATGATATAACCAAAAGGAGCCTACAAGAAAGTACGAGATTTAGTCAACTTGTTGAAGAGCTATTGAAAATCATTTGTGCTTTTCAGCTTGACACAGGTTTGGAGTATGCAAACAGCTATAATTTTGCAAAAAAGGAAAATAACTCTCCTGAACATCTAAAAGATGAAGTTTCTATCATCCAAAGTATGGGCTACAGAAACCGTGCCAAAAGACTTCTACAGAGTGAACCCGAAAATCCTTCCTTG',
    'CAG Repeat Test': 'ATGAAGGCCTTCGAGTCCCTCAAGTCCTTCCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAACAGCCGCCACCGCCGCCGCCGCCGCCGCCTCCTCAGCTTCCTCAGCCGCCGCCG',
    'E. coli lacZ': 'ATGACCATGATTACGGATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCTTTGCCTGGTTTCCGGCACCAGAAGCGGTGCCGGAAAGCTGGCTGGAGTGCGATCTTCCTGAGGCCGATACTGTCGTCGTCCCCTCAAACTGGCAGATGCACGGTTACGATGCGCCCATCTACACCAACGTGACCTATCCCATTACGGTCAATCCGCCGTTTGTTCCCACGGAGAATCCGACGGGTTGTTACTCGCTCACATTTAATGTTGATGAAAGCTGGCTACAGGAAGGCCAGACGCGAATTATTTTTGATGGCGTTAATAAA'
}

# ============== MAIN APP ==============

def main():
    # Header
    st.markdown('<h1 class="main-header">🧬 DNA Pattern Matcher</h1>', unsafe_allow_html=True)
    st.markdown('<p style="text-align: center; color: #888;">Using Finite Automata for Genetic Sequence Analysis | Theory of Computation Project</p>', unsafe_allow_html=True)
    st.markdown('---')
    
    # Sidebar
    with st.sidebar:
        st.header("⚙️ Settings")
        
        st.subheader("📁 Load Sample Sequence")
        sample_choice = st.selectbox("Choose sample:", list(SAMPLE_SEQUENCES.keys()))
        
        st.subheader("🔬 Quick Patterns")
        pattern_buttons = {
            'Start Codon': 'ATG',
            'Stop Codons': 'TAA|TAG|TGA',
            'TATA Box': 'TATA',
            'EcoRI Site': 'GAATTC',
            'CAG Repeat': 'CAG',
            'Poly-A Signal': 'AATAAA'
        }
        
        selected_pattern = None
        for name, pattern in pattern_buttons.items():
            if st.button(name, key=f"btn_{pattern}", use_container_width=True):
                selected_pattern = pattern
        
        st.markdown('---')
        st.subheader("📚 About")
        st.info("""
        **TOC Concepts Used:**
        - Regular Expressions
        - NFA Construction
        - NFA → DFA Conversion
        - String Matching Automata
        """)
    
    # Main content
    col1, col2 = st.columns([2, 1])
    
    with col1:
        st.subheader("🧬 DNA Sequence Input")
        
        # Initialize session state
        if 'dna_sequence' not in st.session_state:
            st.session_state.dna_sequence = SAMPLE_SEQUENCES['Simple Test']
        if 'pattern' not in st.session_state:
            st.session_state.pattern = 'ATG'
        
        # Update from sidebar selections
        if sample_choice:
            dna_input = st.text_area(
                "Enter DNA sequence (A, T, G, C only):",
                value=SAMPLE_SEQUENCES[sample_choice],
                height=150,
                key="dna_input"
            )
        else:
            dna_input = st.text_area(
                "Enter DNA sequence (A, T, G, C only):",
                value=st.session_state.dna_sequence,
                height=150,
                key="dna_input"
            )
    
    with col2:
        st.subheader("🔍 Pattern")
        pattern_input = st.text_input(
            "Enter pattern to search:",
            value=selected_pattern if selected_pattern else st.session_state.pattern,
            key="pattern_input",
            help="Use | for multiple patterns (e.g., TAA|TAG|TGA)"
        )
        
        st.caption("**Examples:** ATG, TATA, TAA|TAG|TGA")
    
    # Search button
    search_clicked = st.button("🔍 Build Automaton & Search", type="primary", use_container_width=True)
    
    if search_clicked and dna_input and pattern_input:
        # Clean inputs
        dna_clean = dna_input.upper().replace('\n', '').replace(' ', '')
        dna_clean = ''.join(c for c in dna_clean if c in 'ATGC')
        pattern_clean = pattern_input.upper().replace(' ', '')
        
        # Validate
        if not dna_clean:
            st.error("❌ Invalid DNA sequence. Please use only A, T, G, C characters.")
            return
        
        if not all(c in 'ATGC|' for c in pattern_clean):
            st.error("❌ Invalid pattern. Please use only A, T, G, C characters (and | for alternation).")
            return
        
        st.markdown('---')
        
        # Create tabs for results
        tab1, tab2, tab3, tab4 = st.tabs(["📊 Results", "🔄 Automaton", "📋 Match Details", "📖 Pattern Info"])
        
        with tab1:
            # Perform search
            if '|' in pattern_clean:
                patterns = pattern_clean.split('|')
                matches = search_multiple_patterns(patterns, dna_clean)
            else:
                dfa_table, pattern_length = build_kmp_automaton(pattern_clean)
                matches = search_pattern(dfa_table, pattern_length, dna_clean)
            
            # Display stats
            col_stat1, col_stat2, col_stat3 = st.columns(3)
            with col_stat1:
                st.metric("Sequence Length", f"{len(dna_clean)} bp")
            with col_stat2:
                st.metric("Pattern", pattern_clean)
            with col_stat3:
                st.metric("Matches Found", len(matches))
            
            # Display highlighted sequence
            st.subheader("Sequence with Highlighted Matches")
            st.markdown(highlight_matches(dna_clean, matches), unsafe_allow_html=True)
        
        with tab2:
            st.subheader("Finite Automaton (DFA)")
            
            # Build and visualize DFA
            if '|' in pattern_clean:
                nfa = regex_to_nfa(pattern_clean)
                dfa = nfa_to_dfa(nfa)
                st.caption("DFA constructed from NFA using subset construction")
                dot = create_automaton_graph(dfa, pattern_clean)
            else:
                st.caption("KMP-style DFA for pattern matching")
                dot = create_simple_dfa_graph(pattern_clean)
                nfa = pattern_to_nfa(pattern_clean)
                dfa = nfa_to_dfa(nfa)
            
            st.graphviz_chart(dot, use_container_width=True)
            
            # DFA stats
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
                # Create dataframe
                import pandas as pd
                df = pd.DataFrame(matches)
                df.index = df.index + 1
                df.columns = ['Start Position', 'End Position', 'Matched Sequence']
                st.dataframe(df, use_container_width=True)
                
                # Download button
                csv = df.to_csv(index=True)
                st.download_button(
                    label="📥 Download Results as CSV",
                    data=csv,
                    file_name="dna_matches.csv",
                    mime="text/csv"
                )
            else:
                st.warning("No matches found for the given pattern.")
        
        with tab4:
            # Get pattern info (for first pattern if multiple)
            display_pattern = pattern_clean.split('|')[0] if '|' in pattern_clean else pattern_clean
            info = get_pattern_info(display_pattern)
            
            st.subheader(f"📖 {info['name']}")
            
            st.markdown(f"""
            <div class="info-box">
                <p><strong>Category:</strong> {info['category']}</p>
                <p><strong>Description:</strong> {info['description']}</p>
                <p><strong>Biological Significance:</strong> {info['significance']}</p>
            </div>
            """, unsafe_allow_html=True)
            
            if display_pattern == 'CAG':
                # Special CAG repeat analysis
                cag_count = dna_clean.count('CAG')
                st.markdown('---')
                st.subheader("🔬 CAG Repeat Analysis")
                
                # Find longest consecutive CAG repeat
                import re
                cag_repeats = re.findall(r'((?:CAG)+)', dna_clean)
                if cag_repeats:
                    longest = max(cag_repeats, key=len)
                    repeat_count = len(longest) // 3
                    
                    st.metric("Longest Consecutive CAG Repeat", f"{repeat_count} repeats")
                    
                    if repeat_count > 36:
                        st.error(f"⚠️ WARNING: {repeat_count} CAG repeats detected. Values >36 are associated with Huntington's disease.")
                    elif repeat_count > 26:
                        st.warning(f"⚡ CAUTION: {repeat_count} CAG repeats is in the intermediate range (27-35).")
                    else:
                        st.success(f"✅ {repeat_count} CAG repeats is within normal range (<26).")

    # Footer
    st.markdown('---')
    st.markdown("""
    <p style="text-align: center; color: #666; font-size: 0.8rem;">
        🧬 DNA Pattern Matcher | Theory of Computation Mini Project<br>
        Built with Streamlit | Concepts: NFA, DFA, Thompson's Construction, Subset Construction
    </p>
