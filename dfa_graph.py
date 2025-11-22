def plot_dna_dfa():
    graph = """
    digraph {
        rankdir=LR
        node [shape=circle, style=filled, fontcolor=black, fontsize=16]

        q0 [fillcolor="lightgray", label="q0"]
        q1 [fillcolor="lightblue", label="q1"]
        q2 [fillcolor="orange", label="q2"]
        q3 [fillcolor="lightgreen", label="q3 (Accept)"]

        q0 -> q1 [label="A"]
        q1 -> q2 [label="T"]
        q2 -> q3 [label="G"]
        q3 -> q1 [label="A"]

        # Other transitions
        q0 -> q0 [label="T,G,C"]
        q1 -> q0 [label="A,C,G"]
        q2 -> q0 [label="A,T,C"]
        q3 -> q0 [label="T,G,C"]
    }
    """
    return graph
