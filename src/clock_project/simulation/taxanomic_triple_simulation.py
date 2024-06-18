from clock_project.simulation.wts import SeqSimulate


def taxanomic_triple_simulation(p0, Q1, Q2, Q3, t1, t2, length, num_repeat, seed):
    simulator1 = SeqSimulate(Q1, length, num_repeat, seed, p0)
    seqs_inter_node = simulator1.main(max_time=t1)[0]
    seq_inter_node = seqs_inter_node[-1]
    simulator2 = SeqSimulate(Q2, length, num_repeat, seed, p0, seq_inter_node)
    simulator3 = SeqSimulate(Q3, length, num_repeat, seed, p0, seq_inter_node)
    internal_t = t2-t1
    seqs_edge_1 = simulator2.main(max_time=internal_t)
    seqs_edge_2 = simulator3.main(max_time=internal_t)
    seqs_edge_3 = simulator1.main(max_time=t2)[0]
    seqs = {'in_group': {'edge1': seqs_edge_1[-1], 'edge2': seqs_edge_2[-1], 'out_group': seqs_edge_3[-1]}}
    return seqs
