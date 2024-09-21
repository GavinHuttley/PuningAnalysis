import multiprocessing
import numpy as np
import json
import os
from cogent3 import get_app, open_data_store
from cogent3.maths.measure import jsd
from collections import Counter
from clock_project.genome_analysis.analytical_correlation_analysis import get_nabla_diff_log_ratio
import click

from clock_project.genome_analysis.analytical_correlation_analysis import (
    get_nabla_diff_log_ratio,
)

# Initialize the load_json_app outside of the function to avoid re-initialization costs
load_json_app = get_app("load_json")
matrix_collection_dir = '/Users/gulugulu/Desktop/honours/data_local_2/valid_matrix_full.json'
ancester_seq_dir = '/Users/gulugulu/Desktop/honours/data_local_2/seed_ancester_sequences.json'
simulated_info_dir = '/Users/gulugulu/Desktop/honours/data_local_2/output_simulated_sequence'

# Load matrices and ancestor sequences once
with open(matrix_collection_dir, 'r') as infile:
    matrices_dict = json.load(infile)
matrix_collection = [item for sublist in matrices_dict.values() for item in sublist]

ancester_seqs = json.load(open(ancester_seq_dir, 'r'))

def calculate_percentages(ancestral_seq):
    count = Counter(ancestral_seq)
    total_count = len(ancestral_seq)
    return [(count.get(str(i), 0) / total_count) for i in range(4)]  # Assuming 0 to 3 are the keys

def change_nucleotide_order(motif_prob):
    nuc_freq = {'A': motif_prob[0], 'C': motif_prob[1], 'G': motif_prob[2], 'T': motif_prob[3]}
    return [nuc_freq[nuc] for nuc in ['T', 'C', 'A', 'G']]

def process_alignment(args):
    path, aln_ens, initial_pi, t = args
    identifier = path.unique_id.split('.')[0]
    matrices = matrix_collection[int(identifier)]
    Q1 = np.array(matrices['ingroup1'])
    Q2 = np.array(matrices['ingroup2'])
    aln = load_json_app(path)
    try:
        motif_probs = aln.probs_per_seq()
        nuc_freq1 = change_nucleotide_order(list(motif_probs['ingroup_edge1']))
        nuc_freq2 = change_nucleotide_order(list(motif_probs['ingroup_edge2']))
        ens1 = aln_ens[identifier]['ingroup_edge1']
        ens2 = aln_ens[identifier]['ingroup_edge2']
        ens_abs_diff = abs(ens1 - ens2)
        jsd_value = np.sqrt(jsd(nuc_freq1, nuc_freq2))
        jsd_diff = abs(jsd(nuc_freq1, initial_pi) - jsd(nuc_freq2, initial_pi))
        ens_diff_log_ratio = abs(np.log(ens1 / ens2))
        nabla_diff_log_ratio = get_nabla_diff_log_ratio(initial_pi, Q1, Q2, t)
        return (identifier, t, ens_diff_log_ratio, ens_abs_diff, nabla_diff_log_ratio, jsd_value, jsd_diff)
    except Exception as e:
        print(f"Error with {identifier}: {e}")
        return None
    
@click.command()
@click.option('--comb', default='low_short', help='Choose the computation scenario.')


def main(comb):
    pool = multiprocessing.Pool(processes=4)
    tasks = []
    results = []
    pi_info, length_info = comb.split('_')
    length = 500 if length_info == 'short' else 5000
    initial_pi = calculate_percentages(ancester_seqs[comb])

    for t in [0.5, 1.0, 1.5, 2.0]:
        aln_path_dir = os.path.join(simulated_info_dir, f'aln_{pi_info}_{length}_{t}')
        print(f'aln_{pi_info}_{length}_{t}')
        aln_paths = open_data_store(aln_path_dir, suffix='json')
        print(len(aln_paths))
        aln_ens_path = os.path.join(simulated_info_dir, f'ens_results_{pi_info}_{t}_{length}.json')
        aln_ens = json.load(open(aln_ens_path, 'r'))

        for path in aln_paths:
            tasks.append((path, aln_ens, initial_pi, t))

    results = pool.map(process_alignment, tasks)
    pool.close()
    pool.join()

    with open(f'/Users/gulugulu/Desktop/honours/data_local_2/simulated_aln_info_{comb}.json', 'w') as outfile:
            json.dump([result for result in results if result], outfile, indent=4)

if __name__ == "__main__":
    main()


