import numpy as np
import os
from clock_project.simulation.wts import taxonomic_triple_simulation, generate_ancestor
import json
from cogent3 import get_app, make_tree
import multiprocessing
import click

def model_fitting(aln, opt_args = None):
    tree="('outgroup_edge3',('ingroup_edge1','ingroup_edge2')'internal_node');"
    GN_model = get_app("model", "GN", 
                    tree=tree,
                    opt_args=opt_args,
                    # unique_trees=True,
                    # lf_args=dict(discrete_edges=[['outgroup_edge3']]),
                    time_het="max",
                    optimise_motif_probs=True)
    res = GN_model(aln)
    return res

def get_seed_ancester_seq(pi_low, pi_high, len_short, len_long):
    seq_ls = generate_ancestor(len_short, pi_low)
    seq_ll = generate_ancestor(len_long, pi_low)
    seq_hs = generate_ancestor(len_short, pi_high)
    seq_hl = generate_ancestor(len_long, pi_high)
    seqs = {'low_short': seq_ls, 'high_short': seq_hs, 'low_long': seq_ll, 'high_long': seq_hl}
    return seqs

import os
import numpy as np
from cogent3 import open_data_store
import ast

def simulate_taxanomic_triples(ancestor_seq, Q_collection, pi, seed, t, path_to_dir):
    ens_dict = {}
    length = len(ancestor_seq)
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)
    for matrix_dict in Q_collection:
        Q_ingroup1, Q_ingroup2, Q_outgroup = matrix_dict.values()
        aln, ens = taxonomic_triple_simulation(pi, np.array(Q_ingroup1), np.array(Q_ingroup2), np.array(Q_outgroup), 0, t, length, 1, seed, ancestor_seq)
        # Open a data store for each simulation
        number = Q_collection.index(matrix_dict)
        ens_dict[number] = ens
        write_json_app(aln, identifier=f'{number}.json')
    with open (f'{path_to_dir}_ens_dict.json', 'w') as outfile:
        json.dump(ens_dict, outfile, indent=4)

import itertools
def generate_combinations(length_list, pi_list):
    # Create all combinations of time, length, and pi
    combinations = itertools.product(length_list, pi_list)
    return list(combinations)

def processing_parameters(combo, t_list, Q_collection, seed, output_dir):
    length, pi = combo
    pi_indicate = 'low' if pi == [0.25, 0.25, 0.25, 0.25] else 'high'
    ancestor_seq = generate_ancestor(length, pi)  # Assuming this function is defined elsewhere
    ans_seq_list = [str(num) for num in ancestor_seq]
    with open (f'{pi_indicate}_{length}_.json', 'w') as outfile:
        json.dump(ans_seq_list, outfile)
    for t in t_list: 
        path_to_dir = os.path.join(output_dir, f'aln_{pi_indicate}_{length}_{t}')
        os.makedirs(path_to_dir, exist_ok=True)
        simulate_taxanomic_triples(ancestor_seq, Q_collection, pi, seed, t, path_to_dir)


@click.command()
@click.option('--time_list', '-t', default='[0, 1]', callback=lambda _, __, x: ast.literal_eval(x),
              help='List of times to simulate, e.g., "[0, 1]"')
@click.option('--length_list', '-l', default='[100, 1000]', callback=lambda _, __, x: ast.literal_eval(x),
              help='List of sequence lengths to simulate, e.g., "[100, 1000]"')
@click.option('--pi_list', '-pi', default='[[0.25, 0.25, 0.25, 0.25], [0.49, 0.02, 0.35, 0.14]]', callback=lambda _, __, x: ast.literal_eval(x),
              help='List of initial nucleotide distributions, e.g., "[[0.25, 0.25, 0.25, 0.25], [0.49, 0.02, 0.35, 0.14]]"')
@click.option('--matrix_collection_dir', '-md', type=click.Path(), help='matrix_collection directory')
@click.option('--seed', '-seed', type=int, help='Seed for random number generation')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory')
def main(time_list, length_list, pi_list, matrix_collection_dir, seed, num_processes, output_dir):
    with open (matrix_collection_dir, 'r') as infile:
        matrices_dict = json.load(infile)
    matrix_collection = []
    for gene, matrices_list in matrices_dict.items():
        matrix_collection.extend(matrices_list)
    combinations = generate_combinations(length_list, pi_list)
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(processing_parameters, [(combo, time_list, matrix_collection, seed, output_dir) for combo in combinations])

if __name__ == "__main__":
    main()





