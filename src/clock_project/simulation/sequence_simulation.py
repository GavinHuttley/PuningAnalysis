import numpy as np
import os
from clock_project.simulation.wts import taxonomic_triple_simulation, generate_ancestor
import json
from cogent3 import get_app
import multiprocessing
import click

# with open ('/Users/gulugulu/Desktop/honours/data_local_2/seed_ancester_sequences.json', 'r') as infile:
#     seqs_list_original = json.load(infile)

# seqs_list1 = {'low_short': seqs_list_original['low_short'], 'high_short': seqs_list_original['high_short']}

# def get_seed_ancester_seq(pi_low, pi_high, len_long):
#     seq_ll = generate_ancestor(len_long, pi_low)
#     seq_hl = generate_ancestor(len_long, pi_high)
#     seqs = {'low_long': seq_ll, 'high_long': seq_hl}
#     return seqs

# seqs = get_seed_ancester_seq([0.25, 0.25, 0.25, 0.25], [0.06, 0.47, 0.08, 0.39], 3000)
# seqs_list = {key:[str(num) for num in seqs[key]] for key in seqs.keys()}

# seqs_list1.update(seqs_list)

# with open ('/Users/gulugulu/Desktop/honours/data_local_2/seed_ancester_sequences_3.json', 'w') as outfile:
#     json.dump(seqs_list1, outfile, indent=4)

import os
import numpy as np
from cogent3 import open_data_store
import ast

seed_ancestor_seqs_path = 'input/seed_ancester_sequences.json'
seed_ancestor_seqs_list = json.load(open(seed_ancestor_seqs_path, 'r'))
seed_ancestor_seqs = {key: [int(str) for str in seed_ancestor_seqs_list[key]] for key in seed_ancestor_seqs_list.keys()}


def simulate_taxonomic_triples(ancestor_seq, Q_group, pi, seed, t, path_to_dir, ens_dict):
    length = len(ancestor_seq)
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)
    Q_ingroup1, Q_ingroup2, Q_outgroup = Q_group
    aln, ens = taxonomic_triple_simulation(
        pi, np.array(Q_ingroup1), np.array(Q_ingroup2), np.array(Q_outgroup), 
        0, t, length, 1, seed, ancestor_seq
    )
    
    # Write alignment to a JSON file named after the process number
    process_id = os.getpid()  # Use process ID for uniqueness
    write_json_app(aln, identifier=f'{process_id}.json')
    
    # Store ens in the shared dictionary
    ens_dict[process_id] = ens

import itertools
def generate_combinations(length_list, pi_list):
    # Create all combinations of time, length, and pi
    combinations = itertools.product(length_list, pi_list)
    return list(combinations)

def processing_parameters(ancestor_seq, length, pi, t, Q_collection, seed, output_dir):
    pi_indicate = 'low' if pi == [0.25, 0.25, 0.25, 0.25] else 'high'
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





