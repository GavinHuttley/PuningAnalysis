import numpy as np
import argparse
import os
import json
from cogent3 import get_app
from cogent3.util.deserialise import deserialise_object
from cogent3.maths.measure import jsd
import numpy as np
import multiprocessing
import glob
import scipy
from clock_project.simulation.wts import calculate_non_stationarity
from clock_project.genome_analysis.homolog_analysis import cpos3, codon_aligner



model_fitting_result_dir = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result_long_seq'
alignment_dir = '/Users/gulugulu/repos/PuningAnalysis/data/ensembl_ortholog_sequences/homologies_alignments_common_name_2'


def get_ingroup_pair(sequence):
    species_names = sequence.names
    ingroup_species_pair = np.random.choice(np.array(species_names), 2, replace=False)

    return list(ingroup_species_pair)

def get_relationship(ingroup_species_pair, ens_tree):
    lca_node = ens_tree.lowest_common_ancestor(ingroup_species_pair).get_node_names()[0]
    parents = ens_tree.lowest_common_ancestor(ingroup_species_pair).parent
    if parents:
        parent_node = parents.get_node_names()[0]
        parent_children = ens_tree.get_nodes_dict()[parent_node].children
        relation_dict = {'sibling': [], 'cousion': []}
        for child in parent_children:
            grad_child = child.get_node_names()
            if grad_child[0] == lca_node:
                relation_dict['sibling'].extend(child.get_tip_names())
            else:
                relation_dict['cousion'].extend(child.get_tip_names())
        return relation_dict
    else:
        return None
    
def get_outgroup_species(ingroup_species_pair, ens_tree):
    relationship = get_relationship(ingroup_species_pair, ens_tree)
    if relationship: 
        cousion = relationship['cousion']
        outgroup_species = np.random.choice(np.array(cousion), 1, replace=False)
        return list(outgroup_species)
    else: 
        return None
    
def get_trace_back_nodes(ingroup_species_gene_name, outgroup, lca_node, parent_node, ens_tree):
    trace_back_nodes = {}
    for species in ingroup_species_gene_name:
        path = []
        for tree in ens_tree.get_connecting_edges(lca_node, species):
            internal_root = tree.get_node_names()[0]
            path.append(internal_root)
        trace_back_nodes[species] = path

    outgroup_path = []
    for tree in ens_tree.get_connecting_edges(parent_node, outgroup):
        internal_root = tree.get_node_names()[0]
        outgroup_path.append(internal_root)
    trace_back_nodes[outgroup] = outgroup_path

    return trace_back_nodes

def get_Q_matrices(triads_species_names, lca_node, parent_node, ens_tree, result_lf, ens_dict):
    ingroup_species = [triads_species_names['ingroup1'], triads_species_names['ingroup2']]
    outgroup = triads_species_names['outgroup']
    trace_back_nodes = get_trace_back_nodes(ingroup_species, outgroup, lca_node, parent_node, ens_tree)
    branch_length_dict = {species: {n: ens_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in trace_back_nodes[species]} for species in triads_species_names.values()}    
    matrices_dict = {species: {n: result_lf.get_rate_matrix_for_edge(n, calibrated=True).to_array() for n in trace_back_nodes[species]} for species in triads_species_names.values()}
    p_matrices = {}
    for species in triads_species_names.values(): 
        p_matrices[species] = {}
        for node_name, matrix in matrices_dict[species].items():
            p_matrices[species][node_name] = scipy.linalg.expm(matrix * branch_length_dict[species][node_name])
    result_p_matrices = {}
    for species in triads_species_names.values(): 
        result_p_matrix = list(p_matrices[species].values())[0]
        for matrix in list(p_matrices[species].values())[1:]:
            result_p_matrix = np.dot(result_p_matrix, matrix)
        result_p_matrices[species] = result_p_matrix
    result_Qt_matrices = {species: scipy.linalg.logm(matrix) for species, matrix in result_p_matrices.items()}
    result_Q_matrices = {species: result_Qt_matrices[species]/ens_dict[species] for species in triads_species_names.values()}

    return result_Q_matrices


def get_just3rd_aligned_no_degenerates(seqs):
    """Processes a sequence file to align and filter third codon positions without degenerates."""
    aligned = codon_aligner(seqs)
    just3rd_aligned = cpos3(aligned)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates(motif_length=3)
    return just3rd_aligned_no_degenerates


def get_triads_info(triads_species_names, ens_tree, sequence, result_lf, lca_node, parent_node, iteration_count):
    triad_sequence = sequence.take_seqs(triads_species_names.values())
    triad_alignment = get_just3rd_aligned_no_degenerates(triad_sequence)
    sub_tree = ens_tree.get_sub_tree(triads_species_names.values())
    ingroup_species_gene_name = [triads_species_names['ingroup1'], triads_species_names['ingroup2']]
    outgroup = triads_species_names['outgroup']
    ingroup_seqs = triad_alignment.take_seqs(ingroup_species_gene_name)
    nuc_freqs = ingroup_seqs.probs_per_seq()
    nuc_freqs1 = nuc_freqs[ingroup_species_gene_name[0]]
    nuc_freqs2 = nuc_freqs[ingroup_species_gene_name[1]]
    ingroup_jsd = jsd(nuc_freqs1, nuc_freqs2)
    distribution_internal_root = list(result_lf.get_motif_probs_by_node()[lca_node])
    jsd1 = jsd(nuc_freqs1, distribution_internal_root)
    jsd2 = jsd(nuc_freqs2, distribution_internal_root)
    distribution_root = list(result_lf.get_motif_probs_by_node()[parent_node])
    ens_dict = {n: sub_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triads_species_names.values()}
    ens_difference = np.sqrt((ens_dict[ingroup_species_gene_name[0]]-ens_dict[ingroup_species_gene_name[1]])**2)
    shortest_ens = min(ens_dict.values())  

    matrices = get_Q_matrices(triads_species_names, lca_node, parent_node, ens_tree, result_lf, ens_dict)
    jsd_dict = {'Ingroup_JSD': ingroup_jsd, 'JSD_difference': {ingroup_species_gene_name[0]: jsd1, ingroup_species_gene_name[1]: jsd2}}
    nuc_freqs_dict = {'Ingroup_nuc_freqs': {ingroup_species_gene_name[0]: list(nuc_freqs1), ingroup_species_gene_name[1]: list(nuc_freqs2)}, 'internal_root_distribution': distribution_internal_root,
        'root_distribution': distribution_root}
    nabla_dict = {}
    for species in ingroup_species_gene_name:
        nabla_dict[species] = calculate_non_stationarity(distribution_internal_root, matrices[species], ens_dict[species])
    
    nabla_dict[outgroup] = calculate_non_stationarity(distribution_root, matrices[outgroup], ens_dict[outgroup])
    triad_data = {'identifier': iteration_count, 'triads_species_names': triads_species_names, 'triads_info_big_tree': 
        {'internal_node': lca_node,
        'root': parent_node,
        'jsd_dict': jsd_dict,
        'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'nuc_freqs_dict': nuc_freqs_dict,
        'nabla_values': nabla_dict}}
    
    return triad_data, triad_alignment


def triads_binning_jsd(triad_data, jsd_bins):
    processed_jsd = False
    
    triad_values  = triad_data['triads_info_big_tree']
    jsd_bin_index = int(triad_values['jsd_dict']['Ingroup_JSD'] // jsd_bins['size'])

    if jsd_bin_index in jsd_bins['bins'] and jsd_bins['bins'][jsd_bin_index] is None:
        jsd_bins['bins'][jsd_bin_index] = triad_data
        processed_jsd = True

    return jsd_bins, processed_jsd

def triads_binning_ens_diff(triad_data, ens_diff_bins):
    processed_ens_diff = False
    
    triad_values  = triad_data['triads_info_big_tree']
    ens_diff_bin_index = int(triad_values['ens_difference'] // ens_diff_bins['size'])

    if ens_diff_bin_index in ens_diff_bins['bins'] and ens_diff_bins['bins'][ens_diff_bin_index] is None:
        ens_diff_bins['bins'][ens_diff_bin_index] = triad_data
        processed_ens_diff = True

    return ens_diff_bins, processed_ens_diff

def triads_binning_shortest_ens(triad_data, shortest_ens_bins):
    processed_shor_ens = False
    
    triad_values  = triad_data['triads_info_big_tree']
    shortest_ens_bin_index = int(triad_values['shortest_ingroup_ens'] // shortest_ens_bins['size'])

    if shortest_ens_bin_index in shortest_ens_bins['bins'] and shortest_ens_bins['bins'][shortest_ens_bin_index] is None:
        shortest_ens_bins['bins'][shortest_ens_bin_index] = triad_data
        processed_shor_ens = True

    return shortest_ens_bins, processed_shor_ens

def initialize_bins(bin_size, max_value):
    return {i: None for i in range(int(max_value / bin_size))}




def process_path(path, output_dir, sequence_dir, output_alignment_dir):
    file_name = os.path.basename(path).rsplit('.', 1)[0]
    print(file_name)
    load_json_app = get_app("load_json")
    result_lf = load_json_app(path)
    ens_tree = result_lf.get_ens_tree()
    sequence_path = os.path.join(sequence_dir, f"{file_name}.json")
    sequence = deserialise_object(json.load(open(sequence_path, 'r')))

    jsd_bins = {'size': 0.00004, 'bins': initialize_bins(0.00004, 0.02)}
    ens_diff_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 0.5)}
    shortest_ens_bins = {'size': 0.0004, 'bins': initialize_bins(0.0004, 0.2)}
    output_path = os.path.join(output_dir, file_name)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_alignment_path = os.path.join(output_alignment_dir, file_name)
    if not os.path.exists(output_alignment_path):
        os.makedirs(output_alignment_path)

    triads_species_names_dict = {}
    

    

    max_iterations = 500
    iteration_count = 0
    not_processed = 0
    matrix_error = 0

    while iteration_count < max_iterations:
        ingroup_species_pair = get_ingroup_pair(sequence)
        outgroup_species = get_outgroup_species(ingroup_species_pair, ens_tree)
        if outgroup_species:
            triads_species_names = {'ingroup1': ingroup_species_pair[0], 'ingroup2': ingroup_species_pair[1], 'outgroup': outgroup_species[0]}
            parent_node = ens_tree.lowest_common_ancestor(ingroup_species_pair).parent.get_node_names()[0]
            lca_node = ens_tree.lowest_common_ancestor(ingroup_species_pair).get_node_names()[0]
            try: 
                triat_data, triad_alignment = get_triads_info(triads_species_names, ens_tree, sequence, result_lf, lca_node, parent_node, iteration_count)
                jsd_bins, processed_jsd = triads_binning_jsd(triat_data, jsd_bins)
                ens_diff_bins, processed_ens_diff = triads_binning_ens_diff(triat_data, ens_diff_bins)
                shortest_ens_bins, processed_shor_ens = triads_binning_shortest_ens(triat_data, shortest_ens_bins)

                if any([processed_jsd, processed_shor_ens, processed_ens_diff]):
                    with open(os.path.join(output_alignment_path, f"{iteration_count}.json"), 'w') as f:
                        json.dump(triad_alignment.to_json(), f, indent=4)
                    
                    triads_species_names_dict[iteration_count] = triads_species_names

                if not processed_shor_ens and not processed_ens_diff and not processed_jsd:
                    not_processed += 1

            except Exception as e:
                matrix_error += 1
                

        iteration_count += 1   
    print(matrix_error)
    print(not_processed)


    with open(os.path.join(output_path, 'triads_species_names_dict.json'), 'w') as f:
        json.dump(triads_species_names_dict, f, indent=4)
    with open(os.path.join(output_path, 'jsd_bins.json'), 'w') as f:
        json.dump(jsd_bins, f, indent=4)
    with open(os.path.join(output_path, 'ens_diff_bins.json'), 'w') as f:
        json.dump(ens_diff_bins, f, indent=4)
    with open(os.path.join(output_path, 'shortest_ens_bins.json'), 'w') as f:
        json.dump(shortest_ens_bins, f, indent=4)

def main(paths, num_processes, output_dir, num_file, sequence_dir, output_alignment_dir):
    if num_file is not None:
        paths = paths[:num_file]
    else:
        paths = paths
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir, sequence_dir, output_alignment_dir) for path in paths])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit model for sequences in different paths.')
    parser.add_argument('input', type=str, help='Path to the directory containing model fitting result')
    parser.add_argument('-s', '--sequence_dir', type=str, help='Path to the directory containing sequence')
    parser.add_argument('-n', '--num_processes', type=int, default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
    parser.add_argument('-o', '--output_dir', type=str, default='/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads', help='Output directory for JSON files (default: /Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result)')
    parser.add_argument('-oa', '--output_alignment_dir', type=str, default='/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_alignments')
    parser.add_argument('-l', '--limit', type=int, help='Limit the number of files to process')


    args = parser.parse_args()
    input_path = args.input
    num_processes = args.num_processes
    output_dir = args.output_dir
    num_file = args.limit
    sequence_dir = args.sequence_dir
    output_alignment_dir = args.output_alignment_dir
    
    if os.path.isdir(input_path):
        paths = glob.glob(os.path.join(input_path, "*.json"))

    main(paths, num_processes, output_dir, num_file, sequence_dir, output_alignment_dir)

# @click.command()
# @click.argument('input', type=click.Path(exists=True))
# @click.option('--alignment_dir', '-a', type=click.Path(exists=True), help='Path to the directory containing alignments')
# @click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
# @click.option('--output_dir', '-o', default='/path/to/default/output', type=click.Path(), help='Output directory for files')
# @click.option('--output_alignment_dir', '-oa', default='/path/to/default/alignment_output', type=click.Path(), help='Output directory for alignments')
# @click.option('--limit', '-l', type=int, help='Limit the number of files to process')
# def main(input, alignment_dir, num_processes, output_dir, output_alignment_dir, limit):
#     """Fit model for sequences in different paths."""
#     if os.path.isdir(input):
#         paths = glob.glob(os.path.join(input, "*.json")
#     else:
#         paths = [input]

#     if limit is not None:
#         paths = paths[:limit]

#     with multiprocessing.Pool(processes=num_processes) as pool:
#         pool.starmap(process_path, [(path, output_dir, alignment_dir, output_alignment_dir) for path in paths])

# if __name__ == "__main__":
#     main()