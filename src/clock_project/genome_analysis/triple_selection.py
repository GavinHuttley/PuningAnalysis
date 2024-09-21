import numpy as np
import click
import os
import json
from cogent3 import get_app, open_data_store
from cogent3.maths.measure import jsd
import numpy as np
import multiprocessing
from clock_project.genome_analysis.sequence_alignment_filtering import cpos3, aligner

load_json_app = get_app("load_json")
# model_fitting_result_dir = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result_long_seq'
# alignment_dir = '/Users/gulugulu/repos/PuningAnalysis/data/ensembl_ortholog_sequences/homologies_alignments_common_name_2'


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

# def get_Q_matrices(triples_species_names, lca_node, parent_node, ens_tree, result_lf, ens_dict):
#     ingroup_species = [triples_species_names['ingroup1'], triples_species_names['ingroup2']]
#     outgroup = triples_species_names['outgroup']
#     trace_back_nodes = get_trace_back_nodes(ingroup_species, outgroup, lca_node, parent_node, ens_tree)
#     branch_length_dict = {species: {n: ens_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in trace_back_nodes[species][1:]} for species in triples_species_names.values()}    
#     matrices_dict = {species: {n: result_lf.get_rate_matrix_for_edge(n, calibrated=True).to_array() for n in trace_back_nodes[species][1:]} for species in triples_species_names.values()}
#     p_matrices = {}
#     for species in triples_species_names.values(): 
#         p_matrices[species] = {}
#         for node_name, matrix in matrices_dict[species].items():
#             p_matrices[species][node_name] = scipy.linalg.expm(matrix * branch_length_dict[species][node_name])
#     result_p_matrices = {}
#     for species in triples_species_names.values(): 
#         result_p_matrix = list(p_matrices[species].values())[0]
#         for matrix in list(p_matrices[species].values())[1:]:
#             result_p_matrix = np.dot(result_p_matrix, matrix)
#         result_p_matrices[species] = result_p_matrix
#     result_Qt_matrices = {species: scipy.linalg.logm(matrix) for species, matrix in result_p_matrices.items()}
#     result_Q_matrices = {species: result_Qt_matrices[species]/ens_dict[species] for species in triples_species_names.values()}

#     return result_Q_matrices


def get_just3rd_aligned_no_degenerates(seqs):
    """Processes a sequence file to align and filter third codon positions without degenerates."""
    aligned = aligner(seqs)
    just3rd_aligned = cpos3(aligned)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates()
    return just3rd_aligned_no_degenerates

def change_nucleotide_order(motif_prob):
    nuc_freq = {'A': motif_prob[0],'C': motif_prob[1],'G': motif_prob[2],'T': motif_prob[3]}
    nuc_freq_reordered = {'T': nuc_freq['T'], 'C': nuc_freq['C'], 'A': nuc_freq['A'], 'G': nuc_freq['G']}
    return list(nuc_freq_reordered.values())

def get_triples_info(triples_species_names, ens_tree, sequence, result_lf, lca_node, parent_node, iteration_count):
    triple_sequence = sequence.take_seqs(triples_species_names.values())
    triple_alignment = get_just3rd_aligned_no_degenerates(triple_sequence)
    sub_tree = ens_tree.get_sub_tree(triples_species_names.values())
    ingroup_species_gene_name = [triples_species_names['ingroup1'], triples_species_names['ingroup2']]
    outgroup = triples_species_names['outgroup']
    ingroup_seqs = triple_alignment.take_seqs(ingroup_species_gene_name)
    nuc_freqs = ingroup_seqs.probs_per_seq()
    nuc_freqs1 = nuc_freqs[ingroup_species_gene_name[0]]
    nuc_freqs2 = nuc_freqs[ingroup_species_gene_name[1]]
    nuc_freqs1_reordered = change_nucleotide_order(nuc_freqs1)
    nuc_freqs2_reordered = change_nucleotide_order(nuc_freqs2)
    ingroup_jsd = jsd(nuc_freqs1_reordered, nuc_freqs2_reordered)
    distribution_internal_root = list(result_lf.get_motif_probs_by_node()[lca_node])
    jsd1 = jsd(nuc_freqs1_reordered, distribution_internal_root)
    jsd2 = jsd(nuc_freqs2_reordered, distribution_internal_root)
    distribution_root = list(result_lf.get_motif_probs_by_node()[parent_node])
    ens_dict = {n: sub_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triples_species_names.values()}
    ens_difference = np.sqrt((ens_dict[ingroup_species_gene_name[0]]-ens_dict[ingroup_species_gene_name[1]])**2)
    shortest_ens = min(ens_dict.values())

    # matrices = get_Q_matrices(triples_species_names, lca_node, parent_node, ens_tree, result_lf, ens_dict)
    jsd_dict = {'Ingroup_JSD': ingroup_jsd, 'JSD_difference': {ingroup_species_gene_name[0]: jsd1, ingroup_species_gene_name[1]: jsd2}}
    nuc_freqs_dict = {'Ingroup_nuc_freqs': {ingroup_species_gene_name[0]: list(nuc_freqs1), ingroup_species_gene_name[1]: list(nuc_freqs2)}, 'internal_root_distribution': distribution_internal_root,
        'root_distribution': distribution_root}
    # nabla_dict = {}
    # for species in ingroup_species_gene_name:
    #     nabla_dict[species] = calculate_non_stationarity(distribution_internal_root, matrices[species], ens_dict[species])
    
    # nabla_dict[outgroup] = calculate_non_stationarity(distribution_root, matrices[outgroup], ens_dict[outgroup])
    triple_data = {'identifier': iteration_count, 'triples_species_names': triples_species_names, 'triples_info_big_tree': 
        {'internal_node': lca_node,
        'root': parent_node,
        'jsd_dict': jsd_dict,
        'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'nuc_freqs_dict': nuc_freqs_dict}}
    
    return triple_data, triple_alignment


def triples_binning_jsd(triple_data, jsd_bins):
    processed_jsd = False
    
    triple_values  = triple_data['triples_info_big_tree']
    jsd_bin_index = int(triple_values['jsd_dict']['Ingroup_JSD'] // jsd_bins['size'])

    if jsd_bin_index in jsd_bins['bins'] and jsd_bins['bins'][jsd_bin_index] is None:
        jsd_bins['bins'][jsd_bin_index] = triple_data
        processed_jsd = True

    return jsd_bins, processed_jsd

def triples_binning_ens_diff(triple_data, ens_diff_bins):
    processed_ens_diff = False
    
    triple_values  = triple_data['triples_info_big_tree']
    ens_diff_bin_index = int(triple_values['ens_difference'] // ens_diff_bins['size'])

    if ens_diff_bin_index in ens_diff_bins['bins'] and ens_diff_bins['bins'][ens_diff_bin_index] is None:
        ens_diff_bins['bins'][ens_diff_bin_index] = triple_data
        processed_ens_diff = True

    return ens_diff_bins, processed_ens_diff

# def triples_binning_shortest_ens(triple_data, shortest_ens_bins):
#     processed_shor_ens = False
    
#     triple_values  = triple_data['triples_info_big_tree']
#     shortest_ens_bin_index = int(triple_values['shortest_ingroup_ens'] // shortest_ens_bins['size'])

#     if shortest_ens_bin_index in shortest_ens_bins['bins'] and shortest_ens_bins['bins'][shortest_ens_bin_index] is None:
#         shortest_ens_bins['bins'][shortest_ens_bin_index] = triple_data
#         processed_shor_ens = True

#     return shortest_ens_bins, processed_shor_ens

def initialize_bins(bin_size, max_value):
    return {i: None for i in range(int(max_value / bin_size))}


def process_path(path, output_dir, sequence_dir, output_alignment_dir):
    gene_name = path.unique_id.split('.')[0]
    print(gene_name)
    result_lf = load_json_app(path).lf
    ens_tree = result_lf.get_ens_tree()
    sequence_path = os.path.join(sequence_dir, f"{gene_name}.json")
    sequence = load_json_app(sequence_path)

    jsd_bins = {'size': 0.00004, 'bins': initialize_bins(0.00004, 0.02)}
    ens_diff_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 0.5)}
    # shortest_ens_bins = {'size': 0.0004, 'bins': initialize_bins(0.0004, 0.2)}
    output_path = os.path.join(output_dir, gene_name)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_alignment_path = os.path.join(output_alignment_dir, gene_name)
    if not os.path.exists(output_alignment_path):
        os.makedirs(output_alignment_path)
    
    out_dstore = open_data_store(output_alignment_path, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)

    triples_species_names_dict = {}
    

    max_iterations = 500
    iteration_count = 0
    not_processed = 0
    matrix_error = 0

    while iteration_count < max_iterations:
        ingroup_species_pair = get_ingroup_pair(sequence)
        outgroup_species = get_outgroup_species(ingroup_species_pair, ens_tree)
        if outgroup_species:
            triples_species_names = {'ingroup1': ingroup_species_pair[0], 'ingroup2': ingroup_species_pair[1], 'outgroup': outgroup_species[0]}
            parent_node = ens_tree.lowest_common_ancestor(ingroup_species_pair).parent.get_node_names()[0]
            lca_node = ens_tree.lowest_common_ancestor(ingroup_species_pair).get_node_names()[0]
            try: 
                triat_data, triple_alignment = get_triples_info(triples_species_names, ens_tree, sequence, result_lf, lca_node, parent_node, iteration_count)
                jsd_bins, processed_jsd = triples_binning_jsd(triat_data, jsd_bins)
                ens_diff_bins, processed_ens_diff = triples_binning_ens_diff(triat_data, ens_diff_bins)
                # shortest_ens_bins, processed_shor_ens = triples_binning_shortest_ens(triat_data, shortest_ens_bins)

                if any([processed_jsd, processed_ens_diff]):
                    write_json_app(triple_alignment, identifier = f"{iteration_count}.json")
                    triples_species_names_dict[iteration_count] = triples_species_names

                if not processed_ens_diff and not processed_jsd:
                    not_processed += 1

            except Exception as e:
                matrix_error += 1

        iteration_count += 1   


    with open(os.path.join(output_path, 'triples_species_names_dict.json'), 'w') as f:
        json.dump(triples_species_names_dict, f, indent=4)
    with open(os.path.join(output_path, 'jsd_bins.json'), 'w') as f:
        json.dump(jsd_bins, f, indent=4)
    with open(os.path.join(output_path, 'ens_diff_bins.json'), 'w') as f:
        json.dump(ens_diff_bins, f, indent=4)
    # with open(os.path.join(output_path, 'shortest_ens_bins.json'), 'w') as f:
    #     json.dump(shortest_ens_bins, f, indent=4)

@click.command()
@click.argument('input_path', type=click.Path(exists=True))
@click.option('-s', '--sequence_dir', type=click.Path(exists=True), required=True, help='Path to the directory containing sequence files.')
@click.option('-n', '--num_processes', type=int, default=multiprocessing.cpu_count(), show_default=True, help='Number of processes to use.')
@click.option('-o', '--output_dir', type=click.Path(), default='/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triples', show_default=True, help='Output directory for JSON files.')
@click.option('-oa', '--output_alignment_dir', type=click.Path(), default='/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triples_alignments', show_default=True, help='Output directory for alignment files.')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process.')
def main(input_path, num_processes, output_dir, sequence_dir, output_alignment_dir, limit):
    """
    Fit model for sequences in different paths using multiprocessing.
    """
    # Get all paths from input_path with a .json extension
    result_input_dstore = open_data_store(input_path, suffix='json', mode='r', limit=limit)

    # Use multiprocessing to process paths
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir, sequence_dir, output_alignment_dir) for path in result_input_dstore.completed])

if __name__ == "__main__":
    main()

