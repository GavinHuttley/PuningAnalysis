import numpy as np
import click
import os
import json
from cogent3 import get_app, open_data_store
from cogent3.maths.measure import jsd
import numpy as np
import multiprocessing
from clock_project.data_processing.mafft_align import mafft_aligner
from itertools import combinations

loader = get_app("load_unaligned", format="fasta", moltype="dna")
third_pos = get_app("take_codon_positions", 3)
just_nucs = get_app("omit_degenerates", moltype="dna")
load_json_app = get_app("load_json")

def quantile_calculation(values, quantile):
    values.sort()
    index = int(len(values) * quantile)
    return values[index]

def get_jsd_range(path):
    aln = loader(path)
    nuc_freqs = aln.probs_per_seq()
    jsd_list = []
    species_keys = nuc_freqs.keys()
    for species_1, species_2 in combinations(species_keys, 2):
        jsd_value = jsd(np.array(nuc_freqs[species_1]), np.array(nuc_freqs[species_2]))
        jsd_list.append(jsd_value)

    jsd_range_value = quantile_calculation(jsd_list, 0.95)

    return jsd_range_value


def get_branch_len_diff_range(ens_tree):
    branch_length_range = []
    edge_names = ens_tree.get_tip_names()
    for comb in combinations(edge_names, 2):
        if comb[0] != comb[1]:
            ingroup_species_pair = [comb[0], comb[1]]
            outgroup_species = get_outgroup_species(ingroup_species_pair, ens_tree)
            if outgroup_species:
                triples_species_names = {'ingroup1': ingroup_species_pair[0], 'ingroup2': ingroup_species_pair[1], 'outgroup': outgroup_species[0]}
                sub_tree = ens_tree.get_sub_tree(triples_species_names.values())
                ens_dict = {n: sub_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triples_species_names.values()}
                ens_difference = np.sqrt((ens_dict[ingroup_species_pair[0]]-ens_dict[ingroup_species_pair[1]])**2)
                branch_length_range.append(ens_difference)

    branch_length_range_value = quantile_calculation(branch_length_range, 0.95)
    return branch_length_range_value


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


def get_just3rd_aligned_no_degenerates(seqs):
    """Processes a sequence file to align and filter third codon positions without degenerates."""
    aligned = mafft_aligner(seqs)
    just3rd_aligned = third_pos(aligned)
    just3rd_aligned_no_degenerates = just_nucs(just3rd_aligned)
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
    jsd_dict = {'Ingroup_JSD': ingroup_jsd, 'JSD_difference': {ingroup_species_gene_name[0]: jsd1, ingroup_species_gene_name[1]: jsd2}}
    nuc_freqs_dict = {'Ingroup_nuc_freqs': {ingroup_species_gene_name[0]: list(nuc_freqs1), ingroup_species_gene_name[1]: list(nuc_freqs2)}, 'internal_root_distribution': distribution_internal_root,
        'root_distribution': distribution_root}
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


def initialize_bins(bin_size, max_value):
    return {i: None for i in range(int(max_value / bin_size))}


def process_path(path, output_dir, sequence_dir,alignment_dir, output_alignment_dir):
    gene_name = path.unique_id.split('.')[0]
    print(gene_name)
    result_lf = load_json_app(path).lf
    ens_tree = result_lf.get_ens_tree()
    sequence_path = os.path.join(sequence_dir, f"{gene_name}.fa")
    alignment_path = os.path.join(alignment_dir, f"{gene_name}.fa")
    sequence = loader(sequence_path)
    ens_diff_range = get_branch_len_diff_range(ens_tree)
    jsd_range = get_jsd_range(alignment_path)

    jsd_bins = {'size': jsd_range/500, 'bins': initialize_bins(jsd_range/500, jsd_range)}
    ens_diff_bins = {'size': ens_diff_range/500, 'bins': initialize_bins(ens_diff_range/500, ens_diff_range)}
    output_path = os.path.join(output_dir, gene_name)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    output_alignment_path = os.path.join(output_alignment_dir, gene_name)
    if not os.path.exists(output_alignment_path):
        os.makedirs(output_alignment_path)
    
    out_dstore = open_data_store(output_alignment_path, mode="w", suffix="fa")
    writer = get_app("write_seqs", data_store=out_dstore, format="fasta")

    triples_species_names_dict = {}
    

    max_iterations = 300
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

                if any([processed_jsd, processed_ens_diff]):
                    writer(triple_alignment, identifier = f"{iteration_count}.fa")
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


@click.command()
@click.argument('input_path', type=click.Path(exists=True))
@click.option('-s', '--sequence_dir', type=click.Path(exists=True), required=True, help='Path to the directory containing sequence files.')
@click.option('-a', '--alignment_dir', type=click.Path(exists=True), required=True, help='Path to the directory containing alignment files.')
@click.option('-n', '--num_processes', type=int, default=multiprocessing.cpu_count(), show_default=True, help='Number of processes to use.')
@click.option('-o', '--output_dir', type=click.Path(), help='Path to the out directory for taxanomi triples')
@click.option('-oa', '--output_alignment_dir', type=click.Path(), help='Path to the out directory for taxanomi triples alignments')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process.')
def main(input_path, num_processes, output_dir, sequence_dir, alignment_dir, output_alignment_dir, limit):
    """
    Fit model for sequences in different paths using multiprocessing.
    """
    # Get all paths from input_path with a .json extension
    result_input_dstore = open_data_store(input_path, suffix='json', mode='r', limit=limit)

    # Use multiprocessing to process paths
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir, sequence_dir, alignment_dir, output_alignment_dir) for path in result_input_dstore.completed])

if __name__ == "__main__":
    main()


