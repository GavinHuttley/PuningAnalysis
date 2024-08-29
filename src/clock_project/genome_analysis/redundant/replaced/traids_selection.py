import argparse
import os
import json
import random
from cogent3 import get_app, make_tree
from cogent3.util.deserialise import deserialise_object
from clock_project.genome_analysis.sequence_alignment_filtering import codon_aligner, cpos3
from clock_project.simulation.wts import calculate_non_stationarity
from cogent3.maths.measure import jsd
import numpy as np
from multiprocessing import Pool


def extract_taxonomic_order_info(taxa_reference, order_reference, species_names):
    taxonomic_info = {species: {'order': taxa_reference[species]['order'],
                                'family': taxa_reference[species]['family'],
                                'genus': taxa_reference[species]['genus']}
                      for species in species_names if species in taxa_reference}

    order_info = {}
    for order, species_list in order_reference.items():
        filtered_species = [species for species in species_list if species in species_names]
        if filtered_species:
            order_info[order] = filtered_species

    return taxonomic_info, order_info

def match_species_gene_name(species_name, name_dict):
    return name_dict.get(species_name)

def get_sequence_alignment(triad_species_name, sequences, name_dict):
    species_gene_names = [match_species_gene_name(name, name_dict) for name in triad_species_name.values()]
    triad_sequences = sequences.take_seqs(species_gene_names)
    triad_sequences_no_stop_codon = triad_sequences.trim_stop_codons(strict=True)
    alignment = codon_aligner(triad_sequences_no_stop_codon)
    just3rd_aligned = cpos3(alignment)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates(motif_length=3)
    return just3rd_aligned_no_degenerates

def set_model(triats_species_name, model_name, name_dict):
    species_gene_names = {}
    for key, species_name in triats_species_name.items():
        species_gene_name = match_species_gene_name(species_name, name_dict)
        species_gene_names[key] = species_gene_name
    tree_topology = make_tree(f"(({species_gene_names['ingroup1']},{species_gene_names['ingroup2']}),{species_gene_names['outgroup']})")
    model = get_app("model", sm=model_name, time_het="max", tree=tree_topology, show_progress=False)
    #  model = get_app("model", sm='GN', show_progress=False)
    return model, species_gene_names

def extract_species_names(seqs):
    seq_names = seqs.names
    species_names = []
    for name in seq_names:
        species_info, gene_info = name.split('-', 1)
        species_names.append(species_info)
    return species_names

def pseudorandom_sequence_selection_between_order(species_names, taxa_info, order_info):
    species_array = np.array(species_names)
    ingroup_species_pair = np.random.choice(species_array, 2, replace=False)
    
    order1 = taxa_info[ingroup_species_pair[0]]['order']
    order2 = taxa_info[ingroup_species_pair[1]]['order']
    different_order_species = [s for s, taxa in taxa_info.items() if taxa['order'] not in {order1, order2}]

    if different_order_species:
        outgroup_species = np.random.choice(different_order_species)

    triad_species_name = {
        'ingroup1': ingroup_species_pair[0],
        'ingroup2': ingroup_species_pair[1],
        'outgroup': outgroup_species
    }
    return triad_species_name

def pseudorandom_sequence_selection_within_order(species_names,taxa_reference, order_reference):
    _, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)
    suitable_orders = [order for order in order_info if len(order_info[order]) > 1]
    if not suitable_orders:
        return None
    random_order = random.choice(suitable_orders)
    ingroup_species_pair = random.sample(order_info[random_order], 2)

    other_orders = {key: value for key, value in order_info.items() if key != random_order}

    other_species = [species for species_list in other_orders.values() for species in species_list]

    outgroup_species = random.choice(other_species) if other_species else None

    triad_species_name = {
        'ingroup1': ingroup_species_pair[0],
        'ingroup2': ingroup_species_pair[1],
        'outgroup': outgroup_species
    }

    return triad_species_name


def get_triads_info(triads_species_name, model_name, sequences, name_dict):
    model, species_gene_names = set_model(triads_species_name, model_name, name_dict)
    triads_alignment = get_sequence_alignment(triads_species_name, sequences, name_dict)
    ingroup_species_gene_name = [species_gene_names['ingroup1'], species_gene_names['ingroup2']]
    ingroup_seqs = triads_alignment.take_seqs(ingroup_species_gene_name)
    nuc_freqs = ingroup_seqs.probs_per_seq()
    ingroup_jsd = jsd(nuc_freqs[ingroup_species_gene_name[0]], nuc_freqs[ingroup_species_gene_name[1]])
    model_fitting_result = model(triads_alignment)  
    matrices = {n: model_fitting_result.lf.get_rate_matrix_for_edge(n, calibrated=True) for n in species_gene_names.values()}
    initial_distribution = list(model_fitting_result.lf.get_motif_probs())
    ens_dict = model_fitting_result.lf.get_lengths_as_ens()
    nabla_dict = {n: calculate_non_stationarity(initial_distribution, matrices[n], ens_dict[n]) for n in species_gene_names.values()}
    ens_difference = np.sqrt((ens_dict[ingroup_species_gene_name[0]]-ens_dict[ingroup_species_gene_name[1]])**2)
    shortest_ens = min(ens_dict.values())  
    triad_data = {'triads_species_name': triads_species_name, 'triads_info': 
        {'ingroup_jsd': ingroup_jsd,
        'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'nabla_info': nabla_dict,
        'ens': ens_dict,
        'initial_distribution': initial_distribution}
    }
    return triad_data, model_fitting_result


def triads_binning_jsd(triad_data, jsd_bins):
    
    triad_values  = triad_data['triads_info']
    # Calculate bin indices for each property
    jsd_bin_index = int(triad_values['ingroup_jsd'] // jsd_bins['size'])

    if jsd_bin_index in jsd_bins['bins'] and jsd_bins['bins'][jsd_bin_index] is None:
        jsd_bins['bins'][jsd_bin_index] = triad_data

    return jsd_bins

def triads_binning_ens_diff(triad_data, ens_diff_bins):
    
    triad_values  = triad_data['triads_info']
    ens_diff_bin_index = int(triad_values['ens_difference'] // ens_diff_bins['size'])

    if ens_diff_bin_index in ens_diff_bins['bins'] and ens_diff_bins['bins'][ens_diff_bin_index] is None:
        ens_diff_bins['bins'][ens_diff_bin_index] = triad_data

    return ens_diff_bins

def triads_binning_shortest_ens(triad_data, shortest_ens_bins):
    
    triad_values  = triad_data['triads_info']
    shortest_ens_bin_index = int(triad_values['shortest_ingroup_ens'] // shortest_ens_bins['size'])

    if shortest_ens_bin_index in shortest_ens_bins['bins'] and shortest_ens_bins['bins'][shortest_ens_bin_index] is None:
        shortest_ens_bins['bins'][shortest_ens_bin_index] = triad_data

    return shortest_ens_bins


def initialize_bins(bin_size, max_value):
    return {i: None for i in range(int(max_value / bin_size))}


def load_taxa_order_reference():
    with open('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_taxa_info.json', 'r') as infile1:
        taxa_reference = json.load(infile1)
    with open('/Users/gulugulu/repos/PuningAnalysis/results/output_data/species_order_info.json', 'r') as infile2:
        order_info_reference = json.load(infile2)
    return taxa_reference, order_info_reference

def process_metric(data):
    available_species, taxa_info, order_info, model_name, sequences, bin_function, bins, name_dict = data
    if len(available_species) >= 3:
        triads_species_names = pseudorandom_sequence_selection_between_order(available_species, taxa_info, order_info)
        if triads_species_names and triads_species_names['outgroup']:
            triat_data, _ = get_triads_info(triads_species_names, model_name, sequences, name_dict)
            bins = bin_function(triat_data, bins)
            return True, bins
    return False, bins

def process_gene_file(gene_file_path, model_name, output_dir, num_cpus, taxa_reference, order_reference):
    sequences = deserialise_object(json.load(open(gene_file_path, 'r')))
    species_names = extract_species_names(sequences)
    name_dict = {name.split('-', 1)[0]: name for name in sequences.names}
    taxa_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)

    file_name = os.path.basename(gene_file_path).rsplit('.', 1)[0]
    output_path = os.path.join(output_dir, file_name)

    if not os.path.exists(output_path):
        os.makedirs(output_path)

    jsd_bins = {'size': 0.0001, 'bins': initialize_bins(0.0001, 0.01)}
    ens_diff_bins = {'size': 0.0001, 'bins': initialize_bins(0.0001, 0.5)}
    shortest_ens_bins = {'size': 0.001, 'bins': initialize_bins(0.001, 0.1)}

    max_iterations = 100
    iteration_count = 0
    pool = Pool(processes=num_cpus)

    while iteration_count < max_iterations:
        process_data = [
            (species_names, taxa_info, order_info, model_name, sequences, triads_binning_jsd, jsd_bins, name_dict),
            (species_names, taxa_info, order_info, model_name, sequences, triads_binning_ens_diff, ens_diff_bins, name_dict),
            (species_names, taxa_info, order_info, model_name, sequences, triads_binning_shortest_ens, shortest_ens_bins, name_dict)
        ]

        results = pool.map(process_metric, process_data)

        process_any = False
        for result, data in zip(results, process_data):
            processed, bins = result
            _, _, _, _, _, _, current_bins, _ = data

            if processed:
                current_bins['bins'] = bins['bins']  # Update the bins
                process_any = True

        iteration_count += 1   
        print(iteration_count)

        if not process_any:  # Break the loop if no processing was possible for any metric
            break     

    pool.close()
    pool.join()

    with open(os.path.join(output_path, 'jsd_bins.json'), 'w') as f:
        json.dump(jsd_bins, f, indent=4)
    with open(os.path.join(output_path, 'ens_diff_bins.json'), 'w') as f:
        json.dump(ens_diff_bins, f, indent=4)
    with open(os.path.join(output_path, 'shortest_ens_bins.json'), 'w') as f:
        json.dump(shortest_ens_bins, f, indent=4)

def main_process(directory_path, model_name, output_dir, num_cpus):
    taxa_reference, order_reference = load_taxa_order_reference()

    gene_files = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.json')]

    for gene_file in gene_files:
        process_gene_file(gene_file, model_name, output_dir, num_cpus, taxa_reference, order_reference)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process multiple gene sequences.')
    parser.add_argument('-d', '--directory_path', type=str, required=True, help='Path to the directory containing gene files')
    parser.add_argument('-m', '--model_name', type=str, required=True, help='Name of the evolutionary model to use')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save output JSON files')
    parser.add_argument('-np', '--num_cpus', type=int, default=1, help='Number of CPUs to use for multiprocessing')

    args = parser.parse_args()

    main_process(args.directory_path, args.model_name, args.output_dir, args.num_cpus)




# def triads_selecting_binning(triad_data, jsd_bins, ens_diff_bins, shortest_ens_bins):
#     triad_added = False  # Flag to indicate if the triad was added to any bin
    
#     triad_values  = triad_data['triats_info']
#     # Calculate bin indices for each property
#     jsd_bin_index = int(triad_values['ingroup_jsd'] // jsd_bins['size'])
#     ens_diff_bin_index = int(triad_values['ens_difference'] // ens_diff_bins['size'])
#     shortest_ens_bin_index = int(triad_values['shortest_ingroup_ens'] // shortest_ens_bins['size'])

#     # Check and add to JSD bin
#     if jsd_bin_index in jsd_bins['bins'] and jsd_bins['bins'][jsd_bin_index] is None:
#         jsd_bins['bins'][jsd_bin_index] = triad_data
#         triad_added = True

#     # Check and add to ENS difference bin
#     if ens_diff_bin_index in ens_diff_bins['bins'] and ens_diff_bins['bins'][ens_diff_bin_index] is None:
#         ens_diff_bins['bins'][ens_diff_bin_index] = triad_data
#         triad_added = True

#     # Check and add to shortest ENS bin
#     if shortest_ens_bin_index in shortest_ens_bins['bins'] and shortest_ens_bins['bins'][shortest_ens_bin_index] is None:
#         shortest_ens_bins['bins'][shortest_ens_bin_index] = triad_data
#         triad_added = True

#     return jsd_bins, ens_diff_bins, shortest_ens_bins, triad_added


# def pseudorandom_sequence_selection_between_order(species_names,taxa_reference, order_reference):
#     ingroup_species_pair = random.sample(species_names, 2)

#     taxa_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, species_names)

#     species1_taxa = taxa_info[ingroup_species_pair[0]]
#     species2_taxa = taxa_info[ingroup_species_pair[1]]

#     if species1_taxa['order'] == species2_taxa['order']:
#         same_order = True
#         same_family = 'family' in species1_taxa and 'family' in species2_taxa and species1_taxa['family'] == species2_taxa['family']
#     else:
#         same_order = False
#         same_family = False

#     if same_family:
#         family_members = [s for s, taxa in taxa_info.items() if taxa.get('family') == species1_taxa['family']]
#         family_members = list(set(family_members) - set(ingroup_species_pair))
#         if family_members:
#             outgroup_species = random.choice(family_members)
#         else:
#             order_members = order_info[species1_taxa['order']]
#             order_members = list(set(order_members) - set(ingroup_species_pair))
#             outgroup_species = random.choice(order_members) if order_members else None
    
#     elif same_order:
#         order_members = order_info[species1_taxa['order']]
#         order_members = list(set(order_members) - set(ingroup_species_pair))
#         if order_members:
#             outgroup_species = random.choice(order_members)
#         else:
#             other_species = [s for s, taxa in taxa_info.items() if taxa['order'] not in [species1_taxa['order'], species2_taxa['order']]]
#             outgroup_species = random.choice(other_species) if other_species else None
#     else:
#         other_species = [s for s, taxa in taxa_info.items() if taxa['order'] not in [species1_taxa['order'], species2_taxa['order']]]
#         outgroup_species = random.choice(other_species) if other_species else None

#     triad_species_name = {
#         'ingroup1': ingroup_species_pair[0],
#         'ingroup2': ingroup_species_pair[1],
#         'outgroup': outgroup_species
#     }

#     return triad_species_name