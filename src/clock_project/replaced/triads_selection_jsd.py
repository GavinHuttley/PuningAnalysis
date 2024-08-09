import argparse
import os
import json
import random
from cogent3.util.deserialise import deserialise_object

from cogent3.maths.measure import jsd
from multiprocessing import Pool

from clock_project.genome_analysis.traids_selection import (
    extract_taxonomic_order_info, 
    extract_species_names, 
    get_triads_info, 
    get_sequence_alignment, 
    initialize_bins, 
    load_taxa_order_reference, 
    match_species_gene_name)


def get_sequence(ingroup_species_name, sequences):
    species_gene_names = []
    for species_name in ingroup_species_name.values():
        species_gene_name = match_species_gene_name(species_name, sequences)
        species_gene_names.append(species_gene_name)
    trait_sequences = sequences.take_seqs(species_gene_names)
    trait_sequences_no_stop_codon = trait_sequences.trim_stop_codons(strict = True)

    return trait_sequences_no_stop_codon

def pseudorandom_ingroup_selection_between_order_jsd(available_species, taxa_reference, order_reference, sequences):
    if len(available_species) < 2:
        return None, None
    ingroup_species_pair = random.sample(available_species, 2)
    ingroup_species_name = {
        'ingroup1': ingroup_species_pair[0],
        'ingroup2': ingroup_species_pair[1]
    }
    ingroup_sequence = get_sequence(ingroup_species_name, sequences)
    ingroup_gene_name = ingroup_sequence.names
    nuc_freq = ingroup_sequence.probs_per_seq()
    jsd_value = jsd(nuc_freq[ingroup_gene_name[0]], nuc_freq[ingroup_gene_name[1]])
    return ingroup_species_name, jsd_value



def triads_binning_jsd_check(ingroup_species_name, jsd_value, jsd_bins):
    jsd_bin_index = int(jsd_value // jsd_bins['size'])
    assigned = False
    if jsd_bin_index in jsd_bins['bins'] and jsd_bins['bins'][jsd_bin_index] is None:
        jsd_bins['bins'][jsd_bin_index] = ingroup_species_name
        assigned = True
    return assigned, jsd_bin_index

def get_outgroup_species(ingroup_species_name, taxa_reference, order_reference, available_species):
    taxa_info, order_info = extract_taxonomic_order_info(taxa_reference, order_reference, available_species)

    species1_taxa = taxa_info[ingroup_species_name['ingroup1']]
    species2_taxa = taxa_info[ingroup_species_name['ingroup2']]

    if species1_taxa['order'] == species2_taxa['order']:
        same_order = True
        same_family = 'family' in species1_taxa and 'family' in species2_taxa and species1_taxa['family'] == species2_taxa['family']
    else:
        same_order = False
        same_family = False

    if same_family:
        family_members = [s for s, taxa in taxa_info.items() if taxa.get('family') == species1_taxa['family']]
        family_members = list(set(family_members) - set(ingroup_species_name))
        if family_members:
            outgroup_species = random.choice(family_members)
        else:
            order_members = order_info[species1_taxa['order']]
            order_members = list(set(order_members) - set(ingroup_species_name))
            outgroup_species = random.choice(order_members) if order_members else None
    
    elif same_order:
        order_members = order_info[species1_taxa['order']]
        order_members = list(set(order_members) - set(ingroup_species_name))
        if order_members:
            outgroup_species = random.choice(order_members)
        else:
            other_species = [s for s, taxa in taxa_info.items() if taxa['order'] not in [species1_taxa['order'], species2_taxa['order']]]
            outgroup_species = random.choice(other_species) if other_species else None
    else:
        other_species = [s for s, taxa in taxa_info.items() if taxa['order'] not in [species1_taxa['order'], species2_taxa['order']]]
        outgroup_species = random.choice(other_species) if other_species else None

    triad_species_name = {
        'ingroup1': ingroup_species_name['ingroup1'],
        'ingroup2': ingroup_species_name['ingroup2'],
        'outgroup': outgroup_species
    }

    return triad_species_name



def process_metric(data):
    available_species, taxa_reference, order_reference, model_name, sequences, jsd_bins, used_species_jsd = data
    if not available_species:
        return False, jsd_bins, used_species_jsd
    ingroup_species_name, jsd_value = pseudorandom_ingroup_selection_between_order_jsd(available_species, taxa_reference, order_reference, sequences)
    if ingroup_species_name:
        assigned, jsd_bin_index = triads_binning_jsd_check(ingroup_species_name, jsd_value, jsd_bins)
        if assigned:
            used_species_jsd.update(ingroup_species_name.values())
            triads_species_names = get_outgroup_species(ingroup_species_name, taxa_reference, order_reference, available_species)
            if triads_species_names and triads_species_names['outgroup']:
                processed = True
                triads_data, model_fitting_result = get_triads_info(triads_species_names, model_name, sequences)
                jsd_bins['bins'][jsd_bin_index] = triads_data
                return processed, jsd_bins, used_species_jsd
    return False, jsd_bins, used_species_jsd

def process_gene_file(sequence_dir, model_name, output_dir, num_cpus):
    taxa_reference, order_reference = load_taxa_order_reference()
    jsd_bins = {'size': 0.000001, 'bins': initialize_bins(0.00001, 0.05)}
    used_species_jsd = set()
    sequences = deserialise_object(json.load(open(sequence_dir, 'r')))
    species_names = extract_species_names(sequences)
    file_name = os.path.basename(sequence_dir).rsplit('.', 1)[0]
    print(file_name)

    output_path = os.path.join(output_dir, file_name)  # Full output directory path
    if not os.path.exists(output_path):
        os.makedirs(output_path)  # Create the directory if it does not exist
            
    max_iterations = 1000  # Set a reasonable upper limit for loop iterations
    iteration_count = 0
    pool = Pool(processes=num_cpus)  # Use specified number of CPUs

    while iteration_count < max_iterations and len(species_names) - len(used_species_jsd) >= 2:
        available_species = list(set(species_names) - used_species_jsd)
        process_data = (available_species, taxa_reference, order_reference, model_name, sequences, jsd_bins, used_species_jsd)
        processed, jsd_bins, used_species_jsd = pool.apply(process_metric, (process_data,))
        iteration_count += 1
        if not processed:
            break

    pool.close()
    pool.join()

    # Saving the bins data to JSON files
    with open(os.path.join(output_path, 'jsd_bins.json'), 'w') as f:
        json.dump(jsd_bins, f, indent=4)

def main_process(directory_path, model_name, output_dir, num_cpus):
    sequence_dirs = [os.path.join(directory_path, f) for f in os.listdir(directory_path) if f.endswith('.json')]

    for sequence_dir in sequence_dirs:
        process_gene_file(sequence_dir, model_name, output_dir, num_cpus)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Process multiple gene sequences for JSD.')
    parser.add_argument('-d', '--directory_path', type=str, required=True, help='Path to the directory containing gene files')
    parser.add_argument('-m', '--model_name', type=str, required=True, help='Name of the evolutionary model to use')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Directory to save output JSON files')
    parser.add_argument('-np', '--num_cpus', type=int, default=1, help='Number of CPUs to use for multiprocessing')
    args = parser.parse_args()

    main_process(args.directory_path, args.model_name, args.output_dir, args.num_cpus)