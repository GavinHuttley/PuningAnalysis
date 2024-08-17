import os
import json
import pandas as pd
import numpy as np

def read_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)


def sample_triples(info_path, align_path, output_path, sample_fraction):
    # Calculate the number of triples in the current gene
    current_gene_triples = len([name for name in os.listdir(align_path) if name.endswith('.json')])
    sampling_interval = int(1 / sample_fraction)
    print(sampling_interval)
    print(current_gene_triples)

    info_data = read_json(info_path)
    print(info_path)
    triples_data = []
    all_triples_info = {}
    
    for key, value in info_data.items():
        if value is not None:
            triples_data.append({
                'identifier': key,
                'ens_difference': value['triads_info_small_tree']['ens_difference'],
                'traid_species_name': value['triads_species_names']
            })
        
    df = pd.DataFrame(triples_data)
    df['ens_difference'] = pd.to_numeric(df['ens_difference'])  # Ensure it's numeric for sorting
    
    # Sort by ens_difference
    df = df.sort_values(by='ens_difference', ascending=True).reset_index(drop=True)

    # Pick every nth triple based on the sampling interval
    sampled_df = df.iloc[::sampling_interval].reset_index(drop=True)

    # Create output directory
    os.makedirs(output_path, exist_ok=True)

    # Copy corresponding alignment files to the output directory
    all_triples_info = {}
    for idx, row in sampled_df.iterrows():
        print(sampled_df.shape[0])
        identifier = row['identifier']
        traid_species_name = row['traid_species_name']
        all_triples_info[identifier] = traid_species_name
        genename = os.path.basename(os.path.dirname(info_path))
        src_file = os.path.join(align_path, f"{identifier}.json")
        dst_file = os.path.join(output_path, f"{genename}_{identifier}.json")
        if os.path.exists(src_file):
            os.system(f"cp {src_file} {dst_file}")
    return all_triples_info
    

# Paths for the directories
info_dir = '/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_model_fitting_350_threshold'
align_dir = '/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_alignment_350_threshold'
output_dir = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset'
output_di2 = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset_info.json'
sample_fraction = 0.1

# Iterate over all genes
all_triples_info_dict = {}
for gene_name in os.listdir(info_dir):
    info_path = os.path.join(info_dir, gene_name, 'triads_info_dict.json')
    align_path = os.path.join(align_dir, gene_name)
    output_path = output_dir
    output_path2 = output_di2


    if os.path.exists(info_path) and os.path.isdir(align_path):
        all_triples_info = sample_triples(info_path, align_path, output_path, sample_fraction)
        all_triples_info_dict[gene_name] = all_triples_info 
        
with open (output_path2, 'w') as outfile:
    json.dump(all_triples_info_dict, outfile, indent=4)


