import os
import json
import pandas as pd
import numpy as np


def bin_based_sampling(df, sample_size):
    # Define bins edges
    df['bin'] = pd.cut(df['ens_difference'], bins=sample_size, labels=False)

    # Group by bin and sample one item per bin
    sampled_df = df.groupby('bin').apply(lambda x: x.sample(1), include_groups=False).reset_index(drop=True)

    # Drop the bin column as it's no longer needed

    return sampled_df

def read_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)


def sample_triples(info_path, align_path, output_path, total_triples):
    # Calculate the number of triples in the current gene
    current_gene_triples = len([name for name in os.listdir(align_path) if name.endswith('.json')])
    print(current_gene_triples)
    # Calculate the fraction of triples to sample
    sample_fraction = current_gene_triples / total_triples
    # Calculate how many samples to take
    num_samples = int(np.round(sample_fraction * total_triples * 0.125))
    print(num_samples)
    info_data = read_json(info_path)
    triples_data = []
    
    for key, value in info_data['bins'].items():
        if value is not None:
            triples_data.append({
                'identifier': key,
                'ens_difference': value['triads_info_big_tree']['ens_difference']
            })

    # Convert to DataFrame for easier handling
    df = pd.DataFrame(triples_data)
    df['ens_difference'] = pd.to_numeric(df['ens_difference'])  # Ensure it's numeric for sorting
    
    sampled_df = bin_based_sampling(df, num_samples)

    # Create output directory
    os.makedirs(output_path, exist_ok=True)

    # Copy corresponding alignment files to the output directory
    for idx, row in sampled_df.iterrows():
        print(sampled_df.shape[0])
        identifier = row['identifier']
        genename = os.path.basename(os.path.dirname(info_path))
        src_file = os.path.join(align_path, f"{identifier}.json")
        dst_file = os.path.join(output_path, f"{genename}_{identifier}.json")
        if os.path.exists(src_file):
            os.system(f"cp {src_file} {dst_file}")
    return current_gene_triples

# Paths for the directories
info_dir = '/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_350_threshold'
align_dir = '/Users/gulugulu/Desktop/honours/data_local/whole_genome_mammal87/triads_alignment_350_threshold'
output_dir = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset'
total_triples = 21176  # Total number of triples

# Iterate over all genes
all_gene_triples = 0
for gene_name in os.listdir(info_dir):
    info_path = os.path.join(info_dir, gene_name, 'ens_diff_bins.json')
    align_path = os.path.join(align_dir, gene_name)
    output_path = output_dir

    if os.path.exists(info_path) and os.path.isdir(align_path):
        current_gene_triples = sample_triples(info_path, align_path, output_path, total_triples)
        all_gene_triples += current_gene_triples

print(all_gene_triples)