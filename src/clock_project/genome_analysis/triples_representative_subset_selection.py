import os
import json
import pandas as pd
import numpy as np
from cogent3 import get_app, open_data_store


load_json_app = get_app("load_json")

def read_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)


def sample_triples(info_path, align_path, sample_fraction, write_json_app):
    # Calculate the number of triples in the current gene
    current_gene_triples = len([name for name in os.listdir(align_path) if name.endswith('.json')])
    sampling_interval = int(1 / sample_fraction)
    print(sampling_interval)
    print(current_gene_triples)

    info_data = read_json(info_path)
    print(info_path)
    triples_data = []
    
    for key, value in info_data.items():
        if value is not None:
            triples_data.append({
                'identifier': key,
                'ens_difference': value['triples_info_small_tree']['ens_difference'],
                'triple_species_name': value['triples_species_names']
            })
        
    df = pd.DataFrame(triples_data)
    df['ens_difference'] = pd.to_numeric(df['ens_difference'])  # Ensure it's numeric for sorting
    
    # Sort by ens_difference
    df = df.sort_values(by='ens_difference', ascending=True).reset_index(drop=True)

    # Pick every nth triple based on the sampling interval
    sampled_df = df.iloc[::sampling_interval].reset_index(drop=True)


    # Copy corresponding alignment files to the output directory
    for _, row in sampled_df.iterrows():
        print(sampled_df.shape[0])
        identifier = row['identifier']
        triples_species_name = row['triple_species_name']
        gene_name = os.path.basename(os.path.dirname(info_path))
        src_file = os.path.join(align_path, f"{identifier}.json")
        aln = load_json_app(src_file)
        aln.info['triples_species_name'] = triples_species_name
        file_identifier = f"{gene_name}_{identifier}.json"
        write_json_app(aln, identifier=file_identifier)

    

def main():
    # Paths for the directories
    info_dir = '/Users/gulugulu/Desktop/honours/data_local_2/triples_model_fitting_550_threshold'
    align_dir = '/Users/gulugulu/Desktop/honours/data_local_2/triples_alignments_550_threshold'
    output_dir = '/Users/gulugulu/Desktop/honours/data_local_2/triples_aln_representative_subset'
    sample_fraction = 0.1

    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)

    # Iterate over all genes and sample triples
    for gene_name in os.listdir(info_dir):
        info_path = os.path.join(info_dir, gene_name, 'triples_info_dict.json')
        align_path = os.path.join(align_dir, gene_name)
        if os.path.exists(info_path) and os.path.isdir(align_path):
            sample_triples(info_path, align_path, sample_fraction, write_json_app)
            


if __name__ == "__main__":
    main()