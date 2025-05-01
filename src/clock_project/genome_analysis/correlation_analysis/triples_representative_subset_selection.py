import os
import json
import pandas as pd
from cogent3 import get_app, open_data_store


load_json_app = get_app("load_json")
loader_aligned = get_app("load_aligned", format="fasta", moltype="dna")


def read_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)

def sample_triples(info_path, align_path, sample_fraction, writer):
    sampling_interval = int(1 / sample_fraction)

    info_data = read_json(info_path)
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
    
    df = df.sort_values(by='ens_difference', ascending=True).reset_index(drop=True)

    sampled_df = df.iloc[::sampling_interval].reset_index(drop=True)


    for _, row in sampled_df.iterrows():
        identifier = row['identifier']
        triples_species_name = row['triple_species_name']
        gene_name = os.path.basename(os.path.dirname(info_path))
        src_file = os.path.join(align_path, f"{identifier}.fa")
        aln = loader_aligned(src_file)
        aln.info['triples_species_name'] = triples_species_name
        writer(aln, identifier = f"{gene_name}_{identifier}.fa")
            

def main():
    info_dir = '/Users/gulugulu/clock/mammal_orthologs_hsap_1/triples_model_fitting'
    align_dir = '/Users/gulugulu/clock/mammal_orthologs_hsap_1/taxanomic_triples_alignments'
    output_dir = '/Users/gulugulu/clock/mammal_orthologs_hsap_1/triples_representative_subset'
    sample_fraction = 0.1

    out_dstore = open_data_store(output_dir, mode="w", suffix="fa")
    writer = get_app("write_seqs", data_store=out_dstore, format="fasta")


    for gene_name in os.listdir(info_dir):
        info_path = os.path.join(info_dir, gene_name, 'triples_info_dict_new.json')
        align_path = os.path.join(align_dir, gene_name)
        if os.path.exists(info_path) and os.path.isdir(align_path):
            sample_triples(info_path, align_path, sample_fraction, writer)

if __name__ == "__main__":
    main()