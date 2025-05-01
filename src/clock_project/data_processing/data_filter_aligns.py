from cogent3 import get_app, open_data_store
import click
import json

from cogent3.app.typing import SeqsCollectionType
import multiprocessing
from cogent3.app.composable import define_app
from cogent3.app.composable import NotCompleted
# todo: add the code for wirte the taxanomic triple sequence into data store 


@define_app
def replace_common_species_names(alignment: SeqsCollectionType	
) -> SeqsCollectionType:
    with open ('/Users/gulugulu/clock/mammal_orthologs_hsap_1/available_species_name.json', 'r') as infile1:
        available_species_names = json.load(infile1)

    with open ('/Users/gulugulu/clock/mammal_orthologs_hsap_1/common_names.json', 'r') as infile2:
        common_names = json.load(infile2)
    name_map = dict(zip(available_species_names, common_names))
    
    new_names = {}
    
    # Iterate through the sequence names in the alignment
    for name in alignment.names:
        species_info, _ = name.split('-', 1)
        
        # Check if the species part exists in the name mapping
        if species_info in name_map:
            # Create a new name using the common name and the gene info
            new_name = name_map[species_info]
            new_names[name] = new_name
        else:
            # If no mapping exists, keep the original name
            print(species_info)
            new_names[name] = species_info
    
    # Update the sequence names in the alignment with new names
    return alignment.rename_seqs(lambda x: new_names[x])

common_name_renamer = replace_common_species_names()



@define_app
def filter_by_num_seqs(data: SeqsCollectionType, min_num_seqs: int) -> SeqsCollectionType:
    """
    Filters sequence collections based on the minimum number of sequences.

    """
    num_seqs = data.num_seqs

    if num_seqs < min_num_seqs:
        msg = f"{num_seqs} < min_num_seqs {min_num_seqs}"
        return NotCompleted("FALSE", None, msg, source=data)

    return data


min_sample_size = filter_by_num_seqs(min_num_seqs = 30)
third_pos = get_app("take_codon_positions", 3)
just_nucs =get_app("omit_degenerates", moltype="dna")
min_length = get_app("min_length", 550)

@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('-o', '--output_dir', help='Output directory for alignment.fa files')
@click.option('-n', '--num_processes', default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process')

def main(input, num_processes, output_dir, limit):
    """
    Fit model for sequences in different paths.
    
    """
    input_data_store = open_data_store(input, suffix='fa', limit=limit)
    out_dstore = open_data_store(output_dir, mode="w", suffix="fa")
    loader = get_app("load_aligned", format="fasta", moltype="dna")
    writer = get_app("write_seqs", out_dstore, format="fasta")
    aln_processer = loader + common_name_renamer + third_pos + just_nucs + min_length + min_sample_size + writer

    aln_processer.apply_to(input_data_store.completed, parallel=True, par_kw=dict(max_workers=num_processes))

if __name__ == "__main__":
    main()


            
