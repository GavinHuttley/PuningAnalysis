from cogent3 import get_app, open_data_store
import click

from clock_project.data_processing.data_filter_seqs import cpos3
import json

from cogent3 import get_app
from cogent3.app.typing import SeqsCollectionType
import multiprocessing
from cogent3.app.composable import define_app


alignment_lengths = {}
sample_size_dict = {}
small_sample_size = {}
short_alignment = {}
valid_homology_list = []

load_json_app = get_app("load_json")

@define_app
def replace_common_species_names(alignment: SeqsCollectionType	
) -> SeqsCollectionType:
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/genome_information/available_species_name.json', 'r') as infile1:
        available_species_names = json.load(infile1)

    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/genome_information/common_names.json', 'r') as infile2:
        common_names = json.load(infile2)
    name_map = dict(zip(available_species_names, common_names))
    
    new_names = {}
    
    # Iterate through the sequence names in the alignment
    for name in alignment.names:
        species_info, gene_info = name.split('-', 1)
        
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

def filter_alignment(path, write_json_app_full_aln, write_json_app_filtered):
    file_name = path.unique_id.split('.')[0]
    aln = load_json_app(path)
    renamed_aln = common_name_renamer(aln)
    aln_3rd = cpos3(renamed_aln)
    aln_3rd_no_degenerates = aln_3rd.no_degenerates()
    try:
        alignment_length = aln_3rd_no_degenerates.get_lengths()[0]
        print(alignment_length)
        sample_size = aln_3rd_no_degenerates.num_seqs
        print(sample_size)
        alignment_lengths[file_name] = alignment_length
        sample_size_dict[file_name] = sample_size
        if sample_size >= 30 and alignment_length >= 550:
            renamed_full_aln = common_name_renamer(aln)
            write_json_app_full_aln(renamed_full_aln, identifier=f'{file_name}.json')
        else:
            renamed_full_aln = common_name_renamer(aln)
            write_json_app_filtered(renamed_full_aln, identifier=f'{file_name}.json')
    except Exception as e:
        renamed_full_aln = common_name_renamer(aln)
        write_json_app_filtered(renamed_full_aln, identifier=f'{file_name}.json')
        
        



@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('-o1', '--output_dir_full_aln', help='Output directory for JSON files')
@click.option('-o2', '--output_dir_filtered', help='Output directory for JSON files')
@click.option('-n', '--num_processes', default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process')
def main(input, output_dir_full_aln, output_dir_filtered, num_processes, limit):
    """
    Fit model for sequences in different paths.
    
    INPUT: Path to the directory containing sequence files or a JSON file listing the paths.
    """
    input_data_store = open_data_store(input, suffix='json', limit=limit)
    out_dstore_full = open_data_store(output_dir_full_aln, mode="w", suffix="json")
    out_dstore_filtered = open_data_store(output_dir_filtered, mode="w", suffix="json")
    write_json_app_full_aln = get_app("write_json", data_store=out_dstore_full)
    write_json_app_filtered = get_app("write_json", data_store=out_dstore_filtered)
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(filter_alignment, [(path, write_json_app_full_aln, write_json_app_filtered) for path in input_data_store.completed])

    # with open('/Users/gulugulu/clock/mammal_orthologs_hsap_1/sampled_aligned/alignment_lengths.json', 'w') as aln_len_file:
    #     json.dump(alignment_lengths, aln_len_file, indent=4)
    
    # with open('/Users/gulugulu/clock/mammal_orthologs_hsap_1/sampled_aligned/sample_sizes.json', 'w') as sample_size_file:
    #     json.dump(sample_size_dict, sample_size_file, indent=4)

if __name__ == "__main__":
    main()


            
