import multiprocessing
import json
import os
import glob
from cogent3 import get_app, open_data_store
import click



load_json_app = get_app("load_json")
loader_aligned = get_app("load_aligned", format="fasta", moltype="dna")


def triples_model_fitting(triples_aln, triples_species_name, ens_tree):
    tree_topology_triples = ens_tree.get_sub_tree(triples_species_name.values())
    model_triples = get_app("model", 'GN', tree = tree_topology_triples, time_het="max", optimise_motif_probs=True)
    model_fitting_result = model_triples(triples_aln)
    return model_fitting_result


def process_path(path, result_lf_dir, triples_info_dir, output_dir):
    file_name = os.path.basename(path.rstrip('/'))
    print(file_name)
    triples_alignment_paths = open_data_store(path, mode="r", suffix="fa")
    triples_info_path = os.path.join(triples_info_dir, file_name, "triples_species_names_dict.json")
    triples_speccies_infos = json.load((open(triples_info_path, 'r')))
    result_lf_path = os.path.join(result_lf_dir, f"{file_name}.json")
    result_lf = load_json_app(result_lf_path).lf

    file_output_dir = os.path.join(output_dir, file_name)
    os.makedirs(file_output_dir, exist_ok=True)  
    result_dir = os.path.join(file_output_dir, 'model_fitting_result')
    os.makedirs(result_dir, exist_ok=True)  
    
    for alignment_path in triples_alignment_paths:
        identifier_number = alignment_path.unique_id.split('.')[0]
        triples_aln = loader_aligned(alignment_path)
        triples_species_name = triples_speccies_infos[identifier_number]
        ens_tree = result_lf.get_ens_tree()
        model_fitting_result = triples_model_fitting(triples_aln, triples_species_name, ens_tree)
        out_dstore = open_data_store(result_dir, mode="w", suffix="json")
        write_json_app = get_app("write_json", data_store=out_dstore)
        write_json_app(model_fitting_result, identifier= f'{identifier_number}.json')

    

@click.command()
@click.argument('triples_alignment_dir', type=click.Path(exists=True))
@click.option('--triples_info_dir', '-td', help='Path to the directory containing model fitting result')
@click.option('--result_lf_dir', '-rd', help='Path to the directory containing model fitting result')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory for model fitting result of taxanomic triples')
@click.option('--limit', '-l', type=int, help='Limit the number of files to process')
def main(triples_alignment_dir, triples_info_dir, result_lf_dir, num_processes, output_dir, limit):
    """Fit model for sequences in different paths."""
    paths = glob.glob(os.path.join(triples_alignment_dir, "*/"))

    if limit is not None:
        paths = paths[:limit]

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, result_lf_dir, triples_info_dir, output_dir) for path in paths])

if __name__ == "__main__":
    main()

