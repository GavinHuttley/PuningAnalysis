from cogent3.util.deserialise import deserialise_object
import multiprocessing
import json
import os
import glob
from cogent3 import get_app, open_data_store
import click
from clock_project.simulation.wts import calculate_non_stationarity
import numpy as np
from cogent3.maths.measure import jsd


load_json_app = get_app("load_json")


def set_model(triads_species_name, ens_tree, model_name):
    tree_topology_triads = ens_tree.get_sub_tree(triads_species_name.values())
    model_triads = get_app("model", model_name, tree = tree_topology_triads, time_het="max", optimise_motif_probs=True)
    return model_triads

def get_triads_info(triads_species_name, triads_aln, model_name, ens_tree):
    ingroup_species_gene_name = [triads_species_name['ingroup1'], triads_species_name['ingroup2']]
    outgroup = triads_species_name['outgroup']
    model = set_model(triads_species_name, ens_tree, model_name)
    model_fitting_result = model(triads_aln)
    matrices = {n: model_fitting_result.lf.get_rate_matrix_for_edge(n, calibrated=True) for n in triads_species_name.values()}
    sub_ens_tree = model_fitting_result.lf.get_ens_tree()
    lca_node = sub_ens_tree.lowest_common_ancestor([triads_species_name['ingroup1'], triads_species_name['ingroup2']]).get_node_names()[0]
    ingroup_seqs = triads_aln.take_seqs(ingroup_species_gene_name)
    nuc_freqs = ingroup_seqs.probs_per_seq()
    nuc_freqs1 = nuc_freqs[ingroup_species_gene_name[0]]
    nuc_freqs2 = nuc_freqs[ingroup_species_gene_name[1]]
    ingroup_jsd = jsd(nuc_freqs1, nuc_freqs2)
    distribution_internal_root = list(model_fitting_result.lf.get_motif_probs_by_node()[lca_node])
    jsd1 = jsd(nuc_freqs1, distribution_internal_root)
    jsd2 = jsd(nuc_freqs2, distribution_internal_root)
    distribution_root = list(model_fitting_result.lf.get_motif_probs_by_node()['root'])
    ens_dict = {n: sub_ens_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triads_species_name.values()}
    ens_difference = abs(ens_dict[ingroup_species_gene_name[0]]-ens_dict[ingroup_species_gene_name[1]])
    jsd_dict = {'Ingroup_JSD': ingroup_jsd, 'JSD_difference': {ingroup_species_gene_name[0]: jsd1, ingroup_species_gene_name[1]: jsd2}}
    nabla_dict = {}
    for species in ingroup_species_gene_name:
        nabla_dict[species] = calculate_non_stationarity(distribution_internal_root, matrices[species], ens_dict[species])
    nabla_dict[outgroup] = calculate_non_stationarity(distribution_root, matrices[outgroup], ens_dict[outgroup])
    shortest_ens = min(ens_dict.values())  
    matrices_list = {n: matrices[n].to_array().tolist() for n in triads_species_name.values()}
    distributions = {'internal_node': distribution_internal_root, 'root': distribution_root, 'ingroup1': list(nuc_freqs1), 'ingroup2': list(nuc_freqs2)}
    triad_data = {'triads_species_names': triads_species_name, 'triads_info_small_tree': 
        {'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'jsd_dict': jsd_dict,
        'nuc_freqs_dict': distributions,
        'nabla_values': nabla_dict,
        'matrices': matrices_list}}
    return triad_data, model_fitting_result.lf


def process_path(path, model_name, result_lf_dir, triads_info_dir, output_dir):
    file_name = os.path.basename(path.rstrip('/'))
    print(file_name)
    triads_alignment_paths = glob.glob(os.path.join(path, "*.json"))
    triads_info_path = os.path.join(triads_info_dir, file_name, "triads_species_names_dict.json")
    triads_speccies_infos = json.load((open(triads_info_path, 'r')))
    result_lf_path = os.path.join(result_lf_dir, f"{file_name}.json")
    result_lf = load_json_app(result_lf_path)

    file_output_dir = os.path.join(output_dir, file_name)
    os.makedirs(file_output_dir, exist_ok=True)  
    result_dir = os.path.join(file_output_dir, 'model_fitting_result')
    os.makedirs(result_dir, exist_ok=True)  
    
    triads_info_dict = {}
    for alignment_path in triads_alignment_paths:
        identifier_number = os.path.basename(alignment_path).rsplit('.', 1)[0]
        triads_aln = deserialise_object(json.load(open(alignment_path, 'r')))
        triads_species_name = triads_speccies_infos[identifier_number]
        print(triads_species_name)
        ens_tree = result_lf.get_ens_tree()
        triads_info, model_fitting_result = get_triads_info(triads_species_name, triads_aln, model_name, ens_tree)
        triads_info_dict[identifier_number] = triads_info
        out_dstore = open_data_store(result_dir, mode="w", suffix="json")
        file_path = identifier_number + ".json"
        write_json_app = get_app("write_json", data_store=out_dstore)
        write_json_app(model_fitting_result, identifier=file_path)
    
    with open(os.path.join(file_output_dir, 'triads_info_dict.json'), 'w') as f:
        json.dump(triads_info_dict, f, indent=4)

@click.command()
@click.argument('triads_alignment_dir', type=click.Path(exists=True))
@click.option('--triads_info_dir', '-td', help='Path to the directory containing model fitting result')
@click.option('--result_lf_dir', '-rd', help='Path to the directory containing model fitting result')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--model_name', '-m', default='GN', type=click.Path(), help='model name')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory')
@click.option('--limit', '-l', type=int, help='Limit the number of files to process')
def main(triads_alignment_dir, triads_info_dir, result_lf_dir, num_processes, model_name, output_dir, limit):
    """Fit model for sequences in different paths."""
    paths = glob.glob(os.path.join(triads_alignment_dir, '*/'))
    

    if limit is not None:
        paths = paths[:limit]

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, model_name, result_lf_dir, triads_info_dir, output_dir) for path in sub_paths])

if __name__ == "__main__":
    main()

