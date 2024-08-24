
import multiprocessing
import json
import os
import glob
from cogent3 import get_app, open_data_store
import click
from clock_project.simulation.wts import calculate_non_stationarity
from cogent3.maths.measure import jsd
from cogent3.app.composable import define_app


from cogent3.app.typing import AlignedSeqsType, SerialisableType, IdentifierType

load_json_app = get_app("load_json")

@define_app
def model_fitting(aln: AlignedSeqsType, opt_args = None) -> SerialisableType:
    tree="('outgroup_edge3',('ingroup_edge1','ingroup_edge2')'internal_node');"
    GN_model = get_app("model", "GN", 
                    tree=tree,
                    opt_args=opt_args,
                    # unique_trees=True,
                    # lf_args=dict(discrete_edges=[['outgroup_edge3']]),
                    time_het="max",
                    optimise_motif_probs=True)
    res = GN_model(aln)
    return res

model_fitting_app = model_fitting()

@define_app
def customised_load_json(DataMember: IdentifierType) -> AlignedSeqsType:
    unique_id = DataMember.unique_id
    aln = load_json_app(DataMember)
    aln.source = unique_id
    return aln

def get_id(result):
    return result.alignment.source

load_json_app_customised = customised_load_json()

def get_triads_info(path):
    aln = load_json_app_customised(path)
    model_fitting_result = model_fitting_app(aln)
    seq_names = aln.names
    matrices = {n: model_fitting_result.lf.get_rate_matrix_for_edge(n, calibrated=True) for n in seq_names}
    sub_ens_tree = model_fitting_result.lf.get_ens_tree()
    lca_node = sub_ens_tree.lowest_common_ancestor(['ingroup_edge1', 'ingroup_edge2']).get_node_names()[0]
    ingroup_seqs = aln.take_seqs(['ingroup_edge1', 'ingroup_edge2'])
    nuc_freqs = ingroup_seqs.probs_per_seq()
    nuc_freqs1 = nuc_freqs['ingroup_edge1']
    nuc_freqs2 = nuc_freqs['ingroup_edge2']
    ingroup_jsd = jsd(nuc_freqs1, nuc_freqs2)
    distribution_internal_root = list(model_fitting_result.lf.get_motif_probs_by_node()[lca_node])
    jsd1 = jsd(nuc_freqs1, distribution_internal_root)
    jsd2 = jsd(nuc_freqs2, distribution_internal_root)
    distribution_root = list(model_fitting_result.lf.get_motif_probs_by_node()['root'])
    ens_dict = {n: sub_ens_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in seq_names}
    ens_difference = abs(ens_dict['ingroup_edge1']-ens_dict['ingroup_edge2'])
    jsd_dict = {'Ingroup_JSD': ingroup_jsd, 'JSD_difference': {'ingroup_edge1': jsd1, 'ingroup_edge2': jsd2}}
    nabla_dict = {}
    for species in ['ingroup_edge1', 'ingroup_edge2']:
        nabla_dict[species] = calculate_non_stationarity(distribution_internal_root, matrices[species], ens_dict[species])
    nabla_dict['outgroup_edge3'] = calculate_non_stationarity(distribution_root, matrices['outgroup_edge3'], ens_dict['outgroup_edge3'])
    shortest_ens = min(ens_dict.values())  
    matrices_list = {n: matrices[n].to_array().tolist() for n in seq_names}
    distributions = {'internal_node': distribution_internal_root, 'root': distribution_root, 'ingroup1': list(nuc_freqs1), 'ingroup2': list(nuc_freqs2)}
    triad_data = {'triads_species_names': seq_names, 'triads_info_small_tree': 
        {'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'jsd_dict': jsd_dict,
        'nuc_freqs_dict': distributions,
        'nabla_values': nabla_dict,
        'matrices': matrices_list}}
    return triad_data, model_fitting_result


def process_path(path, output_dir):
    file_name = os.path.basename(path.rstrip('/'))
    print(file_name)
    triads_alignment_input_data_store = open_data_store(path, mode="w", suffix="json")

    file_output_dir = os.path.join(output_dir, file_name)
    os.makedirs(file_output_dir, exist_ok=True)  
    result_dir = os.path.join(file_output_dir, 'model_fitting_result')
    os.makedirs(result_dir, exist_ok=True)  
    
    triads_info_dict = {}
    for datamember in triads_alignment_input_data_store:
        print(datamember)
        triads_info, model_fitting_result = get_triads_info(datamember)
        identifier_name = model_fitting_result.alignment.source
        identifier_number = identifier_name.split('.')[0]
        triads_info_dict[identifier_number] = triads_info
        out_dstore = open_data_store(result_dir, mode="w", suffix="json")
        write_json_app = get_app("write_json", data_store=out_dstore, id_from_source = get_id)
        write_json_app(model_fitting_result)
    
    with open(os.path.join(file_output_dir, 'triads_info_dict.json'), 'w') as f:
        json.dump(triads_info_dict, f, indent=4)

@click.command()
@click.argument('simulated_alignment_dir', type=click.Path(exists=True))
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory')
@click.option('--limit', '-l', type=int, help='Limit the number of files to process')
def main(simulated_alignment_dir, num_processes, output_dir, limit):
    """Fit model for sequences in different paths."""
    paths = glob.glob(os.path.join(simulated_alignment_dir, '*/'))

    if limit is not None:
        paths = paths[:limit]

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir) for path in paths])

if __name__ == "__main__":
    main()