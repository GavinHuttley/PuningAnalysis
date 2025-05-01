from clock_project.simulation.wts import calculate_non_stationarity
from cogent3.maths.measure import jsd
from cogent3.util.deserialise import deserialise_object
import multiprocessing
from cogent3.app.composable import NotCompleted
import json
import os
import glob
from cogent3 import get_app, open_data_store
import click

load_json_app = get_app("load_json")


def change_nucleotide_order(motif_prob):
    nuc_freq = {'A': motif_prob[0],'C': motif_prob[1],'G': motif_prob[2],'T': motif_prob[3]}
    nuc_freq_reordered = {'T': nuc_freq['T'], 'C': nuc_freq['C'], 'A': nuc_freq['A'], 'G': nuc_freq['G']}
    return list(nuc_freq_reordered.values())

def get_triples_info(model_fitting_result, triples_species_name):
    triples_aln = model_fitting_result.alignment
    ingroup_species_gene_name = [triples_species_name['ingroup1'], triples_species_name['ingroup2']]
    outgroup = triples_species_name['outgroup']
    matrices = {n: model_fitting_result.lf.get_rate_matrix_for_edge(n, calibrated=True) for n in triples_species_name.values()}
    sub_ens_tree = model_fitting_result.lf.get_ens_tree()
    lca_node = sub_ens_tree.lowest_common_ancestor([triples_species_name['ingroup1'], triples_species_name['ingroup2']]).get_node_names()[0]
    ingroup_seqs = triples_aln.take_seqs(ingroup_species_gene_name)
    nuc_freqs = ingroup_seqs.probs_per_seq()
    nuc_freqs1 = nuc_freqs[ingroup_species_gene_name[0]]
    nuc_freqs2 = nuc_freqs[ingroup_species_gene_name[1]]
    nuc_freqs1_reordered = change_nucleotide_order(nuc_freqs1)
    nuc_freqs2_reordered = change_nucleotide_order(nuc_freqs2)
    ingroup_jsd = jsd(nuc_freqs1_reordered, nuc_freqs2_reordered)
    distribution_internal_root = list(model_fitting_result.lf.get_motif_probs_by_node()[lca_node])
    distribution_root = list(model_fitting_result.lf.get_motif_probs_by_node()['root'])
    triples_species_name_list = list(triples_species_name.values())
    triples_species_name_list.extend([lca_node])
    ens_dict = {n: sub_ens_tree.to_rich_dict()['edge_attributes'][n]['length'] for n in triples_species_name_list}
    ens_difference = abs(ens_dict[ingroup_species_gene_name[0]]-ens_dict[ingroup_species_gene_name[1]])
    nabla_dict = {}
    for species in ingroup_species_gene_name:
        nabla_dict[species] = calculate_non_stationarity(distribution_internal_root, matrices[species], ens_dict[species])
    nabla_dict[outgroup] = calculate_non_stationarity(distribution_root, matrices[outgroup], ens_dict[outgroup])
    shortest_ens = min(ens_dict.values())  
    matrices_list = {n: matrices[n].to_array().tolist() for n in triples_species_name.values()}
    distributions = {'internal_node': distribution_internal_root, 'root': distribution_root, 'ingroup1': list(nuc_freqs1_reordered), 'ingroup2': list(nuc_freqs2_reordered)}
    triad_data = {'triples_species_names': triples_species_name, 'triples_info_small_tree': 
        {'ens_difference': ens_difference,
        'shortest_ingroup_ens': shortest_ens,
        'ens': ens_dict,
        'ingroup_jsd': ingroup_jsd,
        'nuc_freqs_dict': distributions,
        'nabla_values': nabla_dict,
        'matrices': matrices_list}}
    return triad_data

def process_path(path, triples_info_dir):
    triples_info_dict = {}
    file_name = os.path.basename(path.rstrip('/'))
    print(file_name)
    triples_info_path = os.path.join(triples_info_dir, file_name, "triples_species_names_dict.json")
    triples_info = deserialise_object(triples_info_path)
    model_fitting_result_path = os.path.join(path, "model_fitting_result")
    model_fitting_result_input_dstore = open_data_store(model_fitting_result_path, suffix='json')
    for sub_path in model_fitting_result_input_dstore:
        identifier = sub_path.unique_id.split('.')[0]
        triples_species_name = triples_info[identifier]
        model_fittqting_result = load_json_app(sub_path)
        if not isinstance(model_fittqting_result, NotCompleted):
            triple_info = get_triples_info(model_fittqting_result, triples_species_name)
            triples_info_dict[identifier] = triple_info
    with open(os.path.join(path, 'triples_info_dict_new.json'), 'w') as f:
        json.dump(triples_info_dict, f, indent=4)


@click.command()
@click.argument('model_fitting_result_dir', type=click.Path(exists=True))
@click.option('--triples_info_dir', '-td', help='Path to the directory containing info of all triples')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--limit', '-l', type=int, help='Limit the number of files to process')
def main(model_fitting_result_dir, triples_info_dir, num_processes, limit):
    """Fit model for sequences in different paths."""
    paths = glob.glob(os.path.join(model_fitting_result_dir, "*/"))

    if limit is not None:
        paths = paths[:limit]

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, triples_info_dir) for path in paths])

if __name__ == "__main__":
    main()

