from ensembl_lite._species import Species
from cogent3 import load_tree
from cogent3.util.deserialise import deserialise_object
import argparse
import multiprocessing
import json
import os
import glob
from cogent3 import get_app, open_data_store

tree_path = '/Users/gulugulu/repos/PuningAnalysis/data/dataset2_ensemble_trees/raw_data/tree_nh_file/vertebrates_species-tree_Ensembl.nh'
tree_vertebrates= load_tree(tree_path, format=None, underscore_unmunge=False)
outgroup_species = ['aquila_chrysaetos_chrysaetos']
from ensembl_lite._species import Species

def make_db_prefixes():
    return [n.lower().replace(" ", "_") for n in Species.get_species_names()]

db_prefixes = make_db_prefixes()

def find_db_prefix(name):
    for db_prefix in db_prefixes:
        if name.startswith(db_prefix):
            return db_prefix
    return None

import json
with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/available_species_name.json', 'r') as infile:
    available_species_names = json.load(infile)

with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/common_names_mapping.json',  'r') as common_name_infile:
    common_names_mapping = json.load(common_name_infile)

db_names = []
probs = []
for tip_name in tree_vertebrates.get_tip_names():
    db_name = find_db_prefix(tip_name.lower())
    if db_name is None:
        probs.append(tip_name)
    else:
        db_names.append(db_name)

db_name_dict = {}
for tip_name in tree_vertebrates.get_tip_names():
    db_name = find_db_prefix(tip_name.lower())
    if db_name is None:
        db_name_dict[tip_name] = tip_name
    else:
        db_name_dict[tip_name] = db_name

for key in db_name_dict.keys():
    if key.split('_')[0:2] == ['Mus', 'musculus']:
        if key != 'Mus_musculus_reference_CL57BL6_strain':
            db_name_dict[key] = 'mus_musculus_strain'
        else:
            db_name_dict[key] = 'mus_musculus'

for key in db_name_dict.keys():
    if key.split('_')[0:2] == ['Sus','scrofa']:
        if key != 'Sus_scrofa_reference_breed':
            db_name_dict[key] = 'sus_scrofa_strain'
        else:
            db_name_dict[key] = 'sus_scrofa'

tree_vertebrates.reassign_names(db_name_dict)
tree_vertebrates.reassign_names(common_names_mapping)


def process_path(path, output_dir):
    file_name = os.path.basename(path).rsplit('.', 1)[0] + ".json"
    print(file_name)
    alignment = deserialise_object(json.load(open(path, 'r')))
    species_names = alignment.names
    tip_names = species_names + outgroup_species
    tree_sub = tree_vertebrates.get_sub_tree(tip_names)
    tree_topology = tree_sub.get_sub_tree(species_names)
    model = get_app("model", "GN", tree = tree_topology, time_het="max", optimise_motif_probs=True)
    model_fitting_result = model(alignment)
    print('model_fitting_done')
    lf_result = model_fitting_result.lf
    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)
    write_json_app(lf_result, identifier=file_name)
    print(f"Result saved to {os.path.join(output_dir, file_name)}")

def main(paths, num_processes, output_dir, num_file):
    if num_file is not None:
        paths = paths[:num_file]
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir) for path in paths])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit model for sequences in different paths.')
    parser.add_argument('input', type=str, help='Path to the directory containing sequence files or a JSON file listing the paths')
    parser.add_argument('-n', '--num_processes', type=int, default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
    parser.add_argument('-o', '--output_dir', type=str, default='/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result', help='Output directory for JSON files (default: /Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result)')
    parser.add_argument('-l', '--limit', type=int, help='Limit the number of files to process')


    args = parser.parse_args()
    input_path = args.input
    num_processes = args.num_processes
    output_dir = args.output_dir
    num_file = args.limit
    
    if os.path.isdir(input_path):
        paths = glob.glob(os.path.join(input_path, "*.json"))

    main(paths, num_processes, output_dir, num_file)
