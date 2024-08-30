import json
import os
import glob
import numpy as np
from yapeng_check_BV import get_bounds_violation, load_param_values
from cogent3 import get_app, open_data_store
from cogent3.app.composable import NotCompleted

# Define global variables
load_json_app = get_app("load_json")
base_dir = '/Users/gulugulu/Desktop/honours/data_local_2/triples_model_fitting_550_threshold'

def get_gene_paths(base_dir):
    """Get a list of paths to gene directories under the base directory."""
    return glob.glob(os.path.join(base_dir, '*/'))

def process_model_fitting_results(gene_path, bounary_violation_function):
    """Process model fitting results for a given gene path."""
    model_fitting_result_dir = os.path.join(gene_path, 'model_fitting_result')
    model_fitting_results_input_dstore = open_data_store(model_fitting_result_dir, suffix='json')
    
    parameter_proximities = {'proximity_lower': [], 'proximity_upper': [], 'ens': []}
    valid_triads_identifier = []

    for path in model_fitting_results_input_dstore:
        identifier = path.unique_id.split('.')[0]
        model_fitting_result = load_json_app(path)
        if not isinstance(model_fitting_result, NotCompleted):
            param = load_param_values(model_fitting_result)
            exclude_params = ("length", "mprobs")
            list_of_params = param.params
            ens_list = model_fitting_result.lf.get_lengths_as_ens()

            for param in list_of_params:
                if param["par_name"] not in exclude_params:
                    proximity_lower = abs(param["init"] - param["lower"])
                    proximity_upper = abs(param["init"] - param["upper"])
                    ens = ens_list[param['edge']]
                    parameter_proximities['ens'].append(ens)
                    parameter_proximities['proximity_lower'].append(proximity_lower) 
                    parameter_proximities['proximity_upper'].append(proximity_upper)
            
            bounary_violation_check = bounary_violation_function.main(model_fitting_result)
            if bounary_violation_check.vio == []:
                valid_triads_identifier.append(identifier)
    
    return valid_triads_identifier

def save_valid_identifiers(valid_triads_identifier_dict):
    """Save the valid triads identifiers to a JSON file."""
    with open('/Users/gulugulu/Desktop/honours/data_local_2/valid_triads_identifier_new.json', 'w') as outfile:
        json.dump(valid_triads_identifier_dict, outfile, indent=4)

def process_gene_info(base_dir, valid_triads_identifier_dict):
    """Process gene information and store matrices and ens data."""
    all_ens_list = []
    matrices_dict_full = {}

    for gene_name, valid_triads_identifier in valid_triads_identifier_dict.items():
        triples_info_dir = os.path.join(base_dir, gene_name, 'triples_info_dict.json')
        triples_info_dict = json.load(open(triples_info_dir, 'r'))

        ens_pairs = []
        matrices_full_list = []

        for identifier in valid_triads_identifier:
            triples_species_name = triples_info_dict[identifier]['triples_species_names']
            ens_dict = triples_info_dict[identifier]['triples_info_small_tree']['ens']
            matrices = triples_info_dict[identifier]['triples_info_small_tree']['matrices']

            matrices_full = {
                'ingroup1': np.array(matrices[triples_species_name['ingroup1']]) * ens_dict[triples_species_name['ingroup1']],
                'ingroup2': np.array(matrices[triples_species_name['ingroup2']]) * ens_dict[triples_species_name['ingroup2']],
                'outgroup': np.array(matrices[triples_species_name['outgroup']]) * ens_dict[triples_species_name['outgroup']]
            }
            matrices_full_list.append(matrices_full)
            
            ens_pairs.append({
                triples_species_name['ingroup1']: ens_dict[triples_species_name['ingroup1']],
                triples_species_name['ingroup2']: ens_dict[triples_species_name['ingroup2']]
            })

        all_ens_list.extend(ens_pairs)
        matrices_dict_full[gene_name] = matrices_full_list

    return matrices_dict_full

def convert_arrays_to_lists(obj):
    """Recursively convert numpy arrays to lists for JSON serialization."""
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    elif isinstance(obj, dict):
        return {key: convert_arrays_to_lists(value) for key, value in obj.items()}
    elif isinstance(obj, list):
        return [convert_arrays_to_lists(element) for element in obj]
    else:
        return obj

def save_matrices_dict_full(matrices_dict_full):
    """Save the full matrices dictionary to a JSON file."""
    matrices_dict_full_converted = convert_arrays_to_lists(matrices_dict_full)
    with open('/Users/gulugulu/Desktop/honours/data_local_2/valid_matrix_full.json', 'w') as outfile:
        json.dump(matrices_dict_full_converted, outfile, indent=4)

def main():
    """Main function to run the entire script."""
    gene_paths = get_gene_paths(base_dir)
    bounary_violation_function = get_bounds_violation()

    valid_triads_identifier_dict = {}
    for path in gene_paths:
        file_name = os.path.basename(path.rstrip('/'))
        valid_triads_identifier = process_model_fitting_results(path, bounary_violation_function)
        valid_triads_identifier_dict[file_name] = valid_triads_identifier

    save_valid_identifiers(valid_triads_identifier_dict)
    matrices_dict_full = process_gene_info(base_dir, valid_triads_identifier_dict)
    save_matrices_dict_full(matrices_dict_full)

if __name__ == "__main__":
    main()
