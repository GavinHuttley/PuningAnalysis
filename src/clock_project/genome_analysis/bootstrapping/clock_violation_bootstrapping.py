from cogent3 import get_app, open_data_store
from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_moltype
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, HypothesisResultType, BootstrapResultType, IdentifierType
from cogent3.app import evo
import json


def get_param_rules_upper_limit(model_name, upper):
    """rules to set the upper value for rate matrix terms"""
    from cogent3 import get_model

    sm = get_model(model_name)
    return [{"par_name": par_name, "upper": upper} for par_name in sm.get_param_list()]

RATE_PARAM_UPPER = 50

@define_app
def test_hypothesis_model(aln_path: IdentifierType, result_lf_path, tree=None, opt_args=None) -> HypothesisResultType:
    aln = json.load(open(aln_path, 'r'))
    result_lf = load_json_app(result_lf_path)

    tree = result_lf.tree

    model_kwargs = dict(
    tree=tree,
    opt_args=opt_args,
    # unique_trees=True,
    optimise_motif_probs=True,
    )
    null = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GSN", RATE_PARAM_UPPER),
            **model_kwargs,
        )
    alt = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            time_het = 'max',
            **model_kwargs,
        )
    
    hyp = evo.hypothesis(null, alt, sequential=True)
    result = hyp(aln)

    return result.to_json()



load_json_app = get_app("load_json")

result_lf_dir = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result_350_threshold'

aln_dir = '/Users/gulugulu/repos/PuningAnalysis/data/ensembl_ortholog_sequences/homologies_alignment_common_name_350_threshold'

aln_paths = open_data_store(aln_dir, suffix='json')

result_lf_path = open_data_store(result_lf_dir, suffix='json')

output_path = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/clock_violation_test_result.json'

import click
import multiprocessing
import glob
import os
@click.command()
@click.argument('aln_dir', type=click.Path(exists=True))
@click.option('--result_lf_dir', '-rd', help='Path to the directory containing model fitting result')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_path', '-o', type=click.Path(), help='Output directory')

def main(aln_dir, result_lf_dir, num_processes, output_path):
    """Fit model for sequences in different paths and aggregate the results into a single JSON file."""
    result_dir = {}
    paths = glob.glob(os.path.join(aln_dir, '*.json'))  # corrected to match .json files

    with multiprocessing.Pool(processes=num_processes) as pool:
        for path in paths:
            gene_name = os.path.basename(path).rsplit('.', 1)[0]
            result_lf_path = os.path.join(result_lf_dir, gene_name + ".json")
            result = pool.apply_async(test_hypothesis_model, (path, result_lf_path))
            result_dir[gene_name] = result.get()  # Get the result from the async call

    # Write the aggregated results to a JSON file
    with open(output_path, 'w') as outfile:
        json.dump(result_dir, outfile, indent=4)
    
    print(f"Results saved to {output_path}")

if __name__ == '__main__':
    main()