
import json
from cogent3 import get_app, open_data_store
from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_moltype
from cogent3.app.composable import define_app, NotCompleted
from cogent3.app.typing import AlignedSeqsType, HypothesisResultType, SerialisableType, IdentifierType
from cogent3.util.deserialise import deserialise_object
import multiprocessing
import click


from cogent3.app import evo

def get_param_rules_upper_limit(model_name, upper):
    """rules to set the upper value for rate matrix terms"""
    from cogent3 import get_model

    sm = get_model(model_name)
    return [{"par_name": par_name, "upper": upper} for par_name in sm.get_param_list()]

load_json_app = get_app("load_json")

RATE_PARAM_UPPER = 50


def test_hypothesis_clock_model_whole_gene(aln_path: IdentifierType, result_lf_dir, output_dir, opt_args=None) -> HypothesisResultType:
    aln = deserialise_object(json.load(open(aln_path, 'r')))
    print('start')
    file_name = os.path.basename(aln_path.rstrip('/'))
    print(file_name)
    result_lf_path = os.path.join(result_lf_dir, file_name)

    tree = load_json_app(result_lf_path).get_ens_tree()

    model_kwargs = dict(
    tree=tree,
    opt_args=opt_args,
    # unique_trees=True,
    optimise_motif_probs=True,
    )
    null = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            time_het = None,
            **model_kwargs,
        )
    alt = evo.model(
            "GN",
            name="GN-max-het",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            time_het = "max",
            **model_kwargs,
        )
    print('begin hypothesis test')
    hyp = evo.hypothesis(null, alt, sequential=True)
    result = hyp(aln)
    print('end hypothesis test')
    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)
    write_json_app(result, identifier=file_name)

import glob
import os
result_lf_dir = '/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result_350_threshold'
gene_aln_dir = '/Users/gulugulu/repos/PuningAnalysis/data/ensembl_ortholog_sequences/homologies_alignment_common_name_350_threshold'


@click.command()
@click.argument('alignment_dir', type=click.Path(exists=True))
@click.option('--result_lf_dir', '-rd', help='Path to the directory containing model fitting result')
@click.option('--num_processes', '-n', default=multiprocessing.cpu_count(), type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory')
def main(alignment_dir, result_lf_dir, num_processes, output_dir):
    """Fit model for sequences in different paths."""
    gene_aln_paths = glob.glob(os.path.join(alignment_dir, '*.json'))

    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(test_hypothesis_clock_model_whole_gene, [(path, result_lf_dir, output_dir) for path in gene_aln_paths])

if __name__ == "__main__":
    main()