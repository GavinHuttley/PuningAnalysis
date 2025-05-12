from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_moltype, get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from cogent3.app import evo
import click
RATE_PARAM_UPPER = 50

def get_id(result):
    return result.source.unique_id


def get_param_rules_upper_limit(model_name, upper):
    """rules to set the upper value for rate matrix terms"""
    from cogent3 import get_model

    sm = get_model(model_name)
    return [{"par_name": par_name, "upper": upper} for par_name in sm.get_param_list()]


@define_app
def test_hypothesis_clock_model_N(aln: AlignedSeqsType, tree=None, opt_args=None) -> SerialisableType:
    tree =  make_tree(tip_names=aln.names)
    sp1 = aln.info['triples_species_name']['ingroup1']
    sp2 = aln.info['triples_species_name']['ingroup2']

    outgroup_name = aln.info['triples_species_name']['outgroup']
    print(outgroup_name)
    outgroup_edge = [outgroup_name]

    model_kwargs = dict(
    tree=tree,
    opt_args=None,
    # unique_trees=True,
    lf_args=dict(discrete_edges=[outgroup_edge]),
    optimise_motif_probs=True,
    )

    null = get_app("model", "GN", name="clock", param_rules=[dict(par_name="length", edges=[sp1, sp2], is_independent=False)], **model_kwargs)
    alt = get_app("model", "GN", name="no-clock", **model_kwargs)
    hyp = get_app("hypothesis", null, alt)

    # clock_p_value_observed[gene_name] = (result.pvalue)
    bootstrapper = evo.bootstrap(hyp, num_reps=100, parallel=False)
    result_test = bootstrapper(aln)    
    print('finish')
    return result_test

load_json_app = get_app("load_json")

def p_value(result):
    return sum(result.observed.LR <= null_lr for null_lr in result.null_dist) / len(result.null_dist)

clock_bootstrapper2 = test_hypothesis_clock_model_N()


@click.command()
@click.argument('input_path', type=click.Path(exists=True))
@click.option('--num_processes', '-n', type=int, help='Number of processes to use (default: number of CPUs)')
@click.option('--output_dir', '-o', type=click.Path(), help='Output directory')
def main(input_path, num_processes, output_dir):

    print(input_path)

    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore, id_from_source = get_id)

    input_data_store = open_data_store(input_path, suffix= 'json')

    clock_bootstrap_process = load_json_app + clock_bootstrapper2 + write_json_app


    clock_bootstrap_process.apply_to(
        input_data_store,
        parallel=True,
        par_kw=dict(
            max_workers=num_processes, use_mpi=False
        ),
    )
    print("Finish app call")
    clock_bootstrap_process.disconnect()


if __name__ == "__main__":
    main()