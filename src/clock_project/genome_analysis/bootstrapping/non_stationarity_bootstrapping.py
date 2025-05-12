from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_app, open_data_store, get_moltype
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from cogent3.app import io as io_app
from scitrack import CachingLogger
import uuid
import click
from cogent3.app import evo
from pathlib import Path


def configure_parallel(parallel: bool, mpi: int, num_processes: int) -> dict:
    mpi = None if mpi < 2 else mpi
    parallel = True if mpi else parallel
    par_kw = (
        dict(max_workers=mpi, use_mpi=True)
        if mpi
        else dict(max_workers=num_processes, use_mpi=False)
    )
    return {"parallel": parallel, "par_kw": par_kw}

RATE_PARAM_UPPER = 50

loader_db = io_app.load_db()

def get_id(result):
    return result.source.unique_id

@register_model("nucleotide")
def GSN(**kwargs):
    """A General Stationary Nucleotide substitution model instance."""
    kwargs["optimise_motif_probs"] = kwargs.get("optimise_motif_probs", True)
    kwargs["name"] = kwargs.get("name", "GSN")
    return GeneralStationary(get_moltype("dna").alphabet, **kwargs)

def get_param_rules_upper_limit(model_name, upper):
    """rules to set the upper value for rate matrix terms"""
    from cogent3 import get_model

    sm = get_model(model_name)
    return [{"par_name": par_name, "upper": upper} for par_name in sm.get_param_list()]


@define_app
def test_hypothesis_non_stationary_model(aln: AlignedSeqsType, num_reps = 100, tree=None, opt_args=None) -> SerialisableType:
    outgroup_name = aln.info['triples_species_name']['outgroup']
    tree = make_tree(tip_names=aln.names)
    
    outgroup_edge = [outgroup_name]

    model_kwargs = dict(
    tree=tree,
    opt_args=opt_args,
    # unique_trees=True,
    lf_args=dict(discrete_edges=[outgroup_edge]),
    optimise_motif_probs=True,
    )
    null = evo.model(
            "GSN",
            param_rules=get_param_rules_upper_limit("GSN", RATE_PARAM_UPPER),
            **model_kwargs,
        )
    alt = evo.model(
            "GN",
            param_rules=get_param_rules_upper_limit("GN", RATE_PARAM_UPPER),
            **model_kwargs,
        )
    
    hyp = evo.hypothesis(null, alt, sequential=True)
    bootstrapper = evo.bootstrap(hyp, num_reps=num_reps, parallel=True)
    result = bootstrapper(aln)    
    return result


_click_command_opts = {
    "no_args_is_help": True,
    "context_settings": {"show_default": True},
}


@click.command(**_click_command_opts)
@click.argument("input_path", type=click.Path(exists=True))
@click.option(
    "--num_processes",
    "-n",
    type=int,
    help="Number of processes to use (default: number of CPUs)",
)
@click.option("--mpi", "-m", type=int, default=0, help="Number of MPI processes to use")
@click.option("--output_dir", "-o", type=click.Path(), help="Output directory")
@click.option("--limit", "-l", type=int, help="limit for number of files")
@click.option(
    "--num_reps", "-r", type=int, default=100, help="Number of bootstrap replicates"
)
def main(input_path, num_processes, mpi, output_dir, limit, num_reps):
    outpath = Path(output_dir) / f"{uuid.uuid4().hex}.log"

    LOGGER = CachingLogger(log_file_path=outpath, create_dir=True)
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.log_versions("cogent3")
    LOGGER.log_versions("clock_project")


    ton_bootstrapper = test_hypothesis_non_stationary_model(num_reps =  num_reps)

    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore, id_from_source=get_id)

    input_data_store = open_data_store(input_path)

    app = loader_db + ton_bootstrapper + write_json_app

    parallel_config = configure_parallel(
        parallel=True, num_processes=num_processes, mpi=mpi
    )

    app.apply_to(
        input_data_store[0:limit],
        show_progress=True,
        cleanup=True,
        logger=LOGGER,
        **parallel_config,
    )

    print("finished")


if __name__ == "__main__":
    main()

