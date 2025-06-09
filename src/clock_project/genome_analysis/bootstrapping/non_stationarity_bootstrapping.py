from cogent3.evolve.models import register_model
from cogent3.evolve.ns_substitution_model import GeneralStationary
from cogent3 import make_tree, get_app, open_data_store, get_moltype
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from scitrack import CachingLogger
import uuid
import click
from cogent3.app import evo
from pathlib import Path
import shutil


def configure_parallel(parallel_option: bool, PBS_NCPUS: int) -> dict:
    """returns parallel configuration settings for use as composable.apply_to(**config)"""
    mpi = None if PBS_NCPUS-1 < 2 else PBS_NCPUS-1  # no point in MPI if < 2 processors
    parallel_option = True if mpi else parallel_option
    par_kw = dict(max_workers=mpi, use_mpi=True) if mpi else None

    return {"parallel": parallel_option, "par_kw": par_kw}

RATE_PARAM_UPPER = 50

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
    bootstrapper = evo.bootstrap(hyp, num_reps=num_reps, parallel=False)
    result = bootstrapper(aln)    
    return result


_click_command_opts = {
    "no_args_is_help": True,
    "context_settings": {"show_default": True},
}


@click.command(**_click_command_opts)
@click.argument("input_path", type=Path)
@click.option("--output_dir", "-o", type=Path)
@click.option("--limit", "-l", type=int, help="limit for number of files")
@click.option("--mpi", "-m", type=int, default=0, help="Number of MPI processes to use")
@click.option(
    "--num_reps", "-r", type=int, default=100, help="Number of bootstrap replicates"
)
@click.option(
    "-p",
    "--parallel_option",
    is_flag=True,
    default=False,
    help="run in parallel (on single machine)",
)

@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="Force overwrite output directory by deleting existing content.",
)

def main(input_path, output_dir, mpi, limit, num_reps, parallel_option, force):

    if force and output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    outpath = output_dir / f"{uuid.uuid4().hex}.log"
    outpath = Path(output_dir) / f"{uuid.uuid4().hex}.log"

    LOGGER = CachingLogger(log_file_path=outpath, create_dir=True)
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.log_versions("cogent3")
    LOGGER.log_versions("clock_project")

    loader = get_app("load_json")
    

    writer = get_app("write_json", 
                   data_store=open_data_store(output_dir, mode="w", suffix="json"), id_from_source=get_id)

    input_dstore = open_data_store(input_path, suffix="json")

    app = loader + test_hypothesis_non_stationary_model(num_reps =  num_reps) + writer

    parallel_config = configure_parallel(
        parallel=parallel_option, PBS_NCPUS=mpi
    )

    app.apply_to(
        input_dstore[0:limit],
        show_progress=True,
        logger=LOGGER,
        **parallel_config,
    )

    print("finished")


if __name__ == "__main__":
    main()

