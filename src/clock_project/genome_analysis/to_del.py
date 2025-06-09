import click
from cogent3 import make_tree, get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from scitrack import CachingLogger
import uuid
from pathlib import Path
from cogent3.app import evo
import shutil
import os
from cogent3.util import parallel


def configure_parallel(parallel_option: bool, PBS_NCPUS: int) -> dict:
    """returns parallel configuration settings for use as composable.apply_to(**config)"""
    mpi = None if PBS_NCPUS-1 < 2 else PBS_NCPUS-1  # no point in MPI if < 2 processors
    parallel_option = True if mpi else parallel_option
    par_kw = dict(max_workers=mpi, use_mpi=True) if mpi else None

    return {"parallel": parallel_option, "par_kw": par_kw}

def get_id(result):
    return result.source.unique_id

@define_app
def test_hypothesis_clock_model(aln: AlignedSeqsType, tree=None, opt_args=None, num_reps=10) -> SerialisableType:
    outgroup_name = aln.info["triples_species_name"]["outgroup"]
    tree = make_tree(tip_names=aln.names)
    sp1 = aln.info["triples_species_name"]["ingroup1"]
    sp2 = aln.info["triples_species_name"]["ingroup2"]
    outgroup_edge = [outgroup_name]

    model_kwargs = dict(
        tree=tree,
        opt_args=opt_args,
        lf_args=dict(discrete_edges=[outgroup_edge]),
        optimise_motif_probs=True,
    )
    null = get_app(
        "model",
        "GN",
        name="clock",
        param_rules=[dict(par_name="length", edges=[sp1, sp2], is_independent=False)],
        **model_kwargs,
    )
    alt = get_app("model", "GN", name="no-clock", **model_kwargs)
    hyp = get_app("hypothesis", null, alt)
    bootstrapper = evo.bootstrap(hyp, num_reps=num_reps, parallel=True)
    result = bootstrapper(aln)
    return result

# @define_app
# def minimal_test(aln: AlignedSeqsType) -> SerialisableType:
#     """A minimal test app that returns the input alignment."""
#     # This is a placeholder for the actual analysis you want to perform
#     return aln


_click_command_opts = {
    "no_args_is_help": True,
    "context_settings": {"show_default": True},
}

@click.command(**_click_command_opts)
@click.argument("input_path", type=click.Path(exists=True))
@click.option("--output_dir", "-o", type=click.Path(), help="Output directory")
@click.option("--limit", "-l", type=int, help="limit for number of files")
@click.option(
    "--num_reps", "-r", type=int, default=100, help="Number of bootstrap replicates"
)
@click.option("--mpi", "-m", type=int, default=0, help="Number of MPI processes to use")
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

def main(input_path, output_dir, limit, num_reps, mpi, parallel_option, force):

    print(f"{parallel.is_master_process()}")

    # Convert to Path right away
    output_dir = Path(output_dir)

    if force and output_dir.exists():
        shutil.rmtree(output_dir, ignore_errors=True)

    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=True)

    outpath = output_dir / f"{uuid.uuid4().hex}.log"

    LOGGER = CachingLogger(log_file_path=outpath, create_dir=True)
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.log_versions("cogent3")
    LOGGER.log_versions("clock_project")
    # Configure parallel processing

    # Build minimal pipeline
    loader = get_app("load_json")
    

    writer = get_app("write_json", 
                   data_store=open_data_store(output_dir, mode="w", suffix="json"), id_from_source=get_id)
    
    # pipeline = loader + test_hypothesis_clock_model(num_reps = num_reps) + writer
    pipeline = loader + test_hypothesis_clock_model(num_reps = num_reps) + writer


    parallel_config = configure_parallel(
        parallel_option=parallel_option, 
        PBS_NCPUS = mpi
    )

    # Apply to data
    input_dstore = open_data_store(input_path, suffix="json")
    pipeline.apply_to(
        input_dstore[0:limit],
        show_progress=True,
        logger=LOGGER,
        **parallel_config
    )

    click.echo("Completed edited clock_bootstrape test")

if __name__ == "__main__":
    main()

