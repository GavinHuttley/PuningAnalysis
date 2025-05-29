import click
from cogent3 import get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from scitrack import CachingLogger
import uuid
from pathlib import Path

load_json_app = get_app("load_json")

def configure_parallel(parallel: bool, mpi: int) -> dict:
    """returns parallel configuration settings for use as composable.apply_to(**config)"""
    mpi = None if mpi < 2 else mpi  # no point in MPI if < 2 processors
    parallel = True if mpi else parallel
    par_kw = dict(max_workers=mpi, use_mpi=True) if mpi else None

    return {"parallel": parallel, "par_kw": par_kw}

@define_app
def dummy_processor(data: AlignedSeqsType)  -> SerialisableType:  # add dummy parameter here
    return data

_click_command_opts = {
    "no_args_is_help": True,
    "context_settings": {"show_default": True},
}

@click.command(**_click_command_opts)
@click.argument("input_path", type=click.Path(exists=True))
@click.option("--mpi", "-m", type=int, default=0, help="Number of MPI processes to use")
@click.option("--output_dir", "-o", type=click.Path(), help="Output directory")
@click.option("--limit", "-l", type=int, help="limit for number of files")

@click.option(
    "-p",
    "--parallel",
    is_flag=True,
    default=False,
    help="run in parallel (on single machine)",
)

def main(input_path, output_dir, mpi, limit, parallel):

    output_dir = Path(output_dir)
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
                   data_store=open_data_store(output_dir, mode="w", suffix="json"))
    
    pipeline = loader + dummy_processor() + writer

    parallel_config = configure_parallel(parallel=parallel, mpi=mpi)

    # Apply to data
    input_dstore = open_data_store(input_path, suffix="json")
    pipeline.apply_to(
        input_dstore[0:limit],
        show_progress=True,
        **parallel_config
    )

    click.echo("Completed minimal MPI test")

if __name__ == "__main__":
    main()