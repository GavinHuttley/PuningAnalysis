
import os
import time
import uuid
from collections import Counter
from pathlib import Path
from cogent3.app import evo
import click
from cogent3 import make_tree, get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3.app.typing import AlignedSeqsType, SerialisableType
from cogent3.util import parallel
from scitrack import CachingLogger

def get_id(result):
    return result.source.unique_id

@define_app
def test_hypothesis_clock_model(
    aln: AlignedSeqsType, tree=None, opt_args=None, num_reps=10
) -> SerialisableType:
    """
    This is your existing Cogent3‐app function. It builds two models (clock vs no‐clock),
    runs the likelihood ratio hypothesis via evo.bootstrap, and returns the result.
    """
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


@click.command(no_args_is_help=True, context_settings={"show_default": True})
@click.argument("input_path", type=click.Path(exists=True, file_okay=False))
@click.argument("output_dir", type=click.Path())
@click.option(
    "--limit", "-l", type=int, default=None, help="Limit number of files to process."
)
@click.option(
    "--num_reps",
    "-r",
    type=int,
    default=100,
    help="Number of bootstrap replicates for each alignment.",
)
@click.option(
    "--force",
    is_flag=True,
    default=False,
    help="If set, delete and recreate output_dir before writing results.",
)
def main(input_path, output_dir, limit, num_reps, force):
    """
    Usage (once in your PBS script):
        mpirun -n $PBS_NCPUS python to_del_mpi.py \
          triples_representative_subset_json \
          mpi_results_dir \
          --limit 4 \
          --num_reps 2 \
          --force
    """
    output_dir = Path(output_dir)
    if force and output_dir.exists():
        for child in output_dir.iterdir():
            if child.is_file():
                child.unlink()
            else:
                import shutil

                shutil.rmtree(child, ignore_errors=True)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = output_dir / f"{uuid.uuid4().hex}.log"
    LOGGER = CachingLogger(log_file_path=log_file, create_dir=True)
    if parallel.get_rank() == 0:
        LOGGER.log_args()
        LOGGER.log_versions("numpy")
        LOGGER.log_versions("cogent3")
        LOGGER.log_versions("clock_project")

    # 3) Determine how many MPI ranks / “workers” we have
    #    PBS_NCPUS is the number of ranks you requested via “mpirun -n $PBS_NCPUS”
    PBS_NCPUS = int(os.environ.get("PBS_NCPUS", "1"))
    world_size = parallel.SIZE  # number of MPI ranks in COMM_WORLD
    my_rank = parallel.get_rank()

    if world_size != PBS_NCPUS:
        # It’s common for users to do “mpirun -n 8” and have PBS_NCPUS=8,
        # but just in case they mismatch, warn:
        print(
            f"[Rank {my_rank}] Warning: PBS_NCPUS={PBS_NCPUS} != parallel.SIZE={world_size}"
        )

    if my_rank == 0:
        print(f"[Rank 0] Running with {world_size} MPI ranks, PBS_NCPUS={PBS_NCPUS}")

    # 4) Build the Cogent3 apps we need:
    loader = get_app("load_json")
    writer = get_app("write_json", 
                   data_store=open_data_store(output_dir, mode="w", suffix="json"), id_from_source=get_id)
    


    all_items = list(open_data_store(input_path, suffix="json"))
    if limit is not None:
        all_items = all_items[:limit]

    def process_one(item):
        aln = loader(item)
        result = test_hypothesis_clock_model(aln, num_reps=num_reps)
        writer(result)

        return result.source.unique_id
    

    for _ in parallel.as_completed(
        process_one,
        all_items,
        use_mpi=True,
        max_workers=PBS_NCPUS,
    ):
        pass


if __name__ == "__main__":
    # The “if __name__” guard is required for mpi4py futures to work properly.
    main()
