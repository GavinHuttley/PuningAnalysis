import os
import click
from clock_project.genome_analysis.bootstrapping.non_stationarity_bootstrapping import bootstrapper
from cogent3 import get_app, open_data_store

from scitrack import CachingLogger

load_json_app = get_app("load_json")

# The following is a script to be run on Gadi
LOGGER = CachingLogger()

# the following environment variable is created by PBS on job execution
PBS_NCPUS = os.environ.get("PBS_NCPUS", None)
if PBS_NCPUS is None:
    PBS_NCPUS = 1  # Set a default for local testing
    # raise RuntimeError("did not get cpu number from environment")

PBS_NCPUS = int(PBS_NCPUS)


def get_id(result):
    return result.source.unique_id

@click.command()
@click.option("-l", "--log", type=str, default=None, help="name of log script")
def main(log):
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.log_versions("cogent3")
    LOGGER.log_versions("clock_project")
    LOGGER.log_file_path = f"/Users/gulugulu/Desktop/honours/data_local/logs/{log}.log"

    
    aln_dir_new = '/Users/gulugulu/Desktop/honours/data_local/triples_aln_subset_info_added'

    path_to_dir = '/Users/gulugulu/Desktop/honours/data_local/bootstrapping_test_non'
    out_dstore = open_data_store(path_to_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore, id_from_source = get_id)

    input_data_store = open_data_store(aln_dir_new, suffix= 'json', limit=10)

    bootstrap_process = load_json_app + bootstrapper + write_json_app

    bootstrap_process.apply_to(
        input_data_store,
        parallel=True,
        par_kw=dict(
            max_workers=PBS_NCPUS, use_mpi=True if PBS_NCPUS > 1 else False
        ),
    )

    bootstrap_process.disconnect()

if __name__ == "__main__":
    main()
