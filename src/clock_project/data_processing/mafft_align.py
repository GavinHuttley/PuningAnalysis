
from cogent3 import get_app, open_data_store
from cogent3.app.composable import define_app
from cogent3 import get_app
from cogent3.app import typing
import pathlib
import contextlib
import tempfile
import subprocess
import sys
import click
import multiprocessing

@contextlib.contextmanager
def tempdir(working_dir: pathlib.Path | str | None = None) -> pathlib.Path:
    """context manager returns a temporary directory in working_dir"""
    with tempfile.TemporaryDirectory(dir=working_dir) as temp_dir:
        yield pathlib.Path(temp_dir)

def exec_command(
    cmnd: str,
    stdout: int = subprocess.PIPE,
    stderr: int = subprocess.PIPE,
) -> str | None:
    """Executes shell command and returns stdout if completed with exit code 0."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        sys.stderr.writelines(f"FAILED: {cmnd}\n{msg}")
        sys.exit(proc.returncode)
    return out.decode("utf8") if out is not None else None

@define_app
def get_mafft_aligned_seq(seqs: typing.SeqsCollectionType, gc=1) -> typing.AlignedSeqsType:
    """
    Loads sequences from the input directory, translates them to amino acids,
    aligns using MAFFT, and returns the aligned DNA sequence collection.

    Parameters
    ----------
    seqs_dir: str

    Returns
    -------
    str
        Path to the aligned amino acid FASTA file.
    """
    # Temporary directory context
    try:
        if seqs.num_seqs < 3:
            print('species number less than 3')
            return None
        else: 
            with tempdir() as temp_dir:
                aa_fasta_path = temp_dir / "aa_sequences.fasta"
                aligned_aa_path = temp_dir / "aligned_aa.fasta"
                translater = get_app("translate_seqs", gc=gc)

                # Load and translate the first FASTA file
                aa_seqs = translater(seqs)

                # Write translated amino acid sequences to temporary FASTA file
                aa_seqs.write(aa_fasta_path, format="fasta")

                # Build the MAFFT command
                mafft_command = f"mafft --amino {aa_fasta_path} > {aligned_aa_path}"
                print(f"Running MAFFT: {mafft_command}")

                # Execute the MAFFT command
                exec_command(mafft_command)

                # Load the aligned amino acid sequences
                loader_aligned = get_app("load_aligned", format="fasta")
                aligned_seq_collection = loader_aligned(str(aligned_aa_path)).to_type(array_align=True)        

                aligned_seqs = aligned_seq_collection.replace_seqs(seqs)
                aligned = True
                print('successfully aligned')

            return aligned_seqs
        
    except Exception as e:
        print("Cannot be aligned")
        return None

aligner = get_mafft_aligned_seq()   
load_json_app = get_app("load_json")


def process_path(path, write_json_app):
    file_name = path.unique_id
    print(file_name)
    seqs = load_json_app(path)
    aligned_seqs = aligner(seqs)
    if aligned_seqs is not None:
        write_json_app(aligned_seqs, identifier=file_name)
    else:
        print(f"Skipping file {file_name} due to alignment failure or insufficient sequences.")

@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('-o', '--output_dir', help='Output directory for JSON files')
@click.option('-n', '--num_processes', default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process')
def main(input, num_processes, output_dir, limit):
    """
    Fit model for sequences in different paths.
    
    INPUT: Path to the directory containing sequence files or a JSON file listing the paths.
    """
    input_data_store = open_data_store(input, suffix='json', limit=limit)
    out_dstore = open_data_store(output_dir, mode="w", suffix="json")
    write_json_app = get_app("write_json", data_store=out_dstore)

    
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, write_json_app) for path in input_data_store.completed])

if __name__ == "__main__":
    main()



# @define_app
# def replace_common_species_names(alignment: SeqsCollectionType	
# ) -> SeqsCollectionType:
#     with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/genome_information/available_species_name.json', 'r') as infile1:
#         available_species_names = json.load(infile1)

#     with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/genome_information/common_names.json', 'r') as infile2:
#         common_names = json.load(infile2)
#     name_map = dict(zip(available_species_names, common_names))
    
#     new_names = {}
    
#     # Iterate through the sequence names in the alignment
#     for name in alignment.names:
#         species_info, gene_info = name.split('-', 1)
        
#         # Check if the species part exists in the name mapping
#         if species_info in name_map:
#             # Create a new name using the common name and the gene info
#             new_name = name_map[species_info]
#             new_names[name] = new_name
#         else:
#             # If no mapping exists, keep the original name
#             new_names[name] = name
    
#     # Update the sequence names in the alignment with new names
#     return alignment.rename_seqs(lambda x: new_names[x])

# common_name_renamer = replace_common_species_names()
