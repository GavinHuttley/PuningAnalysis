import argparse
import glob
import json
import os
from cogent3 import get_app, open_data_store
from clock_project.genome_analysis.homolog_analysis import common_name_renamer, filter

# App initializations
loader = get_app("load_unaligned", format="fasta", moltype="dna")
codon_aligner = get_app("progressive_align", "codon", unique_guides=True)
cpos3 = get_app("take_codon_positions", 3)

def get_just3rd_aligned_no_degenerates(path):
    """Processes a sequence file to align and filter third codon positions without degenerates."""
    seqs = loader(path)
    seqs_filtered = filter(seqs)
    seqs_filtered_no_stop_codon = seqs_filtered.trim_stop_codons(strict=True)
    aligned = codon_aligner(seqs_filtered_no_stop_codon)
    just3rd_aligned = cpos3(aligned)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates(motif_length=3)
    return just3rd_aligned_no_degenerates, seqs_filtered

def process_alignment_and_save(path, output_dir1, output_dir2):
    """Processes the file, renames, gets the alignment, and saves it to JSON."""
    try:
        # Get alignment
        alignment, seqs_filtered = get_just3rd_aligned_no_degenerates(path)
        
        # Rename sequences in the alignment
        renamed_alignment = common_name_renamer(alignment)
        renamed_seqs = common_name_renamer(seqs_filtered)
        
        # Save alignment to JSON
        file_name = os.path.basename(path).rsplit('.', 1)[0] + ".json"

        output_path_alignment = os.path.join(output_dir1, file_name)
        output_path_sequence = os.path.join(output_dir2, file_name)

        with open (output_path_alignment, 'w') as outfile1:
            json.dump(renamed_alignment.to_json(), outfile1)
        print(f"Alignment saved to {output_path_alignment}")

        with open (output_path_sequence, 'w') as outfile2:
            json.dump(renamed_seqs.to_json(), outfile2)
        print(f"Alignment saved to {output_path_sequence}")

    except Exception as e:
        print(f"Failed to process {path}: {e}")

def main(input_dir, output_dir1, output_dir2):
    """Processes each file in the directory, getting alignments, renaming, and saving them as JSON."""
    paths = glob.glob(os.path.join(input_dir, "*.fa"))
    for path in paths:
        process_alignment_and_save(path, output_dir1, output_dir2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process sequence files, rename, and save alignments in JSON format.')
    parser.add_argument('input_dir', type=str, help='Directory containing sequence files.')
    parser.add_argument('-o1', '--output_dir1', type=str, default='./alignment_results', help='Output directory for JSON files (default: current directory)')
    parser.add_argument('-o2', '--output_dir2', type=str, default='./sequence_results', help='Output directory for JSON files (default: current directory)')


    args = parser.parse_args()
    main(args.input_dir, args.output_dir1, args.output_dir2)
