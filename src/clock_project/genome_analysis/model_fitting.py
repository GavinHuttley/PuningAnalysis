import argparse
import multiprocessing
from cogent3 import get_app
from clock_project.genome_analysis.homolog_analysis import filter, common_name_renamer
import json
import os
import glob
from cogent3 import get_app, open_data_store

loader = get_app("load_unaligned", format="fasta", moltype="dna")
codon_aligner = get_app("progressive_align", "codon", unique_guides=True)
cpos3 = get_app("take_codon_positions", 3)

def get_just3rd_aligned_no_degenerates(path):
    seqs = loader(path)
    seqs_filtered = filter(seqs)
    seqs_filtered_no_stop_codon = seqs_filtered.trim_stop_codons(strict=True)
    aligned = codon_aligner(seqs_filtered_no_stop_codon)
    just3rd_aligned = cpos3(aligned)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates(motif_length=3)
    return just3rd_aligned_no_degenerates

dist_cal = get_app("fast_slow_dist", fast_calc="tn93", moltype="dna")
est_tree = get_app("quick_tree", drop_invalid=True)
tree_func = dist_cal + est_tree
model = get_app("model", "GN", tree_func=tree_func, time_het="max")

def fit_model_for_path(path):
    just3rd_aligned_no_degenerates = get_just3rd_aligned_no_degenerates(path)
    alignment_common_name = common_name_renamer(just3rd_aligned_no_degenerates)
    res = model(alignment_common_name)
    return res.lf

def process_path(path, output_dir):
    try:
        lf_result = fit_model_for_path(path)
        file_name = os.path.basename(path).rsplit('.', 1)[0] + ".json"
        out_dstore = open_data_store(output_dir, mode="w", suffix="json")
        write_json_app = get_app("write_json", data_store=out_dstore)
        write_json_app(lf_result, identifier=file_name)
        print(f"Result saved to {os.path.join(output_dir, file_name)}")
    except Exception as e:
        print(f"Failed to process {path}: {e}")

def main(paths, num_processes, output_dir):
    with multiprocessing.Pool(processes=num_processes) as pool:
        pool.starmap(process_path, [(path, output_dir) for path in paths])

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Fit model for sequences in different paths.')
    parser.add_argument('input', type=str, help='Path to the directory containing sequence files or a JSON file listing the paths')
    parser.add_argument('-n', '--num_processes', type=int, default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
    parser.add_argument('-o', '--output_dir', type=str, default='/Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result', help='Output directory for JSON files (default: /Users/gulugulu/repos/PuningAnalysis/results/output_data/model_fitting_result)')

    args = parser.parse_args()
    input_path = args.input
    num_processes = args.num_processes
    output_dir = args.output_dir
    
    if os.path.isdir(input_path):
        paths = glob.glob(os.path.join(input_path, "*.fa"))
    else:
        with open(input_path, 'r') as f:
            paths = json.load(f)

    main(paths, num_processes, output_dir)
