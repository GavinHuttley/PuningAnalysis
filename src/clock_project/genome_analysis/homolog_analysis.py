import os
import glob
import re
from cogent3.maths.measure import jsd
import numpy as np
from cogent3 import get_app, load_aligned_seqs
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json

import pandas as pd

base_dir1 = "/Users/gulugulu/sampled_homology_200"
def gather_fasta_paths(base_dir):
    pattern = os.path.join(base_dir, '*.fasta')
    # Use glob.glob to find all files matching the pattern
    fasta_files = glob.glob(pattern)
    return fasta_files

def gather_fasta_paths2(base_dir):
    pattern = os.path.join(base_dir, '*.fa')
    # Use glob.glob to find all files matching the pattern
    fasta_files = glob.glob(pattern)
    return fasta_files

def extract_info(path):
    match = re.search(r'/([^/]+)\.fasta$', path)
    if match:
        return match.group(1)
    else:
        return "unknown"
    
def pairwise_jsd_matrix(species_data):
    species_keys = list(species_data.keys())
    num_species = len(species_keys)
    jsd_matrix = np.zeros((num_species, num_species))  # Initialize a square matrix

    for i, species_1 in enumerate(species_keys):
        for j, species_2 in enumerate(species_keys):
            if i < j:  # To avoid recomputation, calculate only for i < j
                jsd_value = jsd(species_data[species_1], species_data[species_2])
                jsd_matrix[i, j] = jsd_value
                jsd_matrix[j, i] = jsd_value  # JSD is symmetric

    return jsd_matrix

loader = get_app("load_unaligned", format="fasta", moltype="dna")
codon_aligner = get_app("progressive_align", "codon", unique_guides = True)
cpos3 = get_app("take_codon_positions", 3)
# omit_degens = get_app("omit_degenerates", moltype="dna")
def length_divisible_by_three(row):
    return len(row) % 3 == 0


def prepare_dataframes(info_dict):
    dataframes_dict = {}
    for key, value in info_dict.items():
        data = {
            'JSD Value': value['pairwise_jsd'],
            'Distance Value': value['pairwise_distance'],
            'Correlation Factor': value['col'],
            'P_value':value['p_value'],
        }
        dataframes_dict[key] = pd.DataFrame(data)
    return dataframes_dict

def jsd_genetic_distance_scatters(dataframes_dict):
    keys = list(dataframes_dict.keys())
    rows = int(len(keys) ** 0.5) + 1  # Calculate the number of rows for subplots
    cols = (len(keys) + rows - 1) // rows  # Calculate the number of columns

    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[f'{key}' for key in keys])
    
    # Populate subplots
    for index, (key, df) in enumerate(dataframes_dict.items(), start=1):
        row = (index - 1) // cols + 1
        col = (index - 1) % cols + 1
        
        fig.add_trace(
            go.Scatter(
                x=np.log10(df['JSD Value']),
                y=df['Distance Value'],
                mode='markers'
            ),
            row=row,
            col=col
        )
        # Only add x-axis title to bottom row plots

        fig.update_xaxes(row=row, col=col, range =[-6, 0])
        fig.update_yaxes(row=row, col=col, range =[0, 1.1])

        if row == rows:
            fig.update_xaxes(title_text="Log(JSD Value)", row=row, col=col, range=[-6, 0])
        
        # Only add y-axis title to first column plots
        if col == 1:
            fig.update_yaxes(title_text="Genetic Distance", row=row, col=col)
    
    fig.update_layout(
        height=300 * rows,  # Set a reasonable height based on number of rows
        width=300 * cols,   # Set a reasonable width based on number of columns
        title_text="Scatter Plots of 3rd codon position JSD Vs Genetic Distance",
        showlegend=False
    )
    
    return fig

import numpy as np
import pandas as pd
def get_score_at_quantile(data, target_quantile):
    """
    Returns the score value at the specified quantile based on the scores in the data.

    :param data: Dictionary with gene names as keys and scores as values
    :param target_quantile: The quantile (between 0 and 1) for which to find the score value
    :return: The score value at the specified quantile
    """
    if not 0 <= target_quantile <= 1:
        raise ValueError("Target quantile must be between 0 and 1")

    # Convert the dictionary to a pandas Series
    scores = pd.Series(data.values())

    # Calculate the score value at the target quantile
    score_at_quantile = scores.quantile(target_quantile)

    return score_at_quantile

from c3wrangling.score import seed_and_extend_smith_waterman

def get_matching_score_to_human(seqs):

    reference_name = seqs.names[-1]

    reference = seqs.get_seq(reference_name)
    seq_scores = {}
    for seqid in seqs.names:
        if seqid != reference.name:
            seq_scores[seqid] =  seed_and_extend_smith_waterman(
                    reference,
                    seqs.get_seq(seqid),
                    k=7,
                )
    
    return seq_scores


from cogent3.app.composable import define_app
from cogent3.app.typing import UnalignedSeqsType, SeqsCollectionType	
import numpy as np

loader = get_app("load_unaligned", format="fasta", moltype="dna")
codon_aligner = get_app("progressive_align", "codon", unique_guides=True)
cpos3 = get_app("take_codon_positions", 3)

@define_app
def too_many_ambigs(seqs: UnalignedSeqsType, frac=0.05) -> UnalignedSeqsType:
    w_ambig = seqs.get_lengths(include_ambiguity=True)
    wo_ambig = seqs.get_lengths(include_ambiguity=False)
    keep = {name for name in wo_ambig.keys() if wo_ambig[name] / w_ambig[name] >= 1 - frac}
    return seqs.take_seqs(keep)

@define_app
def low_matching_significance(seqs: UnalignedSeqsType, quantile=0.1) -> UnalignedSeqsType:
    ref_name = seqs.names[-1]
    ref = seqs.get_seq(seqs.names[-1])

    seq_scores = get_matching_score_to_human(seqs)

    quantile_score = get_score_at_quantile(seq_scores, quantile)
    seqs = seqs.take_seqs(ref_name, negate=True)
    seq_match_filtered = seqs.take_seqs_if(
        lambda seq: seq_scores[seq.name] > quantile_score
    )

    seq_match_filtered.add_seqs([ref])

    return seq_match_filtered

@define_app
def length_divisible_by_three(seqs: UnalignedSeqsType) -> UnalignedSeqsType:
    valid_seqs = seqs.take_seqs_if(
        lambda seq: len(seq) % 3 == 0
    )
    return valid_seqs

@define_app
def short_seq(seqs: UnalignedSeqsType, length: int = 600) -> UnalignedSeqsType:
    valid_seqs = seqs.take_seqs_if(
        lambda seq: len(seq) >= length
    )
    return valid_seqs

@define_app
def remove_redundent_seq(seqs: UnalignedSeqsType) -> UnalignedSeqsType:
    seq_names = seqs.names
    genus_list = []
    seq_names_keep = []
    seq_names_remove = []
    for name in seq_names:
        genus = name.split('_')[0]
        if genus not in genus_list:
            genus_list.append(genus)
            seq_names_keep.append(name)
        else:
            seq_names_remove.append(name)

    valid_seqs = seqs.take_seqs(seq_names_keep)
    return valid_seqs


drop_low_matching = low_matching_significance(quantile=0.1)
drop_ambiguous = too_many_ambigs(frac=0.2)
drop_invalid_length = length_divisible_by_three()
drop_short_seq = short_seq(length=600)
seq_without_redundent = remove_redundent_seq()
filter = drop_short_seq + drop_invalid_length + drop_low_matching + drop_ambiguous + seq_without_redundent

def remove_redundent_align(path):
    seqs = loader(path)
    seq_names = seqs.names
    genus_list = []
    seq_names_keep = []
    seq_names_remove = []
    for name in seq_names:
        genus = name.split('_')[0]
        if genus not in genus_list:
            genus_list.append(genus)
            seq_names_keep.append(name)
        else:
            seq_names_remove.append(name)

    new_seqs = seqs.take_seqs(seq_names_keep)
    new_seqs_filtered = filter(new_seqs)
    seqs_filtered_no_stop_codon = new_seqs_filtered.trim_stop_codons(strict=True)
    aligned = codon_aligner(seqs_filtered_no_stop_codon)
    just3rd_aligned = cpos3(aligned)
    just3rd_aligned_no_degenerates = just3rd_aligned.no_degenerates(motif_length=3)

    return just3rd_aligned_no_degenerates


@define_app
def replace_common_species_names(alignment: SeqsCollectionType	
) -> SeqsCollectionType:
    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/available_species_name.json', 'r') as infile1:
        available_species_names = json.load(infile1)

    with open ('/Users/gulugulu/repos/PuningAnalysis/results/output_data/common_names.json', 'r') as infile2:
        common_names = json.load(infile2)
    name_map = dict(zip(available_species_names, common_names))
    
    new_names = {}
    
    # Iterate through the sequence names in the alignment
    for name in alignment.names:
        species_info, gene_info = name.split('-', 1)
        
        # Check if the species part exists in the name mapping
        if species_info in name_map:
            # Create a new name using the common name and the gene info
            new_name = name_map[species_info]
            new_names[name] = new_name
        else:
            # If no mapping exists, keep the original name
            new_names[name] = name
    
    # Update the sequence names in the alignment with new names
    return alignment.rename_seqs(lambda x: new_names[x])

common_name_renamer = replace_common_species_names()

