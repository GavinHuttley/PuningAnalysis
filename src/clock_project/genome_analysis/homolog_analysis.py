import os
import glob
import re
from cogent3.maths.measure import jsd
import numpy as np
from cogent3 import get_app, load_aligned_seqs
import plotly.graph_objects as go
from plotly.subplots import make_subplots

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
            'P_value':value['p_value']
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
                x=np.log(df['JSD Value']),
                y=df['Distance Value'],
                mode='markers'
            ),
            row=row,
            col=col
        )
        # Only add x-axis title to bottom row plots
        if row == rows:
            fig.update_xaxes(title_text="Log(JSD Value)", row=row, col=col)
        
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
