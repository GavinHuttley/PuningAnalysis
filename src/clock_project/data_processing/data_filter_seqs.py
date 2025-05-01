from cogent3 import get_app, open_data_store
import pandas as pd
from c3wrangling.score import seed_and_extend_smith_waterman
from cogent3.app.composable import define_app
from cogent3.app.typing import UnalignedSeqsType
import click
import multiprocessing

trim_stop = get_app("trim_stop_codons", gc=1)
cpos3 = get_app("take_codon_positions", 3)
translater = get_app("translate_seqs")
valid_cds = get_app("select_translatable", frame=1)


#functions used in the sequences filtering process

def length_divisible_by_three(row):
    return len(row) % 3 == 0

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


def get_matching_score_to_human(seqs):
    names = seqs.names
    for name in names:
        if name.split('-')[0] == 'homo_sapiens':
            reference_name = name

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

@define_app
def too_many_ambigs(seqs: UnalignedSeqsType, frac=0.02) -> UnalignedSeqsType:
    w_ambig = seqs.get_lengths(include_ambiguity=True)
    wo_ambig = seqs.get_lengths(include_ambiguity=False)
    keep = [name for name in wo_ambig.keys() if wo_ambig[name] / w_ambig[name] >= 1 - frac]
    return seqs.take_seqs(keep)

@define_app
def low_matching_significance(seqs: UnalignedSeqsType, quantile=0.05) -> UnalignedSeqsType:
    names = seqs.names
    for name in names:
        if name.split('-')[0] == 'homo_sapiens':
            ref_name = name
    ref = seqs.get_seq(ref_name)

    seq_scores = get_matching_score_to_human(seqs)

    quantile_score = get_score_at_quantile(seq_scores, quantile)
    seqs = seqs.take_seqs(ref_name, negate=True)
    seq_match_filtered = seqs.take_seqs_if(
        lambda seq: seq_scores[seq.name] > quantile_score
    )

    seqs = seq_match_filtered.add_seqs([ref])

    return seqs

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

@define_app
def remove_unresolvable_codon_seqeunce(seqs: UnalignedSeqsType) -> UnalignedSeqsType:
    valid_seqs = seqs.take_seqs_if(
        lambda seq: seq.get_translation()
    )
    return valid_seqs

drop_low_matching = low_matching_significance(quantile=0.05)
drop_ambiguous = too_many_ambigs(frac=0.02)
drop_invalid_length = length_divisible_by_three()
drop_short_seq = short_seq(length=500)
seq_without_redundent = remove_redundent_seq()
trim_stop = get_app("trim_stop_codons", gc=1)
remove_untranlatable_seqs = remove_unresolvable_codon_seqeunce()



@click.command()
@click.argument('input', type=click.Path(exists=True))
@click.option('-o', '--output_dir', help='Output directory for JSON files')
@click.option('-n', '--num_processes', default=multiprocessing.cpu_count(), help='Number of processes to use (default: number of CPUs)')
@click.option('-l', '--limit', type=int, help='Limit the number of files to process')
def main(input, num_processes, output_dir, limit):
    input_data_store = open_data_store(input, suffix='fa', limit=limit)
    out_dstore = open_data_store(output_dir, mode="w", suffix="fa")
    writer = get_app("write_seqs", out_dstore, format="fasta")
    loader = get_app("load_unaligned", format="fasta", moltype="dna")
    seq_filter = loader + valid_cds + drop_low_matching + drop_short_seq + drop_invalid_length + drop_ambiguous + seq_without_redundent + trim_stop + writer
    seq_filter.apply_to(input_data_store.completed, parallel=True, par_kw=dict(max_workers=num_processes))


if __name__ == "__main__":
    main()













# def pairwise_jsd_matrix(species_data):
#     species_keys = list(species_data.keys())
#     num_species = len(species_keys)
#     jsd_matrix = np.zeros((num_species, num_species))  # Initialize a square matrix

#     for i, species_1 in enumerate(species_keys):
#         for j, species_2 in enumerate(species_keys):
#             if i < j:  # To avoid recomputation, calculate only for i < j
#                 jsd_value = jsd(species_data[species_1], species_data[species_2])
#                 jsd_matrix[i, j] = jsd_value
#                 jsd_matrix[j, i] = jsd_value  # JSD is symmetric

#     return jsd_matrix



# def extract_info(path):
#     match = re.search(r'/([^/]+)\.fasta$', path)
#     if match:
#         return match.group(1)
#     else:
#         return "unknown"

#def jsd_genetic_distance_scatters(dataframes_dict):
#     keys = list(dataframes_dict.keys())
#     rows = int(len(keys) ** 0.5) + 1  # Calculate the number of rows for subplots
#     cols = (len(keys) + rows - 1) // rows  # Calculate the number of columns

#     fig = make_subplots(rows=rows, cols=cols, subplot_titles=[f'{key}' for key in keys])
    
#     # Populate subplots
#     for index, (key, df) in enumerate(dataframes_dict.items(), start=1):
#         row = (index - 1) // cols + 1
#         col = (index - 1) % cols + 1
        
#         fig.add_trace(
#             go.Scatter(
#                 x=np.log10(df['JSD Value']),
#                 y=df['Distance Value'],
#                 mode='markers'
#             ),
#             row=row,
#             col=col
#         )
#         # Only add x-axis title to bottom row plots

#         fig.update_xaxes(row=row, col=col, range =[-6, 0])
#         fig.update_yaxes(row=row, col=col, range =[0, 1.1])

#         if row == rows:
#             fig.update_xaxes(title_text="Log(JSD Value)", row=row, col=col, range=[-6, 0])
        
#         # Only add y-axis title to first column plots
#         if col == 1:
#             fig.update_yaxes(title_text="Genetic Distance", row=row, col=col)
    
#     fig.update_layout(
#         height=300 * rows,  # Set a reasonable height based on number of rows
#         width=300 * cols,   # Set a reasonable width based on number of columns
#         title_text="Scatter Plots of 3rd codon position JSD Vs Genetic Distance",
#         showlegend=False
#     )
    
#     return fig


# def prepare_dataframes(info_dict):
#     dataframes_dict = {}
#     for key, value in info_dict.items():
#         data = {
#             'JSD Value': value['pairwise_jsd'],
#             'Distance Value': value['pairwise_distance'],
#             'Correlation Factor': value['col'],
#             'P_value':value['p_value'],
#         }
#         dataframes_dict[key] = pd.DataFrame(data)
#     return dataframes_dict


# @define_app
# def align_via_aa(seqs: typing.SeqsCollectionType, gc=1) -> typing.AlignedSeqsType:
#     """Translates a nucleotide align amino acid sequences back to DNA"""
#     translater = get_app("translate_seqs", gc=gc)
#     prot_aln = get_app("progressive_align", "protein", unique_guides=True)
#     app = translater + prot_aln
#     aligned_aa = app(seqs).to_type(array_align=True)
#     return aligned_aa.replace_seqs(seqs)