import json
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.figure_factory as ff
import scipy
from cogent3.maths.measure import jsd
from scipy import stats

from plot_utils.util import update_figure_format


def load_json_data(path):
    with open(path, "r") as file:
        return json.load(file)


def get_ingroup_names(triads_info):
    ingroup_names_dict = {}
    for identifier, info in triads_info.items():
        triads_names = info["triples_species_names"]
        ingroup_names_dict[identifier] = [
            triads_names["ingroup1"],
            triads_names["ingroup2"],
        ]
    return ingroup_names_dict


def get_jsd_diff(triads_info):
    jsd_diff_dict = {}
    for identifier, info in triads_info.items():
        triads_info_value = info["triples_info_small_tree"]
        nuc_freqs_dict = triads_info_value["nuc_freqs_dict"]
        nuc_freq1 = nuc_freqs_dict["ingroup1"]
        nuc_freq2 = nuc_freqs_dict["ingroup2"]
        nuc_freq_internal_node = nuc_freqs_dict["internal_node"]
        jsd1 = jsd(nuc_freq1, nuc_freq_internal_node)
        jsd2 = jsd(nuc_freq2, nuc_freq_internal_node)
        jsd_diff = abs(jsd1 - jsd2)
        jsd_diff_dict[identifier] = jsd_diff
    return jsd_diff_dict


def get_jsd(triads_info):
    jsd_dict = {}
    for identifier, info in triads_info.items():
        triads_info_value = info["triples_info_small_tree"]
        nuc_freqs_dict = triads_info_value["nuc_freqs_dict"]
        nuc_freq1 = nuc_freqs_dict["ingroup1"]
        nuc_freq2 = nuc_freqs_dict["ingroup2"]
        nuc_freq_internal_node = nuc_freqs_dict["internal_node"]
        jsd1 = jsd(nuc_freq1, nuc_freq_internal_node)
        jsd2 = jsd(nuc_freq2, nuc_freq_internal_node)
        jsd_dict[identifier] = {"ingroup1": jsd1, "ingroup2": jsd2}
    return jsd_dict


# def get_ingroup_jsd(triads_info):
#     ingroup_jsd_dict = {}
#     for identifier, info in triads_info.items():
#         triads_info_value = info['triples_info_small_tree']
#         ingroup_jsd = triads_info_value['ingroup_jsd']
#         ingroup_jsd_dict[identifier] = ingroup_jsd
#     return ingroup_jsd_dict


def get_ingroup_ens_absdiff(triads_info):
    ens_ingroup_dict = {}
    for identifier, info in triads_info.items():
        triads_info_value = info["triples_info_small_tree"]
        triads_names = info["triples_species_names"]
        ens_dict = triads_info_value["ens"]
        ens_ingroup = abs(
            ens_dict[triads_names["ingroup1"]] - ens_dict[triads_names["ingroup2"]]
        )
        ens_ingroup_dict[identifier] = ens_ingroup
    return ens_ingroup_dict


def get_ens(triads_info):
    ens_value_dict = {}
    for identifier, info in triads_info.items():
        triads_info_value = info["triples_info_small_tree"]
        triads_names = info["triples_species_names"]
        ens_dict = triads_info_value["ens"]
        ens1 = ens_dict[triads_names["ingroup1"]]
        ens2 = ens_dict[triads_names["ingroup2"]]
        ens_value_dict[identifier] = {"ingroup1": ens1, "ingroup2": ens2}
    return ens_value_dict


# def get_nabla_absdiff(triads_info):
#     nabla_diff_dict = {}
#     for identifier, info in triads_info.items():
#         triads_info_value = info['triples_info_small_tree']
#         triads_names = info['triples_species_names']
#         nabla_dict = triads_info_value['nabla_values']
#         nbala_diff = abs(nabla_dict[triads_names['ingroup1']] - nabla_dict[triads_names['ingroup2']])
#         nabla_diff_dict[identifier] = nbala_diff
#     return nabla_diff_dict


def remove_outliers_iqr(data1, data2):
    def compute_iqr_bounds(data):
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 5 * IQR
        upper_bound = Q3 + 5 * IQR
        return lower_bound, upper_bound

    # Calculate IQR bounds for both lists
    lower_bound1, upper_bound1 = compute_iqr_bounds(data1)
    lower_bound2, upper_bound2 = compute_iqr_bounds(data2)

    # Filter out pairs where either value is an outlier
    filtered_data1 = []
    filtered_data2 = []

    for val1, val2 in zip(data1, data2):
        if (lower_bound1 <= val1 <= upper_bound1) and (
            lower_bound2 <= val2 <= upper_bound2
        ):
            filtered_data1.append(val1)
            filtered_data2.append(val2)

    return filtered_data1, filtered_data2


def get_names(triads_info):
    names_dict = {}
    for identifier, info in triads_info.items():
        triads_names = info["triples_species_names"]
        names_dict[identifier] = triads_names
    return names_dict


def remove_repested_names(gene_paths):
    species_names = {}
    for path in gene_paths:
        gene_name = os.path.basename(path.rstrip("/"))
        triads_data_path = os.path.join(path, "triples_info_dict_new.json")
        triads_info = load_json_data(triads_data_path)
        ingroup_names = get_ingroup_names(triads_info)
        species_names[gene_name] = ingroup_names

    repeated_names_dict = {}
    for gene, species in species_names.items():
        repeated_names_dict[gene] = []  # Initialize a list to store pairs
        for identifier, ingroup in species.items():
            for identifier2, ingroup2 in species.items():
                if identifier < identifier2:  # Ensure each pair is only processed once
                    if ingroup == ingroup2:
                        repeated_names_dict[gene].append((identifier, identifier2))

    # Optional: Remove genes with no repeated pairs
    repeated_names_dict = {
        gene: pairs for gene, pairs in repeated_names_dict.items() if pairs
    }
    removed_identifier = {
        gene: [a[0] for a in repeated_names_dict[gene]]
        if gene in repeated_names_dict
        else []
        for gene in species_names
    }
    return removed_identifier


# Function to compute required values
def compute_values(path, removed_identifier):
    gene_name = os.path.basename(path.rstrip("/"))
    triads_data_path = os.path.join(path, "triples_info_dict_new.json")
    triads_info_original = load_json_data(triads_data_path)
    triads_info = {
        k: v
        for k, v in triads_info_original.items()
        if k not in removed_identifier[gene_name]
    }
    ens_abs_diff_dict = get_ingroup_ens_absdiff(triads_info)
    jsd_diff_dict = get_jsd_diff(triads_info)
    ens_dict = get_ens(triads_info)
    jsd_dict = get_jsd(triads_info)
    species_names_dict = get_names(triads_info)
    ens_abs_diff_list = list(ens_abs_diff_dict.values())
    jsd_diff_list = list(jsd_diff_dict.values())
    ens_list = list(ens_dict.values())
    jsd_list = list(jsd_dict.values())
    species_names_list = list(species_names_dict.values())

    return ens_abs_diff_list, jsd_diff_list, ens_list, jsd_list, species_names_list


def get_gene_value_dict(gene_paths):
    gene_data_dict = {}
    # Populate the dictionary with data for each gene
    removed_identifier = remove_repested_names(gene_paths)
    for path in gene_paths:
        gene_name = os.path.basename(path.rstrip("/"))
        ens_abs_diff_list, jsd_diff_list, ens_list, jsd_list, species_names_list = (
            compute_values(path, removed_identifier)
        )
        gene_data_dict[gene_name] = {
            "ens_abs_diff": ens_abs_diff_list,
            "jsd_diff": jsd_diff_list,
            "ens": ens_list,
            "jsd": jsd_list,
            "species_names": species_names_list,
        }
    return gene_data_dict


def spearman_correlation_analysis(gene_data_dict):
    species_number_dict = {}
    max_group_size_dict = {}
    num_group_dict = {}

    for gene, value in gene_data_dict.items():
        data_f = pd.DataFrame(
            {
                "ens_abs_diff": np.sqrt(value["ens_abs_diff"]),
                "jsd_diff": np.sqrt(value["jsd_diff"]),
                "Species1": [x["ingroup1"] for x in value["species_names"]],
                "Species2": [x["ingroup2"] for x in value["species_names"]],
                "Species3": [x["outgroup"] for x in value["species_names"]],
            }
        )

        data_long = pd.melt(
            data_f,
            id_vars=["ens_abs_diff", "jsd_diff"],
            value_vars=["Species1", "Species2", "Species3"],
            var_name="Species_Position",
            value_name="Species",
        )
        data_long["ens_abs_diff"] = data_long["ens_abs_diff"]
        data_long["jsd_diff"] = data_long["jsd_diff"]
        data_long["Species"] = data_long["Species"].astype(str)
        data_long["Species"] = data_long["Species"].astype("category")

        num_groups = len(set(data_long["Species"]))
        species_number_dict[gene] = len(set(data_long["Species"]))
        num_group_dict[gene] = num_groups
        max_group_size = max(data_long.groupby("Species", observed=False).size())
        max_group_size_dict[gene] = max_group_size

    correlation_list = {}
    p_value_list = {}
    for gene, lists in gene_data_dict.items():
        jsd_diff_list, ens_abs_diff_list = remove_outliers_iqr(
            lists["jsd_diff"], lists["ens_abs_diff"]
        )

        # Add the correlation factor in the list
        cor, p_value = scipy.stats.spearmanr(jsd_diff_list, ens_abs_diff_list)
        correlation_list[gene] = cor
        p_value_list[gene] = p_value

    # Step 1: Correct p-values using the group size
    p_value_list_corrected = {
        gene: p_value_list[gene] * max_group_size_dict[gene]
        for gene in gene_data_dict.keys()
    }

    return p_value_list, p_value_list_corrected, correlation_list


# Ploting functions


def get_correlation_factor_histogram(correlation_list):

    # Data for histogram
    values = list(correlation_list.values())
    # Create the histogram with density normalization
    fig = px.histogram(
        values,
        labels={"x": "Correlation Coefficient", "y": "Count"},
        title=None,
        color_discrete_sequence=["#67a8cd"],  # Set the color to a shade of orange
    )

    # Update layout for presentation
    fig.update_layout(
        yaxis_title="<b>Count</b>", xaxis_title=r"$\Large \hat{\rho}$", showlegend=False
    )

    fig.add_shape(
        type="line",
        x0=0.2,
        y0=0,
        x1=0.2,
        y1=15,
        line=dict(color="#c73d47", width=4, dash="dashdot"),
    )

    # Set transparency level and add a solid line around each bar
    fig.update_traces(
        opacity=1,  # Set the transparency (0 = fully transparent, 1 = fully opaque)
        marker_line_color="black",  # Color of the line around each bar
        marker_line_width=1.5,  # Width of the line around each bar
        xbins=dict(size=0.05),
    )

    fig = update_figure_format(fig)

    return fig


def get_correlation_factors_distplot(correlation_list):
    hist_data = [list(correlation_list.values())]
    group_labels = ["Correlation Coefficient"]

    fig = ff.create_distplot(
        hist_data, group_labels, bin_size=0.02, show_rug=False, colors=["#67a8cd"]
    )

    fig.update_layout(
        xaxis_title="Correlation Coefficient",
        yaxis_title="count",
        title=None,
        showlegend=False,
    )

    fig = update_figure_format(fig)

    fig.update_layout(
        xaxis=dict(
            title=r"$\Large \hat{\rho} (\delta_{\text{JAD}}, \delta_{\mathbb{t}})$",
        ),
        yaxis=dict(title="Density"),
        width=800,
        height=600,
    )

    fig.update_yaxes(title_font=dict(size=30, family="Times New Roman", color="black"))

    fig.update_xaxes(title_font=dict(color="black"))

    return fig


def get_correlation_factor_tos_rejection_correlation_fig(
    proportion_less_than_005_tos, corr_list
):
    corr_list_filtered = {
        gene: corr_list[gene] for gene in proportion_less_than_005_tos.keys()
    }

    fig = px.scatter(
        x=np.array(list(proportion_less_than_005_tos.values())) * 100,
        y=np.array(list(corr_list_filtered.values())),
        labels={
            "x": "Proportion reject the stationarity",
            "y": r"$\text{Spearman } \rho$",
        },
        trendline="ols",
        title=None,
    )

    fig.update_layout(
        xaxis_title="Stationarity rejected %",
        yaxis_title=r"$\Large \hat{\rho}$",
        margin=dict(b=160),
    )

    # Customize the markers (Change the color to #7ed3f6)
    fig.update_traces(
        marker=dict(
            size=8,
            opacity=0.8,
            color="#67a8cd",
        ),
        selector=dict(mode="markers"),
    )

    # Customize the trendline
    fig.update_traces(selector=dict(mode="lines"), line=dict(color="black", width=2))

    cor, p = stats.spearmanr(
        list(proportion_less_than_005_tos.values()), list(corr_list_filtered.values())
    )

    # fig.add_annotation(
    #     x=1,  # Adjust x-position of the annotation (relative to x-axis range)
    #     y=1,  # Adjust y-position of the annotation (relative to y-axis range)
    #     text = rf'$\hat{{\rho}} = {cor:.4f}$',  # R² value with 4 decimal precision
    #     showarrow=False,
    #     font=dict(size=14, color="red"),
    #     xref="paper",  # Positioning relative to the paper (can be set to 'x' and 'y' for data-coordinates)
    #     yref="paper",  # Positioning relative to the paper (can be set to 'x' and 'y' for data-coordinates)
    #     bgcolor="rgba(255, 255, 255, 0.7)",  # Background color for better visibility
    #     bordercolor="black",  # Optional: border around the annotation
    #     borderwidth=1
    # )

    fig = update_figure_format(fig)

    fig.update_layout(
        xaxis=dict(
            title="Stationarity rejected%",
        ),
        yaxis=dict(
            title=r"$\huge \hat{\rho} (\delta_{\text{JAD}}, \delta_{\mathbb{t}})$"
        ),
        title_font=dict(size=50, family="Times New Roman", color="black"),
        width=800,
        height=600,
        margin=dict(l=160),
    )

    fig.update_xaxes(title_font=dict(size=25, family="Times New Roman", color="black"))

    fig.update_yaxes(title_font=dict(color="black"))

    return fig
