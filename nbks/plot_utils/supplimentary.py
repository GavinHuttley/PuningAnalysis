import json
import os

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import statsmodels.formula.api as smf
from cogent3 import get_app
from plotly.subplots import make_subplots

from plot_utils.util import update_figure_format

# Mammal ortholog data alignment length and sample size distirbution

loader = get_app("load_aligned", format="fasta", moltype="dna")


def get_alignment_length_and_sample_size_distirbution(result_dstore):
    alignment_length_dict = {}
    sample_size_dict = {}
    for path in result_dstore.completed:
        file_name = path.unique_id
        aln = loader(path)
        alignment_length = aln.seq_len
        sample_size = aln.num_seqs
        alignment_length_dict[file_name] = alignment_length
        sample_size_dict[file_name] = sample_size

    return alignment_length_dict, sample_size_dict


def get_distirbution_plot(data_dict, x_title, y_title):
    values = list(data_dict.values())
    fig = px.histogram(
        values,
        labels={"x": x_title, "y": y_title},
        title=None,
        color_discrete_sequence=["#67a8cd"],  # Set the color to a shade of orange
        histnorm="density",  # Normalize the histogram to density
    )

    fig.update_layout(
        template="plotly_white",
        margin=dict(l=50, r=50, t=50, b=50),  # Adjust margins for a balanced look
        autosize=True,
        yaxis_title=y_title,  # Explicit y-axis title
        xaxis_title=x_title,  # Explicit x-axis title
        yaxis_title_font=dict(size=18),  # Adjust y-axis font size
        xaxis_title_font=dict(size=18),  # Adjust x-axis font size
        font=dict(size=16),  # General font size for labels and titles
        width=800,  # Set figure width (optional for better control)
        height=500,  # Set figure height (optional for better control)
        showlegend=False,  # Remove the legend
    )

    fig.update_traces(
        opacity=0.7,  # Set the transparency (0 = fully transparent, 1 = fully opaque)
        marker_line_color="black",  # Color of the line around each bar
        marker_line_width=1.5,  # Width of the line around each bar
    )

    fig.update_xaxes(range=[450, max(values) + 50])  # Show the figure

    fig = update_figure_format(fig)

    return fig


# rejected proportion of TOS and TOC


def get_rehected_proportion(p_value_observed_dict):
    from collections import defaultdict

    # Dictionary to store the count of values < 0.05 for each gene
    rejected_cases_count = defaultdict(int)
    # Dictionary to store the total count of entries for each gene
    total_count = defaultdict(int)

    for key, value in p_value_observed_dict.items():
        gene = key.split("_")[0]
        total_count[gene] += 1
        if float(value) < 0.05:
            rejected_cases_count[gene] += 1

    # Calculate the proportion of values < 0.05 for each gene
    proportion_less_than_005_toc = {
        gene: rejected_cases_count[gene] / total for gene, total in total_count.items()
    }
    return proportion_less_than_005_toc


def get_rejected_proportion_plot(rejected_proportion_tos, rejected_proportion_toc):

    # Create the histogram with count normalization (default behavior)
    fig = px.histogram(
        rejected_proportion_tos.values(),
        labels={
            "x": "Rejected proportion",
            "y": "Count",
        },  # Label for y-axis as 'Count'
        title=None,
        color_discrete_sequence=["#67a8cd"],  # Set the color to a shade of blue
    )

    # Update layout for presentation
    fig.update_layout(
        template="plotly_white",
        margin=dict(l=50, r=50, t=50, b=50),  # Adjust margins for a balanced look
        autosize=True,
        yaxis_title="<b>Count</b>",  # Explicit y-axis title
        xaxis_title="<b>Rejected proportion</b>",  # Explicit x-axis title
        yaxis_title_font=dict(size=20),  # Adjust y-axis font size
        xaxis_title_font=dict(size=20),  # Adjust x-axis font size
        font=dict(size=16),  # General font size for labels and titles
        width=800,  # Set figure width (optional for better control)
        height=500,  # Set figure height (optional for better control)
        showlegend=False,  # Remove the legend
    )
    fig.add_shape(
        type="line",
        x0=0.3,
        y0=0,
        x1=0.3,
        y1=35,
        line=dict(color="#c73d47", width=4, dash="dashdot"),
    )

    # Set transparency level and add a solid line around each bar
    fig.update_traces(
        opacity=0.8,  # Set the transparency (0 = fully transparent, 1 = fully opaque)
        marker_line_color="black",  # Color of the line around each bar
        marker_line_width=1.5,  # Width of the line around each bar
    )

    fig2 = px.histogram(
        rejected_proportion_toc.values(),
        labels={
            "x": "Rejected proportion",
            "y": "Count",
        },  # Label for y-axis as 'Count'
        title=None,
        color_discrete_sequence=["#67a8cd"],  # Set the color to a shade of blue
    )

    # Update layout for presentation
    fig2.update_layout(
        template="plotly_white",
        margin=dict(l=50, r=50, t=50, b=50),  # Adjust margins for a balanced look
        autosize=True,
        yaxis_title="<b>Count</b>",  # Explicit y-axis title
        xaxis_title="<b>Rejected proportion</b>",  # Explicit x-axis title
        yaxis_title_font=dict(size=20),  # Adjust y-axis font size
        xaxis_title_font=dict(size=20),  # Adjust x-axis font size
        font=dict(size=16),  # General font size for labels and titles
        width=800,  # Set figure width (optional for better control)
        height=500,  # Set figure height (optional for better control)
        showlegend=False,  # Remove the legend
    )
    fig2.add_shape(
        type="line",
        x0=0.3,
        y0=0,
        x1=0.3,
        y1=35,
        line=dict(color="#c73d47", width=4, dash="dashdot"),
    )

    # Set transparency level and add a solid line around each bar
    fig2.update_traces(
        opacity=0.8,  # Set the transparency (0 = fully transparent, 1 = fully opaque)
        marker_line_color="black",  # Color of the line around each bar
        marker_line_width=1.5,  # Width of the line around each bar
    )

    return fig


def get_multiple_linear_regression(gene_data_dict):
    """
    Perform multiple linear regression on the provided gene data dictionary.
    The dictionary should contain 'ens_abs_diff', 'jsd_diff', and 'species_names' for each gene.
    Returns a DataFrame with regression results.
    """

    multiple_linear_regression_data = {
        "gene": {},
        "p_value": {},
        "marginal_r^2": {},
        "conditional_r^2": {},
    }
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

        # Step 3: Adjust variables
        data_long["ens_abs_diff"] = data_long["ens_abs_diff"] / 3
        data_long["jsd_diff"] = data_long["jsd_diff"] / 3

        # Step 4: Convert species identifiers to strings and factors
        data_long["Species"] = data_long["Species"].astype(str)
        data_long["Species"] = data_long["Species"].astype("category")

        model = smf.mixedlm(
            "ens_abs_diff ~ jsd_diff", data=data_long, groups=data_long["Species"]
        )
        result = model.fit()

        # Random effects variance
        var_random = result.cov_re.iloc[0, 0]

        # Fixed effects variance
        # Calculate variance of the linear predictor (fixed effects)
        fixed_effects = result.fe_params
        X = result.model.exog
        var_fixed = np.var(np.dot(X, fixed_effects))

        # Residual variance
        var_residual = result.scale

        # Marginal R-squared (fixed effects only)
        R_m2 = var_fixed / (var_fixed + var_random + var_residual)

        # Conditional R-squared (fixed + random effects)
        R_c2 = (var_fixed + var_random) / (var_fixed + var_random + var_residual)

        multiple_linear_regression_data["gene"][gene] = gene
        multiple_linear_regression_data["p_value"][gene] = result.pvalues["jsd_diff"]
        multiple_linear_regression_data["marginal_r^2"][gene] = R_m2
        multiple_linear_regression_data["conditional_r^2"][gene] = R_c2

    multiple_linear_regression_data["r^2_difference"] = {
        gene: multiple_linear_regression_data["conditional_r^2"][gene]
        - multiple_linear_regression_data["marginal_r^2"][gene]
        for gene in multiple_linear_regression_data["gene"]
    }
    r_squared_data = pd.DataFrame(multiple_linear_regression_data)

    return r_squared_data


def plot_r2_by_gene(df):
    """
    df: pandas.DataFrame with columns 'gene', 'marginal_r^2', 'conditional_r^2'
    """

    genes = df["gene"].astype(str)  # categorical x-axis

    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=genes,
            y=df["marginal_r^2"],
            mode="lines+markers",
            name="Marginal R²",
            marker=dict(color="#67a8cd", size=4),
            line=dict(width=1),
        )
    )
    fig.add_trace(
        go.Scatter(
            x=genes,
            y=df["conditional_r^2"],
            mode="lines+markers",
            name="Conditional R²",
            marker=dict(color="#6fba4f", size=4),
            line=dict(width=1),
        )
    )

    fig.update_layout(
        xaxis=dict(title="Gene", tickangle=-45),
        yaxis=dict(title="R²", range=[0, 1]),
        template="plotly_white",
        height=500,
        width=1000,
        margin=dict(b=140),  # room for rotated gene labels
        legend=dict(
            title=None,
            x=0.5,  # Center the legend
            y=1.1,
            xanchor="center",
            yanchor="top",  # Anchor the legend to the top to ensure it's below the plot
            orientation="h",  # Horizontal orientation
        ),
    )

    fig = update_figure_format(fig)
    return fig


# Get ENS difference and Ingroup JSD distirbution in selected triples


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


def get_triple_ens_diff_distirbution_plot(path):

    gene_name = str(path).split("/")[-1]
    triads_data_path = os.path.join(path, "triples_info_dict_new.json")
    with open(triads_data_path, "r") as f:
        triads_info = json.load(f)
    ens_ingroup_list = []
    ingroup_jsd_list = []
    for _, info in triads_info.items():
        triads_info_value = info["triples_info_small_tree"]
        triads_names = info["triples_species_names"]
        ingroup_jsd = triads_info_value["ingroup_jsd"]
        ens_dict = triads_info_value["ens"]
        ens_ingroup = abs(
            ens_dict[triads_names["ingroup1"]] - ens_dict[triads_names["ingroup2"]]
        )
        ens_ingroup_list.append(ens_ingroup)
        ingroup_jsd_list.append(ingroup_jsd)

    ens_ingroup_list2, ingroup_jsd_list2 = remove_outliers_iqr(
        ens_ingroup_list, ingroup_jsd_list
    )

    # sort each list for clearer visualization and create indices
    ens_sorted = sorted(ens_ingroup_list2)
    jsd_sorted = sorted(ingroup_jsd_list2)
    idx_ens = list(range(1, len(ens_sorted) + 1))
    idx_jsd = list(range(1, len(jsd_sorted) + 1))

    # create 2-row vertical subplots
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=False,
        vertical_spacing=0.08,
        subplot_titles=(None, None),  # we will use overall title instead
    )

    # top: ENS difference
    fig.add_trace(
        go.Scatter(
            x=idx_ens,
            y=ens_sorted,
            mode="markers",
            marker=dict(size=4, color="#6fba4f"),
            name="ENS difference",
        ),
        row=1,
        col=1,
    )

    # bottom: Ingroup JSD
    fig.add_trace(
        go.Scatter(
            x=idx_jsd,
            y=jsd_sorted,
            mode="markers",
            marker=dict(size=4, color="#67a8cd"),
            name="Ingroup JSD",
        ),
        row=2,
        col=1,
    )

    # axes labels
    fig.update_xaxes(title_text="<b>Index</b>", row=2, col=1)
    fig.update_yaxes(title_text="<b>ENS difference</b>", row=1, col=1)
    fig.update_yaxes(title_text="<b>Ingroup JSD</b>", row=2, col=1)

    # layout and sizing (tweak as needed)
    fig.update_layout(
        title=gene_name,
        showlegend=False,
        width=600,
        height=650,
    )

    fig = update_figure_format(fig)

    return fig
