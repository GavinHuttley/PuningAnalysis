from cogent3 import get_app
from cogent3.app.composable import NotCompleted
import numpy as np
import plotly.graph_objects as go
from scipy.stats import uniform
from collections import defaultdict
from random import sample
import plotly.express as px
from scipy.stats import spearmanr



load_json_app = get_app("load_json")


def p_value_ST(result):
    return sum(result.observed.LR <= null_lr for null_lr in result.null_dist) / len(result.null_dist)

def p_value_null_distirbution(result, value): 
    return sum(value <= null_lr for null_lr in result.null_dist) / len(result.null_dist)



def qq_plot_uniform(data, a=0, b=1):
    """
    Creates a QQ plot of the provided data against a uniform distribution using Plotly.

    Args:
    data (array-like): The dataset to plot.
    a (float): The lower bound of the uniform distribution (default 0).
    b (float): The upper bound of the uniform distribution (default 1).
    """
    # Scale data for the specified uniform range
    data = np.array(data)
    data.sort()  # Sort the data for plotting
    scaled_data = (data - a) / (b - a)

    # Calculate quantiles
    n = len(data)
    theoretical_quantiles = uniform.ppf(np.arange(1, n + 1) / (n + 1))

    # Create a QQ plot
    fig = go.Figure()

    # Adding scatter plot for QQ plot
    fig.add_trace(go.Scatter(x=theoretical_quantiles, y=scaled_data, mode='markers',
                             name=None,
                             showlegend=False,
                             marker=dict(color='#67a8cd', size = 5)))

    # Add line of perfect fit
    fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines',
                             name=None,
                             showlegend=False,
                             line=dict(color='red', dash='dash')))
    
    
    fig.update_layout(title=None, 
                  xaxis_title='Uniform Quantiles', 
                  yaxis_title=r'$\hat{p}-\text{value}$',
                    template='plotly_white',
                    width=500, height=600,
                    margin=dict(l=20, r=20, t=50, b=20),
                    autosize=True,
                    yaxis_title_font={'size': 25},  
                    xaxis_title_font={'size': 25},
                    showlegend=True,
                        xaxis=dict(
            title_font=dict(size=25, family='Arial', color='black'),
            tickfont=dict(size=14),
            gridcolor='lightgrey'
        ),
                    yaxis=dict(
            title_font=dict(size=25, family='Arial', color='black'),
            tickfont=dict(size=14),
            gridcolor='lightgrey')
    )

    return fig


def qq_plot_null_observed(data, data2, a=0, b=1):
    """
    Creates a QQ plot of the provided data against a uniform distribution using Plotly.

    Args:
    data (array-like): The dataset to plot.
    a (float): The lower bound of the uniform distribution (default 0).
    b (float): The upper bound of the uniform distribution (default 1).
    """
    # Scale data for the specified uniform range
    data = np.array(data)
    data.sort()  # Sort the data for plotting
    scaled_data = (data - a) / (b - a)

    # Calculate quantiles
    n = len(data)
    theoretical_quantiles = uniform.ppf(np.arange(1, n + 1) / (n + 1))

    # Scale data for the specified uniform range
    data2 = np.array(data2)
    data2.sort()  # Sort the data for plotting
    scaled_data2 = (data2 - a) / (b - a)

    # Calculate quantiles
    n2 = len(data2)
    theoretical_quantiles2 = uniform.ppf(np.arange(1, n2 + 1) / (n2 + 1))

    # Create a QQ plot
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=theoretical_quantiles, y=scaled_data, mode='markers',
                                name='<b>-ve</b>',
                                marker=dict(color='#67a8cd', size = 4)))

    # Adding scatter plot for QQ plot
    fig.add_trace(go.Scatter(x=theoretical_quantiles2, y=scaled_data2, mode='markers',
                                name='<b>Observed</b>',
                                marker=dict(color='#6fba4f', size = 4)))

    # # Add line of perfect fit
    # fig.add_trace(go.Scatter(x=[0, 1], y=[0, 1], mode='lines',
    #                             name=None,
    #                             line=dict(color='red', dash='dash')))

    fig.update_layout(title=None, xaxis_title='Uniform Quantiles', yaxis_title=r'$\hat{p}-\text{value}$',
                    template='plotly_white',
                    margin=dict(l=20, r=20, t=50, b=20),
                    autosize=True,
                    yaxis_title_font={'size': 25, 'family': 'Arial', 'color': 'black'},  
                    xaxis_title_font={'size': 25}, 
                    width=600, height=600,
                    showlegend=True,
                    legend=dict(font = dict(size=16)),
            xaxis=dict(
            title_font=dict(size=25, family='Arial', color='black'),
            tickfont=dict(size=14),
            gridcolor='lightgrey'
        ),
                    yaxis=dict(
            title_font=dict(size=25, family='Arial', color='black'),
            tickfont=dict(size=14),
            gridcolor='lightgrey'))
    
    
    return fig


def get_proportion_rejected_correlation_fig(proportion_less_than_005_tos, proportion_less_than_005_toc):
    proportion_less_than_005_clock_filtered = {gene: proportion_less_than_005_toc[gene] for gene in proportion_less_than_005_tos.keys()}
    x_values = np.array(list(proportion_less_than_005_tos.values())) * 100
    y_values = np.array(list(proportion_less_than_005_clock_filtered.values())) * 100

    # Create scatter plot
    fig_0 = px.scatter(
        x=x_values,
        y=y_values,
        labels={'x': 'Proportion reject the clock', 'y': 'Proportion reject the stationarity'},
        trendline="ols",
        title=None
    )

    # Update layout for axis titles and legend
    fig_0.update_layout(
        xaxis=dict(
            title='<b>Stationarity rejected %</b>',
        ),
        yaxis=dict(
            title='<b>Clock rejected %</b>',
        )
    )

    # Update traces for markers and trendline
    fig_0.update_traces(
        marker=dict(
            size=8,
            opacity=0.8,
            color='#7ed3f6',
            line=dict(width=1, color='DarkSlateGrey')
        ),
        selector=dict(type='scatter', mode='markers')
    )
    fig_0.update_traces(
        line=dict(color='#c73d47', width=2),
        selector=dict(type='scatter', mode='lines')
    )


    # Calculate Spearman correlation and add annotation
    corr, p = spearmanr(list(proportion_less_than_005_tos.values()), list(proportion_less_than_005_toc.values()))
    fig_0.add_annotation(
        x=1, y=1,
        text=f'ρ = {corr:.4f}',  # Spearman's rho with 4 decimal precision
        showarrow=False,
        font=dict(size=15, color="red"),
        xref="paper", yref="paper",
        bgcolor="rgba(255, 255, 255, 0.7)",
        bordercolor="black",
        borderwidth=1
    )

    return fig_0

def get_p_value_distirbution(input_data_dir):
    p_value_observed_dict = {}
    p_value_tested_dict = {}
    i = 0
    for data in input_data_dir:
        gene_name = data.unique_id.split('.')[0]
        bootstrap_result = load_json_app(data)
        if isinstance(bootstrap_result, NotCompleted):
            i += 1
            print(data) 
        else:
            p_value_observed_dict[gene_name] = bootstrap_result.observed.pvalue
            p_value_tested_dict[gene_name] = p_value_ST(bootstrap_result)

    # get the p-value of null simulated for each algnment 
    #get distirbution of null for tos
    p_value_list_tos_null = {}
    for path in input_data_dir:
        gene_name = path.unique_id.split('.')[0]
        bootstrap_result = load_json_app(path)
        if not isinstance(bootstrap_result, NotCompleted):
            pvalue = None
            # Keep sampling until a non-None pvalue is found
            while pvalue is None:
                null_sample_key = sample(list(bootstrap_result.keys())[1:], 1)  # Sampling from keys, skip first key
                sub_result = bootstrap_result[null_sample_key[0]]
                pvalue = sub_result.pvalue  # Get the pvalue from the sub_result
            p_value = p_value_null_distirbution(bootstrap_result, sub_result.LR)
            p_value_list_tos_null[gene_name] = p_value
    return p_value_observed_dict, p_value_tested_dict, p_value_list_tos_null



def get_rejected_proportion(p_value_observed_dict):
    
    # Dictionary to store the count of values < 0.05 for each gene
    count_less_than_005_clock = defaultdict(int)
    # Dictionary to store the total count of entries for each gene
    total_count = defaultdict(int)

    for key, value in p_value_observed_dict.items():
        gene = key.split('_')[0]
        total_count[gene] += 1
        if value != None:
            if value < 0.05:
                count_less_than_005_clock[gene] += 1

    # Calculate the proportion of values < 0.05 for each gene
    proportion_rejected = {gene: count_less_than_005_clock[gene] / total for gene, total in total_count.items()}
    return proportion_rejected
