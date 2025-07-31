from cogent3 import get_app, open_data_store
from cogent3.app.composable import NotCompleted
import numpy as np
import plotly.graph_objects as go
from scipy.stats import uniform

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


def get_p_value_distirbution(input_data_dir):
    tos_result_dir = '/Users/gulugulu/clock/mammal_orthologs_hsap_1/tos_boostrapping_results'
    input_data_store_tos_result = open_data_store(tos_result_dir, suffix= 'json')

    p_value_observed_dict = {}
    p_value_tested_dict = {}
    i = 0
    for data in input_data_store_tos_result:
        gene_name = data.unique_id.split('.')[0]
        bootstrap_result = load_json_app(data)
        if isinstance(bootstrap_result, NotCompleted):
            i += 1
            print(data) 
        else:
            p_value_observed_dict[gene_name] = bootstrap_result.observed.pvalue
            p_value_tested_dict[gene_name] = p_value_ST(bootstrap_result)

    # get the p-value of null simulated for each algnment 
    from random import sample
    #get distirbution of null for tos
    p_value_list_tos_null = {}
    for path in input_data_store_tos_result:
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