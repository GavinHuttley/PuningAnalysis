import numpy as np
from cogent3.maths.measure import jsd
from clock_project.simulation.magnitude_quantification import calculate_ENS
from cogent3 import get_app
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import statsmodels.api as sm
from scipy.stats import spearmanr
import scipy
import scipy.linalg
from plot_utils.util import update_figure_format


def compute_aginst_Q(matrices_list, pi, t_range):
    result_list = []
    for i in range(len(matrices_list)):
        q_pair = list(matrices_list[i].values())
        Q1 = np.array(q_pair[0])
        Q2 = np.array(q_pair[1])
        
        for t in t_range:
            jad_diff = get_jad_difference(pi, Q1, Q2, t)
            ens_abs_diff = get_ens_abs_diff(pi, Q1, Q2, t)
            result_list.append((i, t, ens_abs_diff, jad_diff))
    
    df = pd.DataFrame(result_list, columns=['Matrix_ID', 'Time', 'ENS_abs_difference', 'JSD_difference'])
    return df

def remove_outliers_iqr(data1, data2):
    def compute_iqr_bounds(data):
        Q1 = np.percentile(data, 25)
        Q3 = np.percentile(data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 6 * IQR
        upper_bound = Q3 + 6 * IQR
        return lower_bound, upper_bound

    # Calculate IQR bounds for both lists
    lower_bound1, upper_bound1 = compute_iqr_bounds(data1)
    lower_bound2, upper_bound2 = compute_iqr_bounds(data2)

    # Filter out pairs where either value is an outlier
    filtered_data1 = []
    filtered_data2 = []

    for val1, val2 in zip(data1, data2):
        if (lower_bound1 <= val1 <= upper_bound1) and (lower_bound2 <= val2 <= upper_bound2):
            filtered_data1.append(val1)
            filtered_data2.append(val2)

    return filtered_data1, filtered_data2



def get_ens_abs_diff(pi, Q1, Q2, t):
    ens1 = calculate_ENS(pi, Q1, t)
    ens2 = calculate_ENS(pi, Q2, t)
    ens_absdiff = abs(ens1-ens2)
    return ens_absdiff


def get_jad_difference(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    jsd_1 = jsd(pi_1, pi)
    jsd_2 = jsd(pi_2, pi)
    jsd_diff = abs(jsd_1 - jsd_2)
    return jsd_diff



def plot_time_grouped_scatter_2x2(df, x_col, y_col):
    times = sorted(df['Time'].unique())
    subplot_titles = [rf'$t_s={t}$' for t in times]
    fig = make_subplots(
        rows=2, cols=2,
        horizontal_spacing=0.15, vertical_spacing=0.2,
        subplot_titles=subplot_titles           # <-- AND THIS
    )


    subplot_positions = [(1, 1), (1, 2), (2, 1), (2, 2)]


    for position, time in zip(subplot_positions, times):
        sub_df = df[df['Time'] == time]
        x_data = sub_df[x_col]
        y_data = sub_df[y_col]

        x_data_no_outliers, y_data_no_outliers = remove_outliers_iqr(x_data, y_data)

        # Scatter plot for the actual data points
        fig.add_trace(go.Scatter(
            x=x_data_no_outliers, 
            y=y_data_no_outliers, 
            mode='markers', 
            name=time,  # This will automatically create legend entries
            marker=dict(color = '#67a8cd',
                size=1  # Adjust size as needed
            )),
            row=position[0], col=position[1])
        
        corr, _ = spearmanr(x_data_no_outliers, y_data_no_outliers)

        # Calculate OLS trendline
        x = sm.add_constant(x_data_no_outliers)  # adding a constant for OLS
        model = sm.OLS(y_data_no_outliers, x).fit()
        trendline = model.predict(x)

        # Add trendline trace
        fig.add_trace(go.Scatter(
            x=x_data_no_outliers, 
            y=trendline,
            mode='lines',
            line=dict(color='firebrick', width=2),  # Deep red for trendlines
            showlegend=False),  # Hide legend for trendline to avoid duplicate entries
            row=position[0], col=position[1])

        # Annotate R^2 value
        fig.add_annotation(
            xref="x domain", yref="y domain",
            x=0.95, y=0.95,
            text = rf'$\hat{{\rho}} = {corr:.4f}$', 
            showarrow=False,
            font=dict(size=10, color="red"),
            align="right",
            ax=0, ay=0,
            bordercolor="black",
            borderwidth=1,
            borderpad=4,
            bgcolor="white",
            opacity=0.8,
            row=position[0], col=position[1])
        
        
        # x_min, x_max = 0, round(np.max(x_data_no_outliers), 1) + 0.1  # keep real max
        # y_min, y_max = 0, round(np.max(y_data_no_outliers), 1) + 0.1

        # x_step = (x_max - x_min) / 5 if x_max > 0 else 0.1      
        # y_step = (y_max - y_min) / 5 if x_max > 0 else 0.1   
        
        # fig.update_xaxes(
        #     range=[x_min, x_max],   # force the range
        #     dtick=x_step,             # fixed spacing between ticks
        #     tickformat=".2f",       # show small decimals clearly
        #     showticklabels=True,    # ensure labels are not hidden in subplots
        #     ticks="outside",
        #     row=position[0], col=position[1]
        # )
                
        # fig.update_yaxes(
        #     range=[y_min, y_max],   # force the range
        #     dtick=y_step,             # fixed spacing between ticks
        #     tickformat=".2f",       # show small decimals clearly
        #     showticklabels=True,    # ensure labels are not hidden in subplots
        #     ticks="outside",
        #     row=position[0], col=position[1]
        # )


    fig.add_shape(
    type='line',
    x0=0.5, y0=0, x1=0.5, y1=1,
    line=dict(color="black", width=0.5),
    xref='paper', yref='paper',
    layer="above"
    )
    fig.add_shape(
        type='line',
        x0=0, y0=0.5, x1=1, y1=0.5,
        line=dict(color="black", width=0.5),
        xref='paper', yref='paper',
        layer="above"
    )

    fig.update_layout( 
        showlegend=False  
    )

    fig.update_traces(
    marker=dict(
        size=3,               
        opacity=1,          
        line=None  
    )
)
    fig.update_xaxes(
    title_text=r'\Large $\delta (JAD)$', 
    row=2
)
    fig.update_yaxes(
        title_text=r'$\Large \delta (ENS)$', 
        col=1
    )

    fig = update_figure_format(fig)


    return fig



def correlation_factor_plot(df_low, df_high, t_range): 

    grouped_low = df_low.groupby('Time')
    grouped_high = df_high.groupby('Time')

    correlation_results_low_jad = {}

    # Loop over each group and calculate the Pearson correlation
    for time, group in grouped_low:
        corr, p = spearmanr(group['JSD_difference'], group['ENS_abs_difference'])
        correlation_results_low_jad[time] = (corr,p)

    # Dictionary to store the correlation results
    correlation_results_high_jad = {}

    # Loop over each group and calculate the Pearson correlation
    for time, group in grouped_high:
        corr, p = spearmanr(group['JSD_difference'], group['ENS_abs_difference'])
        correlation_results_high_jad[time] = (corr, p)

    correlation_factors_low_jad = {time: correlation_results_low_jad[time][0] for time in correlation_results_low_jad.keys()}
    correlation_factors_high_jad = {time: correlation_results_high_jad[time][0] for time in correlation_results_high_jad.keys()}

    # Combine all data into a dataframe
    data = {
        'Time': t_range * 2,
        'Measure': ['JAD'] * 4 + ['JAD'] * 4,
        'High/Low': ['Balanced'] * 4 + ['Imbalanced'] * 4,
        'Correlation': list(correlation_factors_low_jad.values()) + list(correlation_factors_high_jad.values())
    }

    # Create DataFrame
    df_correlation = pd.DataFrame(data)

    # Create a new figure
    fig = go.Figure()


    # Add traces for each measure and condition (Balanced/Imbalanced)
    for measure in df_correlation['Measure'].unique():
        for condition in ['Balanced', 'Imbalanced']:
            filtered_data = df_correlation[(df_correlation['Measure'] == measure) & (df_correlation['High/Low'] == condition)]
            fig.add_trace(go.Scatter(
                x=filtered_data['Time'], 
                y=filtered_data['Correlation'], 
                mode='lines+markers', 
                name=f'{condition}',
                line=dict(dash='dash' if condition == 'Balanced' else 'solid', color='#67a8cd'),
                marker=dict(symbol='circle' if condition == 'Balanced' else 'square'),
                error_y=dict(type="constant", value=0.012, visible=True)
            ))


    fig = update_figure_format(fig)

    fig.update_layout(
    title=None,
    xaxis_title=r'$\LARGE \tau$',
    yaxis_title=r'$ \LARGE \hat{\rho}$',
    showlegend=True,
    legend=dict(
        title_text='π<sub>0</sub>',
        orientation="h",  # Horizontal legend
        yanchor="bottom",
        y=-0.3,  # Position legend below the plot
        xanchor="center",
        x=0.5,
        title_font=dict(family="Times New Roman", size=20),
        font=dict(family="Times New Roman", size=20)
    )
    )

    return fig

# def get_nabla_abs_diff(pi, Q1, Q2, t):
#     nabla1 = calculate_non_stationarity(pi, Q1, t)
#     nabla2 = calculate_non_stationarity(pi, Q2, t)
#     nabla_diff_abs_diff = abs(nabla1-nabla2)
#     return nabla_diff_abs_diff


# def get_ingroup_jsd(pi, Q1, Q2, t):
#     p1 = scipy.linalg.expm(Q1*t)
#     p2 = scipy.linalg.expm(Q2*t)
#     pi_1 = np.dot(pi,p1)
#     pi_2 = np.dot(pi,p2)
#     jsd_value = jsd(pi_1, pi_2)
#     return jsd_value

# correlation_factors_low_nabla = {time: correlation_results_low_nabla[time][0] for time in correlation_results_low_nabla.keys()}
# correlation_factors_high_nabla = {time: correlation_results_high_nabla[time][0] for time in correlation_results_high_nabla.keys()}
# correlation_factors_low_ingroup_jsd = {time: correlation_results_low_ingroup_jsd[time][0] for time in correlation_results_low_ingroup_jsd.keys()}
# correlation_factors_high_ingroup_jsd = {time: correlation_results_high_ingroup_jsd[time][0] for time in correlation_results_high_ingroup_jsd.keys()}

