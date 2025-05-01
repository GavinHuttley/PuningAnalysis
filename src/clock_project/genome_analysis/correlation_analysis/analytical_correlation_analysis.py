import numpy as np
from clock_project.simulation.magnitude_quantification import calculate_non_stationarity, calculate_ENS

from cogent3.maths.measure import jsd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import statsmodels.api as sm

import scipy
import scipy.linalg

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

def get_nabla_abs_diff(pi, Q1, Q2, t):
    nabla1 = calculate_non_stationarity(pi, Q1, t)
    nabla2 = calculate_non_stationarity(pi, Q2, t)
    nabla_diff_abs_diff = abs(nabla1-nabla2)
    return nabla_diff_abs_diff


def get_ingroup_jsd(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    jsd_value = jsd(pi_1, pi_2)
    return jsd_value


def get_jad_difference(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    jsd_1 = jsd(pi_1, pi)
    jsd_2 = jsd(pi_2, pi)
    jsd_diff = abs(jsd_1 - jsd_2)
    return jsd_diff


import statsmodels.api as sm
from scipy.stats import spearmanr

def plot_time_grouped_scatter_2x2(df, x_col, y_col):
    custom_colorscale = {0.5: '#67a8cd',  # Color for time 0.5
    1.0: '#6fba4f',  # Color for time 1.0
    1.5: '#f0a3f8',  # Color for time 1.5
    2.0: '#c66f70'    # Color for time 2.0
}
    times = sorted(df['Time'].unique())

    fig = make_subplots(rows=2, cols=2,
                        horizontal_spacing=0.05, vertical_spacing=0.2)

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
            marker=dict(color = custom_colorscale[time],
                size=1  # Adjust size as needed
            )),
            row=position[0], col=position[1])
        
        corr, p = spearmanr(x_data_no_outliers, y_data_no_outliers)

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
            text=f'\u03C1 = {corr:.2f}',
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

    fig.update_layout(
        template='plotly_white',
        showlegend=True,  # Ensure the legend is shown
        legend_title_text='Time',  # Add a title to the legend to indicate what it represents
        width=1200,               # Increased width for better visibility
        height=800,               # Increased height for better visibility
        margin=dict(l=50, r=50, t=80, b=60),  # Adjusted margins
        coloraxis=dict(
        colorscale='tealrose',   # Warm colorscale: Yellow-Orange-Red
        colorbar=dict(
            title='Time',      # Title for the color bar
            title_font=dict(size=10, family='Arial', color='black'),
            tickfont=dict(size=10),
            len=0.75,          # Length of the color bar
            thickness=20,      # Thickness of the color bar
            x=1,               # Position of the color bar along the x-axis
            y=0.5,             # Position of the color bar along the y-axis
            bgcolor='white',   # Background color of the color bar
        )
    ),

    # Legend customization
    legend=dict(
        title='Time',                 # Legend title
        title_font=dict(size=16, family='Arial', color='black'),
        font=dict(size=14, family='Arial', color='black'),
        bgcolor='rgba(255,255,255,0)',  # Transparent background
        bordercolor='Black',
        borderwidth=1,
        orientation='h',              # Horizontal legend
        x=0.5,                         # Centered horizontally
        y=1.05,                        # Positioned above the plot
        xanchor='center',
        yanchor='bottom'
    ))

    fig.update_traces(
    marker=dict(
        size=3,               # Increased marker size for better visibility
        opacity=1,          # Slight transparency
        line=None  # Subtle border
    )
)

    return fig









































