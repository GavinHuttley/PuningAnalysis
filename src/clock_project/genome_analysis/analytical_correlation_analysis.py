import numpy as np
from clock_project.simulation.magnitude_quantification import calculate_non_stationarity, calculate_ENS

from cogent3.maths.measure import jsd
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import statsmodels.api as sm

import scipy
import scipy.linalg
from scipy.special import kl_div
from scipy.stats import wasserstein_distance

import scipy

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
        if (lower_bound1 <= val1 <= upper_bound1) and (lower_bound2 <= val2 <= upper_bound2):
            filtered_data1.append(val1)
            filtered_data2.append(val2)

    return filtered_data1, filtered_data2

def remove_outliers_iqr_database(data, columns):
    def compute_iqr_bounds(column_data):
        Q1 = np.percentile(column_data, 25)
        Q3 = np.percentile(column_data, 75)
        IQR = Q3 - Q1
        lower_bound = Q1 - 3 * IQR
        upper_bound = Q3 + 3 * IQR
        return lower_bound, upper_bound

    # Make a copy of the data to avoid modifying the original DataFrame
    data_filtered = data.copy()

    for column in columns:
        lower_bound, upper_bound = compute_iqr_bounds(data_filtered[column])
        # Filter out rows where the column value is outside the IQR bounds
        data_filtered = data_filtered[(data_filtered[column] >= lower_bound) & (data_filtered[column] <= upper_bound)]

    return data_filtered

def get_ens_diff_log_ratio(pi, Q1, Q2, t):
    ens1 = calculate_ENS(pi, Q1, t)
    ens2 = calculate_ENS(pi, Q2, t)
    ens_diff = np.log(ens1/ens2)
    return abs(ens_diff)

def get_ens_abs_diff(pi, Q1, Q2, t):
    ens1 = calculate_ENS(pi, Q1, t)
    ens2 = calculate_ENS(pi, Q2, t)
    ens_absdiff = abs(ens1-ens2)
    return ens_absdiff

def get_ens_Hellinger(pi, Q1, Q2, t):
    ens1 = calculate_ENS(pi, Q1, t)
    ens2 = calculate_ENS(pi, Q2, t)
    ens_Hellinger = np.sqrt(2*(np.sqrt(ens1)- np.sqrt(ens2))**2)
    return ens_Hellinger

# def get_nabla_diff_log_ratio(pi, Q1, Q2, t, ens_o_1, ens_o_2):
#     nabla1 = calculate_non_stationarity(pi, Q1, t)/ens_o_1
#     nabla2 = calculate_non_stationarity(pi, Q2, t)/ens_o_2
#     nabla_diff_log_ratio = np.log(nabla1/nabla2)
#     return abs(nabla_diff_log_ratio)

def get_nabla_diff_log_ratio(pi, Q1, Q2, t):
    nabla1 = calculate_non_stationarity(pi, Q1, t)
    nabla2 = calculate_non_stationarity(pi, Q2, t)
    nabla_diff_log_ratio = np.log(nabla1/nabla2)
    return abs(nabla_diff_log_ratio)

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

def get_ingroup_wst(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    wst_value = wasserstein_distance(pi_1, pi_2)
    return wst_value

def get_jsd_difference(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    jsd_1 = jsd(pi_1, pi)
    jsd_2 = jsd(pi_2, pi)
    jsd_diff = abs(jsd_1 - jsd_2)
    return jsd_diff

def get_wasserstein_difference(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    wst_1 = wasserstein_distance(pi, pi_1)
    wst_2 = wasserstein_distance(pi, pi_2)
    wst_diff = abs(wst_1 - wst_2)
    return wst_diff


def get_KL_divergence(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    KL1 = sum(kl_div(pi_1, pi))
    KL2 = sum(kl_div(pi_2, pi))
    return KL1, KL2

def get_KL_divergence_diff(pi, Q1, Q2, t):
    p1 = scipy.linalg.expm(Q1*t)
    p2 = scipy.linalg.expm(Q2*t)
    pi_1 = np.dot(pi,p1)
    pi_2 = np.dot(pi,p2)
    KL1 = sum(kl_div(pi_1, pi))
    KL2 = sum(kl_div(pi_2, pi))
    return np.sqrt(2*(np.sqrt(KL1)-np.sqrt(KL2))**2)

#Plotting function

def bar_data(position3d, size=(1,1,1)):
    # Generate the vertices of a parallelepipedic bar at a specified position and size
    bar = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],
                    [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]], dtype=float)
    bar *= np.array(size)
    bar += np.array(position3d)
    return bar

def triangulate_bar_faces(positions, sizes):
    # Triangulate the faces of multiple bars to generate vertices and indices for Mesh3d
    all_bars = [bar_data(pos, size) for pos, size in zip(positions, sizes)]
    vertices, ixr = np.unique(np.vstack(all_bars), return_inverse=True, axis=0)
    
    I, J, K = [], [], []
    for k in range(len(all_bars)):
        indices = ixr[k * 8:(k + 1) * 8]
        I.extend(indices[[0, 2, 0, 5, 0, 7, 5, 2, 3, 6, 7, 5]])
        J.extend(indices[[1, 3, 4, 1, 3, 4, 1, 6, 7, 2, 4, 6]])
        K.extend(indices[[2, 0, 5, 0, 7, 0, 2, 5, 6, 3, 5, 7]])
    return vertices, I, J, K

def get_plotly_mesh3d(x, y, bins=[10, 10], bargap=0.1):
    # Generate a 3D histogram plot data
    hist, xedges, yedges = np.histogram2d(x, y, bins=bins)
    xpos, ypos = np.meshgrid(xedges[:-1] + np.diff(xedges) / 2,
                             yedges[:-1] + np.diff(yedges) / 2, indexing="ij")
    
    positions = np.column_stack([xpos.ravel(), ypos.ravel(), np.zeros(xpos.size)])
    sizes = np.column_stack([np.full(xpos.size, xedges[1] - xedges[0] - bargap),
                             np.full(ypos.size, yedges[1] - yedges[0] - bargap),
                             hist.ravel()])
    
    vertices, I, J, K = triangulate_bar_faces(positions, sizes)
    return vertices[:, 0], vertices[:, 1], vertices[:, 2], I, J, K


def process_data_for_plot(df, x_col, y_col, x_precision, y_precision):
    # Remove rows with any NaN values
    df = df.dropna()

    # Convert differences to numpy arrays
    x = df[x_col].to_numpy()
    y = df[y_col].to_numpy()

    # Round the differences and update the DataFrame
    df[x_col] = df[x_col].round(x_precision)
    df[y_col] = df[y_col].round(y_precision)

    # Group by the rounded values and count occurrences
    density_data = df.groupby([x_col, y_col]).size().reset_index(name='Density')

    # Call the prepared function for 3D mesh computation
    X, Y, Z, I, J, K = get_plotly_mesh3d(x, y, bins=[20, 20], bargap=0.05)

    return X, Y, Z, I, J, K, density_data


# def create_3d_density_plots(X1, Y1, Z1, I1, J1, K1, X2, Y2, Z2, I2, J2, K2):
#     # Create subplots: 1 row, 2 columns
#     fig = make_subplots(rows=1, cols=2, specs=[[{'type': 'mesh3d'}, {'type': 'mesh3d'}]],
#                         horizontal_spacing=0,
#                         subplot_titles=('Low Information', 'High Information'))

#     # Add first 3D mesh plot to the first subplot
#     fig.add_trace(go.Mesh3d(
#         x=X1, y=Y1, z=Z1,
#         i=I1, j=J1, k=K1,
#         intensity=Z1,  # Typically uses the Z values or another metric for color intensity
#         colorscale='Viridis',  # Reversed Viridis color scale; remove '_r' for normal progression
#         showscale=False,
#         opacity=1,  # Set opacity to make overlaps more discernible
#         coloraxis="coloraxis"
#     ), row=1, col=1)

#     # Add second 3D mesh plot to the second subplot
#     fig.add_trace(go.Mesh3d(
#         x=X2, y=Y2, z=Z2,
#         i=I2, j=J2, k=K2,
#         intensity=Z2,  # Typically uses the Z values or another metric for color intensity
#         colorscale='Viridis',  # Reversed Viridis color scale; remove '_r' for normal progression
#         coloraxis="coloraxis",
#         opacity=1  # Set opacity to make overlaps more discernible
#     ), row=1, col=2)

#     return fig

from scipy.stats import linregress


def create_2d_density_plots(density_data1, density_data2, x_col, y_col):
    # Sort the data by density so that higher density points are plotted last (on top)
    density_data_sorted1_origin = density_data1.sort_values(by='Density', ascending=True)
    density_data_sorted2_origin = density_data2.sort_values(by='Density', ascending=True)

    density_data_sorted1 = remove_outliers_iqr_database(density_data_sorted1_origin, [x_col, y_col])
    density_data_sorted2 = remove_outliers_iqr_database(density_data_sorted2_origin, [x_col, y_col])


    # Create subplots: 1 row, 2 columns
    fig = make_subplots(rows=1, cols=2, subplot_titles=(
        'Balanced',
        'Un-balanced'
    ),
    horizontal_spacing=0.05)

    # Add first scatter plot to the first subplot
    fig.add_trace(go.Scatter(
        x=density_data_sorted1[x_col],
        y=density_data_sorted1[y_col],
        mode='markers',
        marker=dict(
            size=8,  # Adjust size as needed
            color=density_data_sorted1['Density'], 
            colorscale='Viridis', 
            coloraxis="coloraxis",
            opacity=1, 
            showscale=True  
        )
    ), row=1, col=1)

    slope1, intercept1, r_value, p_value, std_err = linregress(density_data_sorted1[x_col], density_data_sorted1[y_col])
    fig.add_trace(go.Scatter(
        x=density_data_sorted1[x_col],
        y=slope1 * density_data_sorted1[x_col] + intercept1,
        mode='lines',
        marker=dict(color='red'),
        name='Fit 1'
    ), row=1, col=1)

    # Add second scatter plot to the second subplot
    fig.add_trace(go.Scatter(
        x=density_data_sorted2[x_col],
        y=density_data_sorted2[y_col],
        mode='markers',
        marker=dict(
            size=8,  # Adjust size as needed
            color=density_data_sorted2['Density'], 
            colorscale='Viridis', 
            coloraxis="coloraxis",
            opacity=1,  
            showscale=True 
        )
    ), row=1, col=2)

    slope2, intercept2, r_value, p_value, std_err = linregress(density_data_sorted2[x_col], density_data_sorted2[y_col])
    fig.add_trace(go.Scatter(
        x=density_data_sorted2[x_col],
        y=slope2 * density_data_sorted2[x_col] + intercept2,
        mode='lines',
        marker=dict(color='red'),
        name='Fit 2'
    ), row=1, col=2)

    # fig.update_yaxes(type='log', row=1, col=2)
    # fig.update_yaxes(type='log', row=1, col=1)

    return fig


import statsmodels.api as sm

def plot_time_grouped_scatter_2x2(df, x_col, y_col):
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
            name=f'Time {time}',  # This will automatically create legend entries
            marker=dict(
                size=2  # Adjust size as needed
            )),
            row=position[0], col=position[1])

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
            text=f"R² = {model.rsquared:.2f}",
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