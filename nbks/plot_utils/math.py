import numpy as np
from cogent3.maths.matrix_exponential_integration import expected_number_subs

from clock_project.maths.evolutionary_rate import (
    calculate_non_stationary_mu_prime,
    calculate_non_stationary_rate,
    calculate_stationary_distribution,
    calculate_stationary_rate,
)
from plot_utils.util import update_figure_format

# t_range = np.linspace(0, 10, 999)
# Q1 = array([[-1.4, 0.1, 0.4, 0.9],
#            [4.0, -6.9, 0.9, 2.0],
#            [6.3, 2.0, -11.3, 3.0],
#            [0.7,0.1, 0.2, -1]])

# Q2 = np.array([[-8.35, 1.98, 4.58, 1.79],
#            [4.69, -6.76, 1.08, 0.99],
#            [0.05, 0.58, -3.24, 2.61],
#            [0.05,1.32, 3.70, -5.07]])

# Q3= np.array([
#     [-4.5,  2.0,  1.0,  1.5],
#     [ 2.0, -3.5,  0.5,  1.0],
#     [ 1.0,  0.5, -2.5,  1.0],
#     [ 0.0,  1.0,  2.0, -3.0]
# ])


# Calculate the range of evolution rate and its derivation given a time intervel and initial nucleotide frequency
def mu_mu_prime_range(Q, pi, t_range):
    """
    Calculate the range of evolution rate and its derivation
    given a time intervel and initial nucleotide frequency

    Parameters:
    Q (numpy.ndarray): The substitution rate matrix.
    pi_0 (numpy.ndarray): The initial nucleotide frequency distribution.
    t_range (numpy.linspace): The time interval for calculating the evolutionary rate and its derivative

    Returns:
    lists: The calculated value of μ(t) and μ'(t) over the time inteval as lists.
    """
    mu_range = []
    mu_prime_range = []
    for t in t_range:
        mu = calculate_non_stationary_rate(Q, pi, t)
        mu_prime = calculate_non_stationary_mu_prime(Q, pi, t)
        mu_range.append(mu)
        mu_prime_range.append(mu_prime)

    return mu_range, mu_prime_range


def generate_ENS(Q, t_range, pi=None):
    """
    Generates the ENS over a range of time points using two different Q matrices before and after a specified time point t1.

    Parameters:
    - pi: A numpy array of shape (1, 4) representing the vector pi.
    - Q1, Q2: Two numpy arrays of shape (4, 4) representing the original and new rate matrices, with different stationary distribution.
    - t_range: numpy.linspace defining the start and end of the time range.
    - t1: The time point at which to switch from using Q1 to Q2.

    Returns:
    - A list of ENS values for each time point in the range.
    """
    ens_values = []

    for t in t_range:
        ens = expected_number_subs(
            pi, Q, t
        )  # Update the accumulated ENS at t1 to continue from this point using Q2
        ens_values.append(ens)

    return ens_values


def get_evolutionary_rate_plot(Q, t_range, pi=None):
    fig = go.Figure()
    # Ploting the evolutionary rate mu over time with rate matrix switch in the middle
    mu_value = []
    # Generate mu values
    if pi is None:
        for t in t_range:
            mu = calculate_stationary_rate(Q)
            mu_value.append(mu)
    else:
        for t in t_range:
            mu = calculate_non_stationary_rate(Q, pi, t)
            mu_value.append(mu)

    # Add Q3 line
    fig.add_trace(
        go.Scatter(
            x=np.linspace(0, 10, 999),
            y=mu_value,
            mode="lines",
            line=dict(color="#58B8D1", width=4),
        )
    )

    # Update layout for better presentation
    fig.update_layout(
        title=None,
        title_font=dict(size=20),  # Increase title size, make it bold
        xaxis_title="<b><i>t</i></b>",
        yaxis_title=r"$\mu$",
        legend=dict(
            title=None,
            x=1.02,  # Adjust position of legend
            y=1,
        ),
    )
    fig = update_figure_format(fig)
    return fig


def get_ens_plot(Q, t_range, pi=None):

    if pi is None:
        pi = calculate_stationary_distribution(Q)
        ens_values = []
        for t in t_range:
            ens = -np.sum(pi * np.diag(Q)) * t
            ens_values.append(ens)
    else:
        ens_values = generate_ENS(Q, t_range, pi)
    # Create the plot
    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=t_range, y=ens_values, mode="lines", line=dict(color="#cd6d2e", width=4)
        )
    )

    # Update layout for better presentation
    fig.update_layout(
        title=None,
        title_font=dict(size=20),  # Increase title size, make it bold
        xaxis_title="<b><i>t</i></b>",
        yaxis_title=r"$ENS$",
    )
    fig = update_figure_format(fig)
    return fig


import plotly.graph_objects as go
from plotly.subplots import make_subplots


def get_evolutionary_rate_change_plot(Q1, Q2, Q3, pi, t_range):

    t = np.linspace(t_range[0], t_range[1], 100)

    Q = [Q1, Q2, Q3]
    titles = [
        "<b>Monotonic-Decrease</b>",
        "<b>Monotonic-Increase</b>",
        "<b>Non-Monotonic</b>",
    ]

    # create 1x3 subplots with subplot titles
    fig = make_subplots(rows=1, cols=3, subplot_titles=titles)

    # add each Q to its own subplot
    for i, Q in enumerate(Q, start=1):
        mu_range, _ = mu_mu_prime_range(Q, pi, t_range)  # using your existing helper
        # ensure mu_range length matches x; if not, try broadcasting/trimming
        y = np.asarray(mu_range)
        if y.size != t.size:
            # fallback: if mu_range itself provides x-like axis, use it directly; else interpolate/trim
            try:
                # if mu_range is returned with same length as internal time grid inside mu_mu_prime_range,
                # assume that's the x-axis and plot against np.arange
                x_plot = np.arange(y.size)
            except Exception:
                x_plot = t
        else:
            x_plot = t

        fig.add_trace(
            go.Scatter(
                x=x_plot,
                y=y,
                mode="lines",
                line=dict(color="#1f77b4", width=4),
                showlegend=False,
            ),
            row=1,
            col=i,
        )

        # set axis labels for this subplot
        fig.update_xaxes(title_text="<b><i>t</i></b>", row=1, col=i)
        fig.update_yaxes(title_text=r"$\mu$", row=1, col=1)

    # style subplot title font size (these are annotations created by plotly)
    for a in fig.layout.annotations:
        a.font = dict(size=16)

    # keep your post-processing function if you have one
    fig = update_figure_format(fig)

    # optionally a main title or layout tweaks
    fig.update_layout(height=380, width=1200)
    # tickvals = [0, 0.05, 0.10, 0.15, 0.20]
    # fig.update_yaxes(tickmode="array", tickvals=tickvals, showgrid=True)  # all subplots
    # fig.update_yaxes(title_text=r'$\mu$', row=1, col=1)

    return fig


def get_ens_difference_change_two_process(Q1, Q2, t_range, pi=None):
    if pi is None:
        pi = calculate_stationary_distribution(Q1)
    # Ploting the ENS over time with rate matrix switch in the middle
    t_range = np.linspace(0, 1, 999)

    ens_values = generate_ENS(Q1, pi, t_range)
    ens_values2 = generate_ENS(Q2, pi, t_range)
    # Create the plot
    fig = go.Figure()
    # Add Q5 line
    fig.add_trace(
        go.Scatter(
            x=t_range,
            y=ens_values2,
            mode="lines",
            name="Species 1",
            line=dict(color="#62BE9F", width=4),
        )
    )

    # Add Q3 line
    fig.add_trace(
        go.Scatter(
            x=t_range,
            y=ens_values,
            mode="lines",
            name="Species 2",
            line=dict(color="#AC8AB3", width=4),
        )
    )

    # Update layout for better presentation
    fig.update_layout(
        title=None,
        xaxis_title="<b><i>t</i></b>",
        yaxis_title=r"$ENS$",
        legend=dict(
            title=None,
            font=dict(size=16),  # Adjust font size of legend
            x=0.5,  # Center the legend horizontally
            y=-0.2,  # Position the legend below the plot
            xanchor="center",
            yanchor="top",  # Anchor the legend to the top to ensure it's below the plot
            orientation="h",  # Horizontal orientation
        ),
    )
    fig = update_figure_format(fig)
    return fig


def get_clock_violation_change_two_process(Q1, Q2, t_range, pi=None):
    if pi is None:
        pi = calculate_stationary_distribution(Q1)

    mu_range1, _ = mu_mu_prime_range(Q1 * 0.2, pi, t_range)
    mu_range2, _ = mu_mu_prime_range(Q2 * 0.2, pi, t_range)
    fig = go.Figure()
    fig.add_trace(
        go.Scatter(
            x=t_range,
            y=mu_range1,
            mode="lines",
            name="Species 2",
            line=dict(color="#AC8AB3", width=3),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=t_range,
            y=mu_range2,
            mode="lines",
            name="Species 1",
            line=dict(color="#62BE9F", width=3),
        )
    )

    fig.add_trace(
        go.Scatter(
            x=t_range,
            y=np.abs(np.array(mu_range1) - np.array(mu_range2)),
            mode="lines",
            name="Difference (Clock violation)",
            line=dict(color="#ff0b00", width=3, dash="dash"),
        )
    )

    # Update layout for better presentation
    fig.update_layout(
        title=None,
        xaxis_title="<b><i>t</i></b>",
        yaxis_title=r"$\mu(t)$",
        legend=dict(
            title=None,
            font=dict(size=16),  # Adjust font size of legend
            x=0.5,  # Center the legend horizontally
            y=-0.2,  # Position the legend below the plot
            xanchor="center",
            yanchor="top",  # Anchor the legend to the top to ensure it's below the plot
            orientation="h",  # Horizontal orientation
        ),
    )

    fig = update_figure_format(fig)
    return fig
