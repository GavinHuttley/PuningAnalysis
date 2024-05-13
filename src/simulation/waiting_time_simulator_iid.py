import numpy as np
import random
from scipy.linalg import expm
from cogent3 import make_seq
from cogent3.util.dict_array import DictArrayTemplate

def join_number_to_base_cogent3(seq):
    number_to_base = {0: 'T', 1: 'C', 2: 'A', 3: 'G'}
    ances_seq_join_alpha = ''.join(number_to_base[number] for number in seq)

    return make_seq(''.join(ances_seq_join_alpha), moltype='dna') 

def convert_sequence_to_numeric(ancestor_seq):
    base_to_number = {'T': 0, 'C': 1, 'A': 2, 'G': 3}
    numeric_seq = [base_to_number[base] for base in str(ancestor_seq)]
    return np.array(numeric_seq)


def generate_ancestor(n, pi=None):
    """
    Generate an ancestor DNA sequence of length n with customizable probability distribution.

    Input: 
        n (integer): The desired length of the ancestor DNA sequence.
        pi (list, optional): List of probabilities for each nucleotide. If None, a random distribution is used.

    Output: 
        list: The randomly generated DNA sequence of length n in list format.
    """
    nucleotides = [0, 1, 2, 3]  # T, C, A, G

    if pi is None:
        # If no probability distribution is provided, use a random distribution, equal probabilities of 0.25
        return list(random.choices(nucleotides, k=n))
    else:
        # If a custom probability distribution is provided, use it
        if len(pi) != 4:
            raise ValueError("Probability distribution must contain exactly 4 values.")
        
        # Check if the probabilities sum to 1
        if abs(sum(pi) - 1) > 1e-10:
            raise ValueError("Probabilities must sum to 1.")
        
        return list(random.choices(nucleotides, weights=pi, k=n))
    
def generate_rate_matrix():
    """
    Generate a single 4 by 4 rate matrix.

    Input: 
        None

    Output: 
        rate_matrix (array): A single rate matrix.

    """
    matrix = np.zeros((4, 4))
    for i in range(4):
        row_sum = 0 # sum of non-diagonal elements of current row
        for j in range(4):
            if i != j: # fill up non-diagonal elements of current row
                element = np.random.uniform(0.01, 1.0)  # Non-diagonal elements are between 0.01 and 1.0
                row_sum += element
                matrix[i, j] = element
        matrix[i,i] = -row_sum # Ensure every row adds up to 0 

    return matrix

def generate_rate_matrix_cogent3() -> DictArrayTemplate: 
    """
    Generate a single 4 by 4 rate matrix.

    Output: 
        DictArray: A single rate matrix.
    """
    matrix = np.zeros((4, 4))
    for i in range(4):
        row_sum = 0  # sum of non-diagonal elements of current row
        for j in range(4):
            if i != j:  # fill up non-diagonal elements of current row
                element = np.random.uniform(0.01, 1.0)  # Non-diagonal elements are between 0.01 and 1.0
                row_sum += element
                matrix[i, j] = element
        matrix[i, i] = -row_sum  # Ensure every row adds up to 0 

    template = DictArrayTemplate(['T', 'C', 'G', 'A'], ['T', 'C', 'G', 'A'])
    return template.wrap(matrix)

    

def transition_matrix(Q, t):
    """
    Calculate the transition matrix P from rate matrix Q for a given time t.

    Input:
        Q: Rate matrix (2D array)
        t: Time passed (float)

    Output:
        P: Transition matrix (2D array)
    """
    Qt = np.dot(Q, t)
    P = expm(Qt)
    return P

def initialize_waiting_times_vectorized(DNA_seq, Q_dict):
    n_bases = len(DNA_seq)
    waiting_times = np.full((n_bases, 4), float('inf'))  # Initialize all waiting times to infinity

    # Iterate through each possible nucleotide transition
    for next_base in range(4):
        rates = np.array([Q_dict[curr_base, next_base] if curr_base != next_base else float('inf') 
                          for curr_base in DNA_seq])

        # Calculate waiting times where rates are not infinity
        valid_rates = rates != float('inf')
        waiting_times[valid_rates, next_base] = np.random.exponential(scale=1.0 / rates[valid_rates])

    # Find the minimum time and its position
    min_position = np.unravel_index(np.argmin(waiting_times), waiting_times.shape)
    min_time = waiting_times[min_position]

    return min_position, min_time
    
def initialize_waiting_times(DNA_seq, Q_dict):
    """
    Input:
        DNA_seq (list): The DNA sequence as a list.
        Q_dict (dict): A dictionary of rate matrices for context-dependent DNA substitution. 
            The key is a tuple of left and right neighbours as strings from permutations with [0, 1, 2, 3] (=[T, C, A, G]). 
            The value is the corresponding rate matrix.
        markov_order (int): The Markov order.

    Output: 
        waiting_times (array): A 2D numpy array with dimensions (length of DNA sequence) x 4. Every row contains the 4 
            waiting times of 4 possible bases. The current base waiting time is set to inf.
        min_position (tuple): A tuple containing the index of the base in the DNA sequence with the smallest
            substitution time and the next base that it will be substituted with. 
        min_time (float): The minimum time needed for next substitution to take place in the DNA sequence.
    """
    waiting_times = np.full((len(DNA_seq),4), float('inf'))
    min_time = float('inf')
    min_position = None

    for seq_index, curr_base in enumerate(DNA_seq):
        curr_base = DNA_seq[seq_index]

        for next_base in range(4):
            if next_base != curr_base:
                rate = 1/(Q_dict[curr_base, next_base])
                time = np.random.exponential(rate)
                waiting_times[seq_index, next_base] = time
                if time < min_time:
                    min_time = time
                    min_position = seq_index, next_base

    return min_position, min_time


def update_waiting_times(DNA_seq, Q_dict, min_position, min_time):
    """
    Input:
        DNA_seq (list): The DNA sequence as a list.
        Q_dict (dict): A dictionary of rate matrices for context-dependent DNA substitution. 
            The key is a tuple of left and right neighbours as strings from permutations with [0, 1, 2, 3] (=[T, C, A, G]). 
            The value is the corresponding rate matrix.
        waiting_times (array): A 2D numpy array with dimensions (length of DNA sequence) x 4. Every row contains the 4 
            waiting times of 4 possible bases. The current base waiting time is set to inf.
        min_position (tuple): A tuple containing the index of the base in the DNA sequence with the smallest
            substitution time and the next base that it will be substituted with. 
        min_time (float): The minimum time needed for next substitution to take place in the DNA sequence.
        markov_order (int): The Markov order.

    Output:
        waiting_times (array): A 2D numpy array with dimensions (length of DNA sequence) x 4. Every row contains the 4 
            waiting times of 4 possible bases. The current base waiting time is set to inf.
        min_position (tuple): A tuple containing the index of the base in the DNA sequence with the smallest
            substitution time and the next base that it will be substituted with. 
        min_time (float): The minimum time needed for next substitution to take place in the DNA sequence.
    """

    seq_index, new_base = min_position
    
    DNA_seq[seq_index] = new_base # Substitue with new base in DNA sequence 

    # Regenerate all waiting times
    min_position, min_time = initialize_waiting_times_vectorized(DNA_seq, Q_dict)

    return min_position, min_time

def simulate_seq(max_time, rate_matrix, ancestor_seq = None):
    """
    Input:
        ancestor_seq (list): The initial non-substituted generated DNA sequence of length n.
        max_time (float): The maximum time allowed for substitution(s) to take place in the simulation.
        rate_matrices (dict): A dictionary of rate matrices for context-dependent DNA substitution. 
            The key is a tuple of left and right neighbours as strings from permutations with [0, 1, 2, 3] (=[T, C, A, G]). 
            The value is the corresponding rate matrix.
        markov_order (int): The Markov order.

    Output:
        tuple: A tuple containing the final DNA sequence at the end of the simulation as well as an array
        containing the history of substituted DNA sequences throughout the simulation.
        
    """
    history = [ancestor_seq.copy(),]

    time_passed = 0

    DNA_seq = ancestor_seq.copy()
    min_position, min_time = initialize_waiting_times_vectorized(DNA_seq, rate_matrix)

    while time_passed + min_time <= max_time:
        time_passed += min_time
        seq_index, new_base = min_position
        min_position, min_time = update_waiting_times(DNA_seq, rate_matrix, min_position, min_time)
        DNA_seq[seq_index] = new_base
        history.append(DNA_seq.copy())

    return history[:-1]
    
def average_substitution(Q, t, repeats, length, pi):
    """This one from Puning, it calculates the average number of substitution per site during the defined time period of evolution for each repeat. 
    """
    
    n = length
    ns_per_site_list = []
    ns_total_list = []
    for i in range(repeats):
        ancestor_sequence = generate_ancestor(length, pi)
        history = simulate_seq(t, Q, ancestor_sequence)
        ns_total = len(history)-1
        ns_per_site = ns_total/n
        ns_per_site_list.append(ns_per_site)
        ns_total_list.append(ns_total)
    
    average_ns_total = np.average(ns_total_list)
    average_ns_per_site = average_ns_total/n

    return ns_per_site_list, average_ns_per_site
    

import plotly.graph_objects as go
from plotly.subplots import make_subplots

def get_histograms2(ns_dict, theoretical_ns_list):
    lengths = list(ns_dict.keys())
    all_times = {time for length in ns_dict for time in ns_dict[length].keys()}
    times = sorted(all_times)  # Sorting to maintain a consistent order

    rows = len(lengths)
    cols = len(times)
    
    # Determine global minimum and maximum x values for axis range consistency
    x_min = min(min(ns_dict[length][time]['ns_per_site_list']) for length in ns_dict for time in ns_dict[length] if time in ns_dict[length])
    x_max = max(max(ns_dict[length][time]['ns_per_site_list']) for length in ns_dict for time in ns_dict[length] if time in ns_dict[length])
    
    # Create subplots
    fig = make_subplots(rows=rows, cols=cols, subplot_titles=[f'Time = {t}, Length = {l}' for l in lengths for t in times])

    # Populate subplots
    for row, length in enumerate(lengths, start=1):
        for col, time in enumerate(times, start=1):
            data = ns_dict[length][time]['ns_per_site_list'] if time in ns_dict[length] else []
            theoretical_ns = theoretical_ns_list[time] if time in theoretical_ns_list else None
            average_ns = ns_dict[length][time]['avg_ns_per_site'] if time in ns_dict[length] else None

            # Add histogram to subplot
            fig.add_trace(
                go.Histogram(
                    x=data,
                    xbins=dict(  # Control the bar widths here
                        start=x_min,
                        end=x_max,
                        size=(x_max - x_min) / 20  # Adjust size for consistent bar width
                    ),
                    marker=dict(line=dict(width=1)),
                    name=f'Length {length}, Time {time}'
                ),
                row=row,
                col=col
            )
            
            # Add vertical lines for average and theoretical values
            fig.add_vline(x=average_ns, line_width=2, line_dash="dash", line_color="red", row=row, col=col)
            fig.add_vline(x=theoretical_ns, line_width=2, line_dash="dash", line_color="green", row=row, col=col)

    # Update layout for all subplots
    fig.update_layout(
        height=300 * rows,
        width=300 * cols,
        showlegend=False,
        bargap=0.05  # Adjust space between bars
    )
    
    # Set consistent x-axis range across all subplots
    for r in range(1, rows+1):
        for c in range(1, cols+1):
            fig.update_xaxes(range=[-0.09, 0.06], row=r, col=c)

    return fig
