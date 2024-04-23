from collections import Counter
import itertools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import random
from scipy.linalg import expm

def generate_ancestor(n, pi=None):
    """
    Generate an ancestor DNA sequence of length n with customizable probability distribution.

    Input: 
        n (integer): The desired length of the ancestor DNA sequence.
        pi (list, optional): List of probabilities for each nucleotide. If None, a random distribution is used.

    Output: 
        list: The randomly generated DNA sequence of length n in list format.
    """
    nucleotides = ['0', '1', '2', '3']  # T, C, A, G

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

    rate_matrix = matrix.tolist()
    return rate_matrix

def generate_rate_matrices(markov_order):
    """
    Input:
        markov_order (int): The Markov order.

    Output: 
        rate_matrices (dict): A dictionary of rate matrices for context-dependent DNA substitution. 
            The key is a tuple of left and right neighbours as strings from permutations with [0, 1, 2, 3] (=[T, C, A, G]). 
            The value is the corresponding rate matrix.

    Example output: 
    generate_rate_matrices(1):
        {('0','0') :[[-0.8, 0.1, 0.4, 0.3],
                    [0.3, -0.7, 0.2, 0.2],
                    [0.1, 0.2, -0.6, 0.3],
                    [0.1, 0.2, 0.4, -0.7]],
        ('0','1') : [[...]...], ...}
    """
    rate_matrices = {}

    if markov_order == 0:
        rate_matrix = generate_rate_matrix()
        rate_matrices['0'] = rate_matrix # For standardizing so that output is always a dictionary of matrices.
        return rate_matrices
    
    else:
        nucleotides = ['0', '1', '2', '3']
        one_side_neighbours = [''.join(p) for p in itertools.product(nucleotides, repeat=markov_order)]
        permutations = list(itertools.product(one_side_neighbours, repeat=2))

        for perm in permutations:
            rate_matrix = generate_rate_matrix()
            rate_matrices[perm] = rate_matrix

        return rate_matrices
    

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

def get_context(seq_index, DNA_seq, markov_order, anchor_base = '0'):
    """
    Given the index of a base in the given DNA sequence list and order, obtain the context i.e. left and right neighbours in a tuple.

    Input: 
        seq_index (int): The index of a base in the given DNA sequence list.
        DNA_seq (list): The DNA sequence as a list.
        markov_order (int): The Markov order.
        anchor_base (string): The anchoring base beyond the DNA sequence to address incomplete nucleotide contexts at the left and right ends of the sequence.

    Output:
        context (tuple): The left and right neighbours of the base at the given seq_index.

    Example output: 
    get_context(1, ['3', '2', '3', '1', '2'], 2) = ('03', '31')
    """

    DNA_seq = ''.join(DNA_seq) # convert list DNA_seq to string 
    
    if markov_order == 0:
        return '0' # needed to access single independent matrix in rate_matrices_dict
    
    elif markov_order > 0:

        if ((seq_index < markov_order) & (len(DNA_seq) - (seq_index+1) < markov_order)): # handles incomplete left and right neighbours
            left_base = ''.join([(anchor_base * (markov_order-seq_index)), DNA_seq[:seq_index]])
            right_base =''.join([DNA_seq[seq_index+1:], (anchor_base * (markov_order - (len(DNA_seq)-seq_index-1)))])
            context = (left_base, right_base)

        elif seq_index < markov_order: # handles incomplete left neighbour
            left_base = ''.join([(anchor_base * (markov_order-seq_index)), DNA_seq[:seq_index]])
            right_base = DNA_seq[seq_index+1:seq_index+1+markov_order]
            context = (left_base, right_base)

        elif (len(DNA_seq) - (seq_index+1) < markov_order): # handles incomplete right neighbour
            left_base = DNA_seq[seq_index-markov_order:seq_index]
            right_base = ''.join([DNA_seq[seq_index+1:], (anchor_base * (markov_order - (len(DNA_seq)-seq_index-1)))])
            context = (left_base, right_base)

        else:
            left_base = DNA_seq[seq_index-markov_order:seq_index]
            right_base = DNA_seq[seq_index+1:seq_index+1+markov_order]
            context = (left_base, right_base)
        return context
    
def initialize_waiting_times(DNA_seq, Q_dict, markov_order):
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

    for seq_index in range(len(DNA_seq)):
        curr_base = int(DNA_seq[seq_index])

        curr_context = get_context(seq_index, DNA_seq, markov_order)

        for next_base in range(4):
            if next_base != curr_base:
                rate = 1/(Q_dict[curr_context][curr_base][next_base])
                time = np.random.exponential(rate)
                waiting_times[seq_index, next_base] = time
                if time < min_time:
                    min_time = time
                    min_position = (seq_index, next_base)

    return waiting_times, min_position, min_time

def get_neighbour_indices(DNA_seq, seq_index, markov_order):

    """
    Input:
        DNA_seq (list): The DNA sequence as a list.
        seq_index (int): The index of a base in the given DNA sequence list.
        markov_order (int): The Markov order.

    Output:
        neighbours (array): An array of indices of the left and right neighbours of the given base at the given
            index in the DNA sequence, including the given index.
    """
    neighbours = []
    sequence_length = len(DNA_seq)

    if ((seq_index < 0) | (seq_index >= sequence_length)):
        print("Error: Index out of bounds")
    
    else:

        for index in range(max(0, seq_index - markov_order), min(sequence_length, seq_index + markov_order + 1)):
            neighbours.append(index)

        return neighbours

def update_waiting_times(DNA_seq, Q_dict, waiting_times, min_position, min_time, markov_order):
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

    seq_index = min_position[0] # obtain the index of the current base that is about to be substituted in the sequence

    new_base = min_position[1] # obtain the new base that is going to substitute the current base
    
    DNA_seq[seq_index] = f"{new_base}" # Substitue with new base in DNA sequence 

    # Regenerate all waiting times
    waiting_times, min_position, min_time = initialize_waiting_times(DNA_seq, Q_dict, markov_order)

    return waiting_times, min_position, min_time

def simulate_seq(ancestor_seq, max_time, rate_matrices_dict, markov_order):
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
    history = [ancestor_seq,]

    time_passed = 0

    DNA_seq = ancestor_seq.copy()
    Q_dict = rate_matrices_dict
    waiting_times, min_position, min_time = initialize_waiting_times(DNA_seq, Q_dict, markov_order)

    while time_passed <= max_time:
        seq_index = min_position[0]
        new_base = min_position[1]
        substitution_time = min_time 
        waiting_times, min_position, min_time = update_waiting_times(DNA_seq, Q_dict, waiting_times, min_position, min_time, markov_order)
        DNA_seq[seq_index] = f"{new_base}" # Substitue with new base in DNA sequence
        history.append(DNA_seq.copy())
        time_passed += substitution_time

    if len(history) <= 2: # meaning time for any substitution to occur exceeded the max_time, no substitution took place
        return (ancestor_seq, history[:-1])
    else:
        return (history[-2], history[:-1])
    

#Simulation study

def generate_stats(theo_df, sim_df):
    """ 
    Input: 
        theo_df (dataframe): Dataframe containing all the theoretical propbabilities of each transition.
        sim_df (dataframe): Dataframe containing all the simulated proportions of each transition, each row representing one repeat.
        
    Output:
        None        
    """
    # Calculate mean, standard deviation, and standard error for each column
    mean_sim = sim_df.mean()
    std_sim = sim_df.std()
    se_sim = std_sim / np.sqrt(len(sim_df))

    # Compare with theoretical values
    for col in sim_df.columns:
        print(f"Mean Simulated Proportion for {col}: {mean_sim[col]}")
        print(f"Theoretical Probability for {col}: {theo_df[col].iloc[0]}")
        print(f"Standard Error for {col}: {se_sim[col]}")
        print(f"Z-score for {col}: {(mean_sim[col] - theo_df[col].iloc[0])/se_sim[col]}")
        print()

def count_change_props(starting_seq, final_seq):
    """ 
    Input:
        starting_seq (list): Starting DNA sequence.
        final_seq (list): Final DNA sequence after one run of simulation.

    Output:
        changes_proportions (dict): A dictionary of proportions of each possible transition.
    """
    changes_counter = Counter(zip(starting_seq, final_seq))
    changes_proportions = {}

    nucleotide_counts = {'0': 0, '1': 0, '2': 0, '3': 0}
    nucleotide_counts.update(Counter(starting_seq)) # Update counts based on the starting sequence

    seq_length = len(starting_seq)

    for from_nucleotide in '0123':
        for to_nucleotide in '0123':
            key = f"{from_nucleotide}_to_{to_nucleotide}"
            count = changes_counter[(from_nucleotide, to_nucleotide)]
            proportion = count / seq_length
            changes_proportions[key] = proportion

    return changes_proportions

def simulation2(repeat, n, pi, max_time, rate_matrices_dict, markov_order):
    """ 
    Simulate sequence for multiple times based on the repeat (int), get the transition probability between each of two nucleotide

    Input:
        repeat (int): Number of times to repeat the simulation.
        n (int): Length of DNA sequence.
        pi (list): List of probabilities for each nucleotide. 
        max_time (float): Maximum time allowed for DNA substitutions to take place for each repeat.
        rate_matrices_dict (dict): A dictionary of rate matrices for context-dependent DNA substitution. 
            The key is a tuple of left and right neighbours as strings from permutations with [0, 1, 2, 3] (=[T, C, A, G]). 
            The value is the corresponding rate matrix.
        markov_order (int): Markov order.

    Output:
        theo_probs_df (dataframe): Dataframe containing all the theoretical propbabilities of each transition.
        sim_props_df (dataframe): Dataframe containing all the simulated proportions of each transition, each row representing one repeat.
    """

    rate_mat_Q = rate_matrices_dict['0']
    trans_mat_P = transition_matrix(rate_mat_Q, max_time)
    pi_diag = np.diag(pi)
    joint_mat_J = np.dot(pi_diag, trans_mat_P) 
    theo_probs_dict = {}
    for i in range(4):
        for j in range(4):
            theo_probs_dict[f'{i}_to_{j}'] = joint_mat_J[i,j]

    simulated_props = {}

    for i in range(repeat):
        ancestor_seq = generate_ancestor(n, pi) # new ancestor every run
        simulation_results = simulate_seq(ancestor_seq,max_time,rate_matrices_dict,markov_order)
        final_seq = simulation_results[0]
        change_props = count_change_props(ancestor_seq, final_seq)
        simulated_props[i+1] = change_props

    theo_probs_df = pd.DataFrame([theo_probs_dict])
    sim_props_df = pd.DataFrame.from_dict(simulated_props, orient='index')
        
    return theo_probs_df, sim_props_df

