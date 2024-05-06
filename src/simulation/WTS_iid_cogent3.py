from typing import Optional
import numpy as np
from scipy.linalg import expm
import cogent3
import random
from cogent3 import make_seq
from cogent3.util.dict_array import DictArrayTemplate
from cogent3.core.sequence import Sequence

from cogent3.app.composable import define_app
from cogent3.app.typing import SeqType, SeqsCollectionType
from cogent3 import make_unaligned_seqs
from numpy.random import default_rng



def generate_ancestor_cogent3(n: int, rng, pi=None) -> Sequence:
    """
    Generate an ancestor DNA sequence of length n with customizable probability distribution.

    Input: 
        n (int): The desired length of the ancestor DNA sequence.
        pi (list, optional): Array of probabilities for each nucleotide. 
            If None, a random distribution is used.

    Output: 
        cogent3.core.alignment.Alignment: The randomly generated DNA sequence of length n.
    """
    nucleotides = ['T', 'C', 'A', 'G']

    if pi is None:
        pi = [0.25, 0.25, 0.25, 0.25]  # Default equal probabilities if none provided
    elif len(pi) != 4:
        raise ValueError("Probability distribution must contain exactly 4 values.")
    elif not np.isclose(sum(pi), 1):
        raise ValueError("Probabilities must sum to 1.")

    # Using np.random.Generator.choice instead of random.choices
    seq = rng.choice(nucleotides, size=n, p=pi)
    ancestor_sequence = make_seq(''.join(seq), moltype='dna')

    return ancestor_sequence


def generate_rate_matrix_cogent3(rng):
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


def transition_matrix_cogent3(Q, t):
    """
    Calculate the transition matrix P from rate matrix Q for a given time t.

    Input:
        Q: Rate matrix (2D array)
        t: Time passed (float)

    Output:
        P: Transition matrix (2D array)
    """
    Qt = Q.array*t
    P = expm(Qt)
    template = DictArrayTemplate(['T', 'C', 'G', 'A'], ['T', 'C', 'G', 'A'])

    return template.wrap(P)


def initialize_waiting_times_iid_cogent3(DNA_seq: list, Q: DictArrayTemplate) -> tuple:
    """
    Initialize waiting times for DNA sequence substitution.

    Input:
        DNA_seq (list): The DNA sequence as a list.
        Q_dict (DictArray): A DictArray of rate matrices for context-dependent DNA substitution.

    Output:
        tuple: A tuple containing the waiting times, the position of the base with the smallest substitution time,
               and the minimum time needed for the next substitution.
    """
    waiting_times = np.full((len(DNA_seq), 4), np.inf)
    min_time = np.inf
    min_position = None
    nucleotides = ['T', 'C', 'G', 'A']

    for seq_index in range(len(DNA_seq)):
        curr_base = str(DNA_seq[seq_index])
        for next_base in nucleotides:
            if next_base != curr_base:
                next_base_index = nucleotides.index(next_base)
                rate = 1 / Q[curr_base][next_base]
                time = np.random.exponential(rate)
                waiting_times[seq_index, next_base_index] = time
                if time < min_time:
                    min_time = time
                    min_position = (seq_index, next_base)

    return waiting_times, min_position, min_time

def update_waiting_times_iid_cogent3(DNA_seq: Sequence, Q: DictArrayTemplate, waiting_times: np.ndarray, min_position: tuple, min_time: float) -> tuple:
    """
    Update waiting times for DNA sequence substitution.

    Input:
        DNA_seq (cogent3.core.sequence.Sequence): The DNA sequence.
        Q (DictArray): A DictArrayTemplate representing the rate matrix for DNA substitution.
        waiting_times (np.ndarray): An array containing waiting times for DNA sequence substitution.
        min_position (tuple): A tuple containing the position of the base with the smallest substitution time.
        min_time (float): The minimum time needed for the next substitution.

    Output:
        tuple: A tuple containing the updated waiting times, the position of the base with the new smallest substitution time,
               and the new minimum time needed for the next substitution.
    """

    seq_index = min_position[0] # obtain the index of the current base that is about to be substituted in the sequence
    new_base = min_position[1] # obtain the new base that is going to substitute the current base
    
    seq_list = list(DNA_seq)
    seq_list[seq_index] = new_base
    DNA_seq_new = make_seq(str(''.join(seq_list)))
    
    waiting_times, min_position, min_time = initialize_waiting_times_iid_cogent3(DNA_seq_new, Q)

    return waiting_times, min_position, min_time

def simulate_seq_iid_cogent3(ancestor_seq: Sequence, max_time : float, Q: DictArrayTemplate):
    """
    Input:
        ancestor_seq (cogent3.Sequence): The initial non-substituted generated DNA sequence of length n.
        max_time (float): The maximum time allowed for substitution(s) to take place in the simulation.
        rate_matrices (DictArray): A substitution rate matrix. 

    Output:
        tuple: A tuple containing the final DNA sequence at the end of the simulation as well as an array
        containing the history of substituted DNA sequences throughout the simulation.
        
    """
    history = [ancestor_seq,]

    time_passed = 0

    DNA_seq = ancestor_seq.copy()
    waiting_times, min_position, min_time = initialize_waiting_times_iid_cogent3(DNA_seq, Q)

    while time_passed <= max_time:
        seq_index = min_position[0]
        new_base = min_position[1]
        substitution_time = min_time 
        waiting_times, min_position, min_time = update_waiting_times_iid_cogent3(DNA_seq, Q, waiting_times, min_position, min_time)
        seq_list = list(DNA_seq)
        seq_list[seq_index] = new_base
        DNA_seq = make_seq(str(''.join(seq_list)))
        history.append(DNA_seq.copy())
        time_passed += substitution_time

    if len(history) <= 2: # meaning time for any substitution to occur exceeded the max_time, no substitution took place
        aln_history = make_unaligned_seqs(history[:-1])
        return aln_history
    else:
        aln_history = make_unaligned_seqs(history[:-1])
        final_seq = history[-2]
        return (final_seq, aln_history)




###define the cogent3 app

from numpy.random import SeedSequence, default_rng


@define_app
class SeqSimulate:
    """
    Input:
        ancestor_seq (cogent3.Sequence): The initial non-substituted generated DNA sequence of length n.
        max_time (float): The maximum time allowed for substitution(s) to take place in the simulation.
        Q (DictArrayTemplate): A substitution rate matrix.

    Output:
        tuple: A tuple containing the final DNA sequence at the end of the simulation as well as an array
        containing the history of substituted DNA sequences throughout the simulation.
        
    """

    def __init__(self, max_time: float, Q: DictArrayTemplate, length: int, num_simulations: int, seed: int = None, pi: list = None):
        self.max_time = max_time
        self.Q = Q
        self.n = length
        self.num_simulations = num_simulations
        self.seed = seed
        self.pi = pi if pi is not None else [0.25, 0.25, 0.25, 0.25]
        self.ss = SeedSequence(self.seed)
        self.child_seeds = self.ss.spawn(self.num_simulations)
        self.rngs = [default_rng(s) for s in self.child_seeds]

    def main(self, ancestor_seq: SeqType = None) -> SeqsCollectionType:
        results = []
        for rng in self.rngs:
            simulation_result = self.simulate_sequence(rng, ancestor_seq)
            results.append(simulation_result)
        return results
    
    def generate_ancestor(self, rng) -> Sequence:
        """
        Generate an ancestor DNA sequence of length n with customizable probability distribution.

        Input: 
            n (int): The desired length of the ancestor DNA sequence.
            pi (list, optional): Array of probabilities for each nucleotide. 

        Output: 
            cogent3.core.alignment.Alignment: The randomly generated DNA sequence of length n.
        """
        nucleotides = ['T', 'C', 'A', 'G']  # Correctly map indices to bases directly

        if len(self.pi) != 4:
            raise ValueError("Probability distribution must contain exactly 4 values.")
        else:
            # If a custom probability distribution is provided, use it
            if len(self.pi) != 4:
                raise ValueError("Probability distribution must contain exactly 4 values.")
            
            # Check if the probabilities sum to 1
            if abs(sum(self.pi) - 1) > 1e-10:
                raise ValueError("Probabilities must sum to 1.")
            
        seq = rng.choices(nucleotides, weights=self.pi, k=self.n)
        ancestor_sequence = make_seq(''.join(seq), moltype='dna')

        return ancestor_sequence
        
    def simulate_sequence(self, rng, ancestor_seq: SeqType = None) -> SeqsCollectionType:
        if ancestor_seq == None:
            ancestor_sequence = generate_ancestor(rng)
        else:
            ancestor_sequence = ancestor_seq
        names = ['ancestor',]
        time_passed = 0
        history = [ancestor_sequence,]
        DNA_seq = ancestor_sequence.copy()
        waiting_times, min_position, min_time = self.initialize_waiting_times_iid(DNA_seq, rng)
        i = 1

        while time_passed <= self.max_time:
            seq_index = min_position[0]
            new_base = min_position[1]
            substitution_time = min_time
            waiting_times, min_position, min_time = self.update_waiting_times_iid(DNA_seq, waiting_times, min_position, min_time, rng)
            seq_list = list(DNA_seq)
            seq_list[seq_index] = new_base
            DNA_seq = make_seq(''.join(seq_list), moltype='dna')
            history.append(DNA_seq.copy())
            names.append(f'seq{i}')
            time_passed += substitution_time
            i += 1

        if len(history) <= 2: # meaning time for any substitution to occur exceeded the max_time, no substitution took place
            seq_dict = {name:seq for (name, seq) in  zip(names, history)}
            aln_history = make_unaligned_seqs(data = seq_dict, moltype='dna')
            
            return aln_history
        else:
            seq_dict = {name:seq for (name, seq) in  zip(names, history)}
            aln_history = make_unaligned_seqs(data = seq_dict, moltype='dna')
            return aln_history


    
    def initialize_waiting_times_iid(self, DNA_seq: SeqsCollectionType, rng) -> tuple:
        """
        Initialize waiting times for DNA sequence substitution.

        Input:
            DNA_seq (Sequence): The DNA sequence.

        Output:
            tuple: A tuple containing the waiting times, the position of the base with the smallest substitution time,
                   and the minimum time needed for the next substitution.
        """
        waiting_times = np.full((len(DNA_seq), 4), np.inf)
        min_time = np.inf
        min_position = None
        nucleotides = ['T', 'C', 'A', 'G']

        for seq_index in range(len(DNA_seq)):
            curr_base = str(DNA_seq[seq_index])
            for next_base in nucleotides:
                if next_base != curr_base:
                    next_base_index = nucleotides.index(next_base)
                    rate = 1 / self.Q[curr_base][next_base]
                    time = rng.exponential(rate)
                    waiting_times[seq_index, next_base_index] = time
                    if time < min_time:
                        min_time = time
                        min_position = (seq_index, next_base)

        return waiting_times, min_position, min_time

    
    
    
    
    def update_waiting_times_iid(self, DNA_seq: SeqsCollectionType, waiting_times: np.ndarray, min_position: tuple, min_time: float, rng) -> tuple:
        """
        Update waiting times for DNA sequence substitution.

        Input:
            DNA_seq (Sequence): The DNA sequence.
            Q (DictArrayTemplate): A rate matrix for DNA substitution.
            waiting_times (np.ndarray): An array containing waiting times for DNA sequence substitution.
            min_position (tuple): A tuple containing the position of the base with the smallest substitution time.
            min_time (float): The minimum time needed for the next substitution.

        Output:
            tuple: A tuple containing the updated waiting times, the position of the base with the new smallest substitution time,
                   and the new minimum time needed for the next substitution.
        """

        seq_index = min_position[0]
        new_base = min_position[1]

        # Update the sequence
        seq_list = list(DNA_seq)
        seq_list[seq_index] = new_base
        DNA_seq = make_seq(''.join(seq_list), moltype='dna')

        waiting_times, min_position, min_time = self.initialize_waiting_times_iid(DNA_seq, rng)

        # Reinitialize the waiting times
        return waiting_times, min_position, min_time
