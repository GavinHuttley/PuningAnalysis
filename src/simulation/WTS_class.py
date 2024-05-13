import numpy as np
from cogent3 import make_seq
from numpy.random import SeedSequence, default_rng
from cogent3 import make_unaligned_seqs
from cogent3 import Sequence
from waiting_time_simulator_iid import generate_rate_matrix_cogent3


def transform_Q_to_array(Q_dict):
    """
    Transforms a dictionary-based Q matrix to a numpy array.
    Assumes Q_dict is a dictionary of dictionaries with keys 'T', 'C', 'G', 'A'.

    Args:
        Q_dict (dict): The rate matrix as a dictionary of dictionaries.

    Returns:
        numpy.ndarray: A 4x4 numpy array representing the rate matrix.
    """
    base_to_index = {'T': 0, 'C': 1, 'G': 2, 'A': 3}
    matrix = np.zeros((4, 4))
    for from_base, transitions in Q_dict.items():
        from_index = base_to_index[from_base]
        for to_base, rate in transitions.items():
            to_index = base_to_index[to_base]
            matrix[from_index, to_index] = rate
    return matrix

def convert_sequence_to_numeric(seq: Sequence):
    base_to_number = {'T': 0, 'C': 1, 'G': 2, 'A': 3}
    numeric_seq = [base_to_number[base] for base in str(seq)]
    return np.array(numeric_seq)

def join_number_to_base_cogent3(seq):
    number_to_base = {0: 'T', 1: 'C', 2: 'G', 3: 'A'}
    ances_seq_join_alpha = ''.join(number_to_base[number] for number in seq)

    return make_seq(''.join(ances_seq_join_alpha), moltype='dna') 

def convert_sequence_history_to_cogent3_sequence_colletion(history):
    seq_names = []
    for i in range(len(history)):
        if i == len(history)-1:
            seq_names.append('ancestor_seq')
        elif i == len(history)-1:
            seq_names.append('final_seq')
        else:
            seq_names.append(f'internediate_seq_{i}')
    
    seq_base = (join_number_to_base_cogent3(seq) for seq in history)
    data = zip(seq_names, seq_base)
    seqs = make_unaligned_seqs(data, 'dna')
    return seqs

class SeqSimulate:
    def __init__(self, Q: np.array, length: int, num_repeat: int, seed: int, pi: list = None):
        self.Q = Q
        self.length = length
        self.num_repeat = num_repeat
        self.pi = pi if pi is not None else [0.25, 0.25, 0.25, 0.25] 
        self.seed = seed
        self.ss = SeedSequence(self.seed)
        self.child_seeds = self.ss.spawn(self.num_repeat)
        self.rngs = [default_rng(s) for s in self.child_seeds]
    
    def main(self, max_time):
        results = []
        for i in range(self.num_repeat):
            rng = self.rngs[i]
            ancestor_sequence = self.generate_ancestor(rng)
            simulation_result = self.simulate_seq(ancestor_sequence, max_time, rng)

            results.append(simulation_result)
        return results


    def generate_ancestor(self, rng):
        nucleotides = [0, 1, 2, 3]  # T, C, A, G
        return list(rng.choice(nucleotides, size=self.length, p=self.pi))

    def initialize_waiting_times_vectorized(self, DNA_seq, rng):
        n_bases = len(DNA_seq)
        waiting_times = np.full((n_bases, 4), float('inf'))  # Initialize all waiting times to infinity

        # Iterate through each possible nucleotide transition
        for next_base in range(4):
            rates = np.array([self.Q[curr_base, next_base] if curr_base != next_base else float('inf') 
                            for curr_base in DNA_seq])

            # Calculate waiting times where rates are not infinity
            valid_rates = rates != float('inf')
            waiting_times[valid_rates, next_base] = rng.exponential(scale=1.0 / rates[valid_rates])

        # Find the minimum time and its position
        min_position = np.unravel_index(np.argmin(waiting_times), waiting_times.shape)
        min_time = waiting_times[min_position]

        return min_position, min_time

    # def initialize_waiting_times(self, DNA_seq):
    #     waiting_times = np.full((len(DNA_seq), 4), float('inf'))
    #     min_time = float('inf')
    #     min_position = None

    #     for seq_index in range(len(DNA_seq)):
    #         curr_base = DNA_seq[seq_index]
    #         for next_base in range(4):
    #             if next_base != curr_base:
    #                 rate = 1 / self.Q[curr_base, next_base]
    #                 time = np.random.exponential(rate)
    #                 waiting_times[seq_index, next_base] = time
    #                 if time < min_time:
    #                     min_time = time
    #                     min_position = (seq_index, next_base)

    #     return waiting_times, min_position, min_time

    def update_waiting_times(self, DNA_seq, min_position, min_time, rng):
        seq_index, new_base = min_position
        DNA_seq[seq_index] = new_base  # Substitue with new base in DNA sequence
        min_position, min_time = self.initialize_waiting_times_vectorized(DNA_seq, rng)
        return min_position, min_time

    def simulate_seq(self, ancestor_seq, max_time, rng):
        history = [ancestor_seq.copy()]
        time_passed = 0
        DNA_seq = ancestor_seq.copy()
        min_position, min_time = self.initialize_waiting_times_vectorized(DNA_seq, rng)

        while time_passed <= max_time:
            seq_index, new_base = min_position
            substitution_time = min_time
            min_position, min_time = self.update_waiting_times(DNA_seq, min_position, min_time, rng)
            DNA_seq[seq_index] = new_base
            history.append(DNA_seq.copy())
            time_passed += substitution_time

        return history[:-1]

    def average_substitution(self, max_time):
        ns_per_site_list = []
        ns_total_list = []
        results = self.main(max_time)
        for simulation in results:
            num_subt = len(simulation) - 1
            ns_per_site = num_subt/self.length
            ns_per_site_list.append(ns_per_site)
            ns_total_list.append(num_subt)
        
        average_ns_total = np.average(ns_total_list)
        average_ns_per_site = average_ns_total/self.length
        
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


#how to use this class with cogent3 data

# seq = make_seq('ACGTTCGACGA', moltype='dna')
# pi1 = [0.1, 0.2, 0.3, 0.4]
# Q1_cogent3 = generate_rate_matrix_cogent3()

# Q1 = Q1_cogent3.array

# simulator = SeqSimulate(Q=Q1, length=100, num_repeat=100, seed = 123, pi = pi1)

# result = simulator.main(max_time= 1)

# cogent3_seqs = convert_sequence_history_to_cogent3_sequence_colletion(result)

