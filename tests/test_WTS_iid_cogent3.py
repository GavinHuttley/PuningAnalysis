import sys
sys.path.append('/Users/gulugulu/repos/PuningAnalysis/src')
import pytest
import numpy as np

from simulation.waiting_time_simulator import SeqSimulate, generate_rate_matrix_cogent3, transform_Q_to_array
from cogent3 import make_seq
from cogent3.util.dict_array import DictArrayTemplate

Q = np.array([[-1.707,  0.537,  0.306,  0.864],
    [ 0.249, -0.889 ,  0.116,  0.525 ],
    [ 0.038,  0.182, -0.555,  0.335],
    [ 0.203,  0.580,  0.234, -1.017]])
template = DictArrayTemplate(['T', 'C', 'G', 'A'], ['T', 'C', 'G', 'A'])
Q_matrix = template.wrap(Q)

@pytest.fixture

def setup_simulation():
    # Set parameters
    max_time = 2
    length = 100
    num_simulations = 5
    seed = 123
    Q = Q_matrix

    simulator = SeqSimulate(max_time=max_time, Q=Q, length=length, num_simulations=num_simulations, seed=seed)
    return simulator


def test_number_of_simulations(setup_simulation):
    results = setup_simulation.main()
    assert len(results) == setup_simulation.num_simulations, "The number of simulations does not match the expected."

def test_sequence_length(setup_simulation):
    results = setup_simulation.main()
    expected_length = setup_simulation.n
    all_lengths_correct = all(len(str(seq)) == expected_length for result in results for seq in result.seqs)
    assert all_lengths_correct, "Not all sequences match the expected length."

def test_deterministic_output(setup_simulation):
    # Setting the same parameter as the fixture to test if the outcome stay the same for the same seed
    max_time = 2
    length = 100
    num_simulations = 5
    seed = 123
    Q = Q_matrix
    simulator_test = SeqSimulate(max_time=max_time, Q=Q, length=length, num_simulations=num_simulations, seed=seed)
    actual_result = simulator_test.main()
    results = setup_simulation.main()
    assert actual_result == results, "The output sequence does not match the expected deterministic sequence."

def test_single_substitution_difference(setup_simulation):
    results = setup_simulation.main()
    for simulation in results:
        for i in range(1, len(simulation)):
            prev_seq = str(simulation.seqs[i-1])
            curr_seq = str(simulation.seqs[i])
            # Check for exactly one character difference
            differences = sum(1 for a, b in zip(prev_seq, curr_seq) if a != b)
            assert differences == 1, f"Found {differences} differences between sequences, expected exactly 1."

# def test_single_substitution_difference(setup_simulation):
#     results = setup_simulation.main()
#     for simulation in results:
#         for i in range(1, len(simulation)):
#             prev_seq = str(simulation[i-1])
#             curr_seq = str(simulation[i])
#             # Check for exactly one character difference
#             differences = sum(1 for a, b in zip(prev_seq, curr_seq) if a != b)
#             assert differences == 1, f"Found {differences} differences between sequences, expected exactly 1."



