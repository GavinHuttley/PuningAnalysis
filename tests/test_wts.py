import pytest
import numpy as np
from numpy.random import default_rng
from clock_project.simulation.wts import SeqSimulate, generate_ancestor


@pytest.fixture
def equal_frequency():
    return [0.25, 0.25, 0.25, 0.25]

@pytest.fixture
def sequence_length():
    return int(10000)

@pytest.fixture
def rng():
    return default_rng(seed=42)

@pytest.fixture
def Q_matrix():
    return np.array([
        [-1.707,  0.537,  0.306,  0.864],
        [ 0.249, -0.889,  0.116,  0.525],
        [ 0.038,  0.182, -0.555,  0.335],
        [ 0.203,  0.580,  0.234, -1.017]
    ])

@pytest.fixture
def max_time():
    return int(1)

@pytest.fixture
def setup_simulation(Q_matrix):
    pi = [0.15, 0.25, 0.35, 0.25]
    return SeqSimulate(Q=Q_matrix, length=100, num_repeat=5, seed=123, pi=pi)


def test_number_of_simulations(setup_simulation, max_time):
    results = setup_simulation.main(max_time)
    expected_num_simulations = setup_simulation.num_repeat
    assert len(results) == expected_num_simulations, "The number of simulations does not match the expected."

def test_sequence_length(setup_simulation, max_time):
    results = setup_simulation.main(max_time)
    expected_length = setup_simulation.length
    all_lengths_correct = all(len(seq) == expected_length for result in results for seq in result)
    assert all_lengths_correct, "Not all sequences match the expected length."

def test_deterministic_output(setup_simulation, Q_matrix, max_time):
    # Create a separate instance with identical parameters to compare deterministic output
    other_simulator = SeqSimulate(Q=Q_matrix, length=100, num_repeat=5, seed=123, pi = [0.15, 0.25, 0.35, 0.25])
    other_results = other_simulator.main(max_time)
    results = setup_simulation.main(max_time)
    assert all(results) == all(other_results), "The output sequence does not match the expected deterministic sequence."

def test_single_substitution_difference(setup_simulation, max_time):
    results = setup_simulation.main(max_time)
    for simulation in results:
        for i in range(1, len(simulation)):
            prev_seq = str(simulation[i-1])
            curr_seq = str(simulation[i])
            # Check for exactly one character difference
            differences = sum(1 for a, b in zip(prev_seq, curr_seq) if a != b)
            assert differences == 1, f"Found {differences} differences between sequences, expected exactly 1."


@pytest.mark.parametrize("expected_count", [2500])
def test_generate_ancestor_nucleotide_distribution(equal_frequency, sequence_length, expected_count):
    seq = generate_ancestor(sequence_length, equal_frequency)
    counts = np.bincount(seq, minlength=4)
    print("Generated sequence counts:", counts)  # Debug output
    assert all(abs(count - expected_count) <= 0.05 * expected_count for count in counts)

