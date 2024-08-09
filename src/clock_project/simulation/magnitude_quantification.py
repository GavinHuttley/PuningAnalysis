from cogent3.maths.matrix_exponential_integration import expected_number_subs
from numpy import dot
from scipy.linalg import expm
import numpy as np


unifor_distribution = [0.25, 0.25, 0.25, 0.25]

def calculate_ENS(p0, Q, t):
    return expected_number_subs(p0, Q, t)

def calculate_non_stationarity(p0, Q, t):
    pi_deriv = dot(p0, dot(Q, expm(Q * t)))
    conv = np.linalg.norm(pi_deriv)
    return conv

def entropy_calculation(pk):
    return -np.sum(pk*np.log2(pk))

def calculate_information(p0):
    return entropy_calculation(unifor_distribution) - entropy_calculation(p0)


def random_nucleotide_distribution():
    # Generate 4 random numbers
    random_numbers = np.random.rand(4)
    # Normalize these numbers so that their sum equals 1
    distribution = random_numbers / random_numbers.sum()
    return distribution