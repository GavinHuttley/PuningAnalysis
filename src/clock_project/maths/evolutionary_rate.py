import numpy as np
from scipy.linalg import expm
from numpy import dot
import math

#stationary process 
def calculate_stationary_distribution(Q):
    """
    Calculate the stationary distribution pi for a given substitution rate matrix Q.

    Parameters:
    Q (numpy.ndarray): The substitution rate matrix.

    Returns:
    numpy.ndarray: The stationary distribution pi.
    """
    # Add an additional equation to account for the sum of pi elements being 1
    A = np.vstack([Q.T, np.ones(Q.shape[0])])

    b = np.zeros(Q.shape[0] + 1)
    b[-1] = 1

    # Solve for pi
    pi = np.linalg.lstsq(A, b, rcond=None)[0]

    return pi


def calculate_stationary_rate(Q):
    """
    Calculate the stationary evolution rate mu_stationary for a given substitution rate matrix Q.

    Parameters:
    Q (numpy.ndarray): The substitution rate matrix.

    Returns:
    float: The stationary evolution rate mu_stationary.
    """
    # Get the stationary disitrbution of Q
    pi_stationary = calculate_stationary_distribution(Q)

    # Calculate stationary evolution rate using the formula mu = - sum_i(pi*Qii)
    mu_stationary = - math.fsum(pi_stationary*np.diagonal(Q))
    return mu_stationary

def matrix_calibration(Q):
    Q_c = Q/ (- math.fsum(calculate_stationary_distribution(Q)*np.diagonal(Q)))
    return Q_c


#stationary process ENS accumulation function
def generate_ENS(Q1, Q2, t_range, t1):
    """
    Generates the ENS over a range of time points using two different Q matrices before and after a specified time point t1.
    
    Parameters:
    - pi: A numpy array of shape (1, 4) representing the vector pi.
    - Q1, Q2: Two numpy arrays of shape (4, 4) representing the original and new rate matrices.
    - t_range: numpy.linspace defining the start and end of the time range.
    - t1: The time point at which to switch from using Q1 to Q2.
    
    Returns:
    - A list of ENS values for each time point in the range.
    """
    pi = calculate_stationary_distribution(Q1)

    ens_values = []
    ens_accumulated = 0  # To keep track of the accumulated ENS value
    
    for t in t_range:
        if t <= t1:
            ens = -np.sum(pi * np.diag(Q1)) * t + ens_accumulated
            ens_values.append(ens)
        else:
            ens_2 = -np.sum(pi * np.diag(Q2)) * (t-t1) + ens # Update the accumulated ENS at t1 to continue from this point using Q2
            ens_values.append(ens_2)
    
    return ens_values

#non-stationary process 
#evolution rate at each time point in non-stationary process

def calculate_non_statioanry_mu(Q, pi_0, t):
    """
    Calculate the value of mu prime (μ(t)) for a given substitution rate matrix Q,
    initial nucleotide frequency pi_0, and time t.

    Parameters:
    Q (numpy.ndarray): The substitution rate matrix.
    pi_0 (numpy.ndarray): The initial nucleotide frequency distribution.
    t (float): The time at which to calculate μ(t).

    Returns:
    float: The calculated value of μ(t).
    """
    # Calculate f(t) = pi_0 * exp(Qt)
    f_t = dot(pi_0, expm(Q * t))
    
    # Calculate mu'(t) as the sum of the element-wise product of f(t) and the diagonal of Q
    mu = - dot(f_t,np.diagonal(Q))
    
    return mu

#derivative of evolution rate
def calculate_non_stationary_mu_prime(Q, pi_0, t):
    """
    Correctly calculate the value of mu prime (μ'(t)) based on the provided formula:
    μ'(t) = - pi_0 * Q * exp(Qt) * diag(Q)

    Parameters:
    Q (numpy.ndarray): The substitution rate matrix.
    pi_0 (numpy.ndarray): The initial nucleotide frequency distribution.
    t (float): The time at which to calculate μ'(t).

    Returns:
    numpy.ndarray: The calculated value of μ'(t) as a vector.
    """
    
    # Calculate μ'(t) using the provided formula
    mu_prime_t = -pi_0.dot(Q).dot(expm(Q * t)).dot(np.diagonal(Q))
    
    return mu_prime_t

#stationary process ENS accumulation function
def generate_ENS(Q1, Q2, t_range, t1):
    """
    Generates the ENS over a range of time points using two different Q matrices before and after a specified time point t1.
    
    Parameters:
    - pi: A numpy array of shape (1, 4) representing the vector pi.
    - Q1, Q2: Two numpy arrays of shape (4, 4) representing the original and new rate matrices.
    - t_range: numpy.linspace defining the start and end of the time range.
    - t1: The time point at which to switch from using Q1 to Q2.
    
    Returns:
    - A list of ENS values for each time point in the range.
    """
    pi = calculate_stationary_distribution(Q1)

    ens_values = []
    ens_accumulated = 0  # To keep track of the accumulated ENS value
    
    for t in t_range:
        if t <= t1:
            ens = -np.sum(pi * np.diag(Q1)) * t + ens_accumulated
            ens_values.append(ens)
        else:
            ens_2 = -np.sum(pi * np.diag(Q2)) * (t-t1) + ens # Update the accumulated ENS at t1 to continue from this point using Q2
            ens_values.append(ens_2)
    
    return ens_values