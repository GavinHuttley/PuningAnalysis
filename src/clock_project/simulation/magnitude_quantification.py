from cogent3.maths.matrix_exponential_integration import expected_number_subs
from numpy import dot
from scipy.linalg import expm
import numpy as np




def calculate_ENS(p0, Q, t):
    return expected_number_subs(p0, Q, t)

def calculate_non_stationarity(p0, Q, t):
    pi_deriv = dot(p0, dot(Q, expm(Q * t)))
    conv = np.linalg.norm(pi_deriv)
    return conv