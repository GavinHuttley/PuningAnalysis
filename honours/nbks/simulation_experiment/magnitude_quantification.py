from cogent3.maths.matrix_exponential_integration import expected_number_subs
from numpy import dot, norm
from scipy.linalg import expm
from numpy.linalg import norm




def calculate_ENS(p0, Q, t):
    return expected_number_subs(p0, Q, t)

def calculate_non_stationarity(p0, Q, t):
    pi_deriv = dot(p0, dot(Q, expm(Q * t)))
    conv = norm(pi_deriv)
    return conv