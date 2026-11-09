"""Section polar and two-lag approximation of Theodorsen aerodynamics."""
import numpy as np
from scipy.special import hankel2

JONES_A = np.array([0.165, 0.335])
JONES_BETA = np.array([0.0455, 0.3])

def make_polar(cl_alpha, cl_max):
    """Super-basic NACA 0012 polar: linear through zero, tanh-saturating at +-cl_max.
    Swap for np.interp on an XFOIL polar if you want a real stall curve."""
    return lambda a: cl_max * np.tanh(cl_alpha * a / cl_max)

def exact_C(k):
    h1, h0 = hankel2(1, k), hankel2(0, k)
    return h1 / (h1 + 1j * h0)


def jones_C(k):
    ik = 1j * np.asarray(k)
    return 1 - sum(A * ik / (ik + B) for A, B in zip(JONES_A, JONES_BETA))
