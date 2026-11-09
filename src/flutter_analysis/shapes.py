"""Tip-normalized structural mode shapes."""
import numpy as np

def _cantilever(eta):
    bl, sig = 1.8751040687, 0.7340955137
    f = lambda x: np.cosh(x) - np.cos(x) - sig * (np.sinh(x) - np.sin(x))
    return f(bl * eta) / f(bl)


def _shape(kind, eta):
    table = {
        "rigid_roll": lambda: eta,
        "cantilever": lambda: _cantilever(eta),
        "rigid_pitch": lambda: np.ones_like(eta),
        "linear": lambda: eta,
        "sine": lambda: np.sin(0.5 * np.pi * eta),
    }
    return table[kind]()
