"""Track oscillatory eigenvalue branches across airspeed."""
import numpy as np
from scipy.optimize import linear_sum_assignment

def sweep(model, Vs):
    """Track the two oscillatory branches across V. Returns arrays (nV, 2) of f [Hz], zeta."""
    f = np.full((len(Vs), 2), np.nan); z = np.full_like(f, np.nan)
    prev = None
    for i, V in enumerate(Vs):
        osc, _ = model.split(model.eigs(V))
        if osc.size < 2:
            continue
        osc = osc[np.argsort(osc.imag)][:2]
        if prev is not None:
            cost = np.abs(prev[:, None] - osc[None, :])
            r, c = linear_sum_assignment(cost)
            osc = osc[c]
        prev = osc
        f[i] = osc.imag / (2 * np.pi)
        z[i] = -osc.real / np.abs(osc)
    return f, z
