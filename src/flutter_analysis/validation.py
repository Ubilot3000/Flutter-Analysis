"""Physical and numerical checks for the reduced-order model."""
from dataclasses import replace
import numpy as np
from .parameters import Params
from .model import HStab
from .aerodynamics import exact_C, jones_C

def selftest(p: Params):
    ok = True
    def report(name, cond, msg):
        nonlocal ok; ok &= bool(cond)
        print(f"[{'PASS' if cond else 'FAIL'}] {name}: {msg}")

    # 1. Jones lag approximation vs exact Theodorsen function
    ks = np.linspace(0.02, 3.0, 200)
    err = np.max(np.abs(exact_C(ks) - jones_C(ks)))
    report("C(k) Jones vs exact Theodorsen", err < 0.02, f"max |error| = {err:.4f}")

    # 2. static torsional divergence: eigenvalue crossing vs closed form
    pd = replace(p, x_ea=0.40, x_cg=0.40, f_twist=8.0, f_roll=12.0)   # soft twist, EA aft of 1/4c
    m = HStab(pd)
    Vd_num, Vd_ana = m.divergence_speed(), m.divergence_speed_analytic()
    report("divergence speed", abs(Vd_num - Vd_ana) / Vd_ana < 0.01,
           f"eigenvalue sweep {Vd_num:.2f} m/s vs closed form {Vd_ana:.2f} m/s")

    # 3. cg position: aft cg must destabilise, forward cg must help
    res = {}
    for lbl, xcg in (("cg aft of EA +10%c", p.x_ea + 0.10), ("cg on EA", p.x_ea),
                     ("cg fwd of EA -10%c", p.x_ea - 0.10)):
        res[lbl] = HStab(replace(p, x_cg=xcg)).flutter_speed(150.0)
    v = list(res.values())
    report("cg location trend", v[0] < v[1] <= v[2], ", ".join(f"{k}: {x:.1f}" for k, x in res.items()) + " m/s")

    # 4. Zero-speed limit: compare to the coupled structure including added mass.
    # Input damping ratios refer to uncoupled dry coordinates, so they need not
    # equal modal damping once mass coupling and apparent mass are included.
    mm = HStab(p)
    osc, _ = mm.split(mm.eigs(1e-4))
    zeta = np.sort(-osc.real / np.abs(osc))
    still_air = np.block([[np.zeros((2, 2)), np.eye(2)],
                          [-np.linalg.solve(mm.Mtot, mm.K), -np.linalg.solve(mm.Mtot, mm.C)]])
    ref_ev = np.linalg.eigvals(still_air)
    ref_ev = ref_ev[ref_ev.imag > 0]
    expected = np.sort(-ref_ev.real / np.abs(ref_ev))
    report("zero-speed damping limit", np.allclose(zeta, expected, atol=1e-5),
           f"modal zeta = {np.round(zeta, 5)} vs coupled still-air {np.round(expected, 5)}")

    # 5. strip-count convergence
    vf = [HStab(replace(p, n_strips=n)).flutter_speed(150.0) for n in (6, 12, 48)]
    if np.all(np.isfinite(vf)):
        report("strip convergence", abs(vf[1] - vf[2]) / vf[2] < 0.01,
               f"Vf(6,12,48 strips) = {np.round(vf, 2)}")
    elif np.all(np.isinf(vf)):
        # Absence of a crossing in the search range does not define an infinite
        # physical boundary. Compare finite modal growth rates instead.
        speed = min(15.0, .5 * mm.divergence_speed_analytic())
        growth = [HStab(replace(p, n_strips=n)).flutter_margin(speed) for n in (6, 12, 48)]
        report("strip convergence", abs(growth[1]-growth[2]) / max(abs(growth[2]), 1e-6) < .01,
               f"No crossing through 150 m/s; growth rates at {speed:g} m/s: {np.round(growth, 4)}")
    else:
        report("strip convergence", False, f"inconsistent crossings: {vf}")

    # 6. time sim reproduces the linear eigen-solution (start on the least-damped oscillatory eigenvector)
    Vt = min(30.0, 0.8 * vf[1], 0.8 * mm.divergence_speed_analytic())
    w, vec = np.linalg.eig(mm.jacobian(Vt))
    idx = np.where(w.imag > 3.0)[0]
    k = idx[np.argmax(w.real[idx])]
    s0 = np.real(vec[:, k] * 1e-3 / abs(vec[0, k]))
    t, y = mm.simulate(Vt, T=0.5, s0=s0)
    ref = np.real(np.outer(vec[:, k] * 1e-3 / abs(vec[0, k]), np.exp(w[k] * t)))
    rel = np.max(np.abs(y[0] - ref[0])) / np.max(np.abs(ref[0]))
    report("time sim vs eigen-solution", rel < 0.02,
           f"max rel. deviation of tip motion = {rel:.4f} (lambda = {w[k].real:.1f} + {w[k].imag/(2*np.pi):.1f} Hz i)")
    return ok
