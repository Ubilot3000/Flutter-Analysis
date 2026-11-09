"""Two-coordinate horizontal-tail strip model and time integration."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from .parameters import Params
from .aerodynamics import JONES_A, JONES_BETA, make_polar
from .shapes import _shape

class HStab:
    def __init__(self, p: Params, polar=None):
        self.p = p
        self.polar = polar or make_polar(p.cl_alpha, p.cl_max)
        N = p.n_strips
        edges = np.linspace(0.0, p.semispan, N + 1)
        self.y = 0.5 * (edges[:-1] + edges[1:])
        self.dy = np.diff(edges)
        eta = self.y / p.semispan
        c = p.c_root + (p.c_tip - p.c_root) * eta
        self.b = 0.5 * c                                   # semichord
        self.a = np.full(N, 2.0 * p.x_ea - 1.0)            # EA from midchord in semichords
        self.e = (p.x_cg - p.x_ea) * c                     # cg aft of EA [m]
        m = p.mass_half * c / np.sum(c * self.dy)          # mass/span, proportional to chord
        I_ea = m * (p.k_gyr * c) ** 2 + m * self.e ** 2
        self.phi = _shape(p.roll_shape, eta)
        self.psi = _shape(p.twist_shape, eta)
        ph, ps, dy = self.phi, self.psi, self.dy

        # structural generalised mass  (T = 1/2 m wdot^2 - m e wdot thdot + 1/2 I_ea thdot^2)
        Ms = np.array([[np.sum(m * ph ** 2 * dy), -np.sum(m * self.e * ph * ps * dy)],
                       [-np.sum(m * self.e * ph * ps * dy), np.sum(I_ea * ps ** 2 * dy)]])
        # apparent (added) mass from Theodorsen non-circulatory terms
        b, a, rho = self.b, self.a, p.rho
        Ma = np.array([[np.sum(np.pi * rho * b ** 2 * ph ** 2 * dy),
                        np.sum(np.pi * rho * b ** 3 * a * ph * ps * dy)],
                       [np.sum(np.pi * rho * b ** 3 * a * ph * ps * dy),
                        np.sum(np.pi * rho * b ** 4 * (1 / 8 + a ** 2) * ps ** 2 * dy)]])
        self.Ms, self.Mtot = Ms, Ms + Ma
        wh, wt = 2 * np.pi * p.f_roll, 2 * np.pi * p.f_twist
        self.K = np.diag([wh ** 2 * Ms[0, 0], wt ** 2 * Ms[1, 1]])
        self.C = np.diag([2 * p.zeta_roll * wh * Ms[0, 0], 2 * p.zeta_twist * wt * Ms[1, 1]])
        self.n = 4 + 2 * N

    # ----- state s = [q_h, q_th, qd_h, qd_th, lag states (2 x N)] -----
    def rhs(self, t, s, V, gust=0.0):
        p, N = self.p, self.p.n_strips
        V = max(V, 1e-6)
        q, v, x = s[0:2], s[2:4], s[4:].reshape(2, N)
        ph, ps, b, a, dy, rho = self.phi, self.psi, self.b, self.a, self.dy, p.rho
        # 3/4-chord downwash per strip (up +):  V*theta - wdot + b(1/2 - a)*thetadot
        # Positive vertical gust increases effective incidence. This uses the
        # existing wake-lag approximation, not a separate Sears gust model.
        Q = V * ps * q[1] - ph * v[0] + b * (0.5 - a) * ps * v[1] + gust
        lam = np.outer(JONES_BETA, V / b)
        xdot = -lam * x + Q
        CQ = (1.0 - JONES_A.sum()) * Q + (JONES_A[:, None] * lam * x).sum(axis=0)   # C(k)*Q
        alpha_eff = np.deg2rad(p.alpha0_deg) + CQ / V
        L = rho * V ** 2 * b * self.polar(alpha_eff)                  # circulatory lift/span, at 1/4 chord
        Fc = np.array([np.sum(L * ph * dy), np.sum(L * b * (a + 0.5) * ps * dy)])
        Fnc = np.array([v[1] * np.sum(np.pi * rho * b ** 2 * V * ph * ps * dy),
                        -v[1] * np.sum(np.pi * rho * b ** 3 * V * (0.5 - a) * ps ** 2 * dy)])
        qdd = np.linalg.solve(self.Mtot, -self.K @ q - self.C @ v + Fc + Fnc)
        return np.concatenate([v, qdd, xdot.ravel()])

    def jacobian(self, V, h=1e-7):
        A = np.zeros((self.n, self.n))
        for j in range(self.n):
            d = np.zeros(self.n); d[j] = h
            A[:, j] = (self.rhs(0, d, V) - self.rhs(0, -d, V)) / (2 * h)
        return A

    def eigs(self, V):
        return np.linalg.eigvals(self.jacobian(V))

    # ----- stability metrics -----
    @staticmethod
    def split(ev, w_min=3.0):
        osc = ev[ev.imag > w_min]                       # oscillatory (one per conj pair)
        real = ev[np.abs(ev.imag) <= w_min].real
        return osc, real

    def flutter_margin(self, V):
        osc, _ = self.split(self.eigs(V))
        return osc.real.max() if osc.size else -1.0

    def divergence_margin(self, V):
        _, real = self.split(self.eigs(V))
        return real.max() if real.size else -1.0

    def find_critical(self, fn, vmax=100.0, dv=0.5):
        Vs = np.arange(1.0, vmax + dv, dv)
        prev = fn(Vs[0])
        for V0, V1 in zip(Vs[:-1], Vs[1:]):
            cur = fn(V1)
            if prev < 0 <= cur:
                return brentq(fn, V0, V1, xtol=1e-3)
            prev = cur
        return np.inf

    def flutter_speed(self, vmax=100.0):
        return self.find_critical(self.flutter_margin, vmax)

    def divergence_speed(self, vmax=200.0):
        return self.find_critical(self.divergence_margin, vmax, dv=1.0)

    def divergence_speed_analytic(self):
        """Static torsional divergence: K_th = rho V^2 cl_alpha sum b^2 (a+1/2) psi^2 dy."""
        p = self.p
        S = np.sum(self.b ** 2 * (self.a + 0.5) * self.psi ** 2 * self.dy)
        return np.sqrt(self.K[1, 1] / (p.rho * p.cl_alpha * S)) if S > 0 else np.inf

    # ----- time simulation -----
    def simulate(self, V, T=1.5, q0=(0.004, np.deg2rad(0.5)), fs=2000, s0=None):
        if not np.isfinite(V) or V <= 0:
            raise ValueError("simulation airspeed must be finite and positive")
        if not np.isfinite(T) or T <= 0 or not np.isfinite(fs) or fs <= 0:
            raise ValueError("duration and sample rate must be finite and positive")
        if s0 is None:
            s0 = np.zeros(self.n); s0[:2] = q0
        t_eval = np.arange(0, T, 1 / fs)
        sol = solve_ivp(self.rhs, (0, T), s0, args=(V,), t_eval=t_eval,
                        method="LSODA", rtol=1e-8, atol=1e-10)
        if not sol.success or not np.all(np.isfinite(sol.y)):
            raise RuntimeError(f"Time integration failed: {sol.message}")
        return sol.t, sol.y

    def summary(self):
        p = self.p
        Mt, Ms = self.Mtot, self.Ms
        print(f"  generalised mass  M_h={Ms[0,0]:.4g} kg   I_th={Ms[1,1]:.4g} kg m^2   "
              f"coupling={Ms[0,1]:.3g}")
        w = np.sqrt(np.abs(np.linalg.eigvals(np.linalg.solve(Mt, self.K))))
        print(f"  natural freqs incl. added mass & coupling: {np.sort(w) / (2*np.pi)} Hz "
              f"(inputs {p.f_roll}, {p.f_twist})")
        print(f"  EA at {p.x_ea:.2f}c, cg at {p.x_cg:.2f}c  ->  cg aft of EA by {(p.x_cg-p.x_ea)*100:.0f}% chord")
