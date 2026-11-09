"""Illustrative continuous gust forcing of two clamped semispans.

Vertical gust speed is an Ornstein-Uhlenbeck (first-order colored noise)
process with stationary RMS sigma and correlation time L/V. This is a useful
broadband demonstration, not a calibrated Dryden/von Karman atmosphere or
an identified propeller wake. No longitudinal airspeed noise is added.
"""
from dataclasses import asdict, dataclass
import math
import numpy as np


@dataclass(frozen=True)
class Excitation:
    kind: str = "gusts"       # gusts | periodic | combined | release
    pattern: str = "mixed"    # mixed | symmetric | antisymmetric
    gust_rms: float = 0.25    # vertical gust m/s per half, ensemble RMS
    length_scale: float = 0.5 # m, illustrative short-scale disturbance
    correlation: float = 0.75 # mixed-mode left/right gust correlation
    tone_amplitude: float = 0.10 # m/s, peak periodic vertical airflow
    tone_frequency: float = 20.0 # Hz
    seed: int = 42

    def __post_init__(self):
        if self.kind not in ("gusts", "periodic", "combined", "release"):
            raise ValueError("invalid excitation kind")
        if self.pattern not in ("mixed", "symmetric", "antisymmetric"):
            raise ValueError("invalid spanwise pattern")
        for name in ("gust_rms", "tone_amplitude"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name in ("length_scale", "tone_frequency"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f"{name} must be finite and positive")
        if not 0 <= self.correlation <= 1:
            raise ValueError("correlation must be between zero and one")
        if not isinstance(self.seed, int) or not 0 <= self.seed <= 0xffffffff:
            raise ValueError("seed must be an unsigned 32-bit integer")


class NormalRandom:
    """LCG + Box-Muller, shared bit-for-bit with the offline JS simulator."""
    def __init__(self, seed):
        self.seed, self.spare = seed, None

    def uniform(self):
        self.seed = (1664525 * self.seed + 1013904223) & 0xffffffff
        return (self.seed + 0.5) / 4294967296

    def normal(self):
        if self.spare is not None:
            value, self.spare = self.spare, None
            return value
        radius = math.sqrt(-2 * math.log(self.uniform()))
        angle = 2 * math.pi * self.uniform()
        self.spare = radius * math.sin(angle)
        return radius * math.cos(angle)


def gust_history(V, duration, dt, excitation):
    """Return t, [left,right] gust speed, starting in stationary noise."""
    e = excitation
    n = max(2, int(math.ceil(duration / dt)) + 1)
    times = np.arange(n) * dt
    gust = np.zeros((n, 2))
    rng = NormalRandom(e.seed)
    process = [rng.normal() for _ in range(3)]
    decay = math.exp(-V * dt / e.length_scale)
    innovation = math.sqrt(1 - decay ** 2)
    common, local = math.sqrt(e.correlation), math.sqrt(1 - e.correlation)
    for i, t in enumerate(times):
        if i:
            process = [decay * v + innovation * rng.normal() for v in process]
        left = right = 0.0
        if e.kind in ("gusts", "combined"):
            if e.pattern == "symmetric":
                left = right = e.gust_rms * process[0]
            elif e.pattern == "antisymmetric":
                right = e.gust_rms * process[0]
                left = -right
            else:
                left = e.gust_rms * (common * process[0] + local * process[1])
                right = e.gust_rms * (common * process[0] + local * process[2])
        if e.kind in ("periodic", "combined"):
            tone = e.tone_amplitude * math.sin(2 * math.pi * e.tone_frequency * t)
            right += tone
            left += -tone if e.pattern == "antisymmetric" else tone
        gust[i] = left, right
    return times, gust


def export_kernel(model):
    """Serialize the same half-tail equations for browser integration."""
    return {"params": asdict(model.p), "b": model.b.tolist(), "a": model.a.tolist(),
            "dy": model.dy.tolist(), "phi": model.phi.tolist(), "psi": model.psi.tolist(),
            "inverseMass": np.linalg.inv(model.Mtot).tolist(),
            "K": model.K.diagonal().tolist(), "C": model.C.diagonal().tolist()}


def simulate_tail(model, V, duration=6.0, fs=600, excitation=None, substeps=4):
    """Integrate both semispans with RK4; stop at the small-motion boundary.

    The halves share a fixed root and correlated airflow, with no flexible-boom,
    fuselage or across-root structural coupling. Each has bending + twist DOFs.
    """
    e = excitation or Excitation()
    if not np.isfinite(V) or V <= 0 or not np.isfinite(duration) or duration <= 0:
        raise ValueError("airspeed and duration must be finite and positive")
    if not np.isfinite(fs) or fs <= 0 or substeps < 1:
        raise ValueError("sample rate and substeps must be positive")
    dt = 1 / (fs * substeps)
    times, gust = gust_history(V, duration, dt, e)
    state = np.zeros((2, model.n))
    if e.kind == "release":
        state[:, :2] = (0.004, np.deg2rad(0.5))
        if e.pattern == "antisymmetric":
            state[0] *= -1
    history, outputs, gust_outputs = [], [], []
    stopped = False
    for i in range(len(times)):
        if i % substeps == 0:
            history.append(float(times[i]))
            outputs.append(state[:, :2].copy())
            gust_outputs.append(gust[i].copy())
            if np.any(np.abs(state[:, 0]) > .1 * model.p.semispan) or np.any(np.abs(state[:, 1]) > np.deg2rad(10)):
                stopped = True
                break
        if i + 1 == len(times):
            break
        for side in (0, 1):
            s, g0, g1 = state[side], gust[i, side], gust[i+1, side]
            gm = .5 * (g0 + g1)
            k1 = model.rhs(times[i], s, V, g0)
            k2 = model.rhs(times[i] + dt/2, s + dt/2*k1, V, gm)
            k3 = model.rhs(times[i] + dt/2, s + dt/2*k2, V, gm)
            k4 = model.rhs(times[i] + dt, s + dt*k3, V, g1)
            state[side] += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        if not np.all(np.isfinite(state)):
            raise RuntimeError("Forced time integration became nonfinite")
    output = np.asarray(outputs)
    return {"speed": V, "time": history, "leftHeave": output[:, 0, 0].tolist(),
            "leftTwist": output[:, 0, 1].tolist(), "rightHeave": output[:, 1, 0].tolist(),
            "rightTwist": output[:, 1, 1].tolist(), "gust": np.asarray(gust_outputs).tolist(),
            "stopped": stopped, "excitation": asdict(e)}
