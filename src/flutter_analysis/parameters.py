"""SI model inputs. Only the planform dimensions come from the supplied table."""
from dataclasses import dataclass
import math

INCH = 0.0254


@dataclass(frozen=True)
class Params:
    # 29.40 in is interpreted as FULL span; no taper or sweep was specified.
    semispan: float = 29.40 * INCH / 2
    c_root: float = 10.11 * INCH
    c_tip: float = 10.11 * INCH
    airfoil: str = "NACA 0012"
    thickness_ratio: float = 0.12
    # Remaining values are original illustrative assumptions, not measurements.
    x_ea: float = 0.30       # elastic axis, fraction of chord from LE
    x_cg: float = 0.38       # section cg, fraction of chord from LE
    mass_half: float = 0.050 # kg, ONE half including elevator and hardware
    k_gyr: float = 0.29      # cg radius of gyration / chord
    f_roll: float = 12.0     # Hz, uncoupled bending frequency (legacy name)
    f_twist: float = 30.0    # Hz, uncoupled torsion frequency
    zeta_roll: float = 0.01
    zeta_twist: float = 0.01
    roll_shape: str = "cantilever"  # clamped center: zero displacement AND bending slope
    twist_shape: str = "linear"
    rho: float = 1.225
    cl_alpha: float = 5.5    # assumed effective slope, 1/rad
    cl_max: float = 1.0     # illustrative tanh polar ceiling
    alpha0_deg: float = 0.0
    n_strips: int = 12

    @property
    def span(self):
        return 2 * self.semispan

    @property
    def area(self):
        return self.semispan * (self.c_root + self.c_tip)

    @property
    def aspect_ratio(self):
        return self.span ** 2 / self.area

    def __post_init__(self):
        for name in ("semispan", "c_root", "c_tip", "mass_half", "k_gyr",
                     "f_roll", "f_twist", "rho", "cl_alpha", "cl_max", "thickness_ratio"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0:
                raise ValueError(f"{name} must be finite and positive")
        for name in ("zeta_roll", "zeta_twist"):
            if not math.isfinite(getattr(self, name)) or getattr(self, name) < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name in ("x_ea", "x_cg"):
            if not math.isfinite(getattr(self, name)) or not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be a chord fraction from 0 to 1")
        if not math.isfinite(self.alpha0_deg):
            raise ValueError("alpha0_deg must be finite")
        if isinstance(self.n_strips, bool) or not isinstance(self.n_strips, int) or self.n_strips < 2:
            raise ValueError("n_strips must be an integer >= 2")
        if self.roll_shape not in ("rigid_roll", "cantilever"):
            raise ValueError("roll_shape must be rigid_roll or cantilever")
        if self.twist_shape not in ("rigid_pitch", "linear", "sine"):
            raise ValueError("twist_shape must be rigid_pitch, linear, or sine")
