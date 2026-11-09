"""NACA 00xx display mesh and antisymmetric modal deformation.

Axes: x aft, y spanwise, z up; elastic axis at x=0. Positive twist
lifts the leading edge. Display thickness does not determine stiffness.
"""
import numpy as np
from .shapes import _shape


def thickness(x, ratio=0.12):
    """Classic NACA four-digit half-thickness / chord, finite trailing edge."""
    x = np.asarray(x)
    return 5 * ratio * (0.2969 * np.sqrt(x) - 0.1260 * x - 0.3516 * x ** 2
                        + 0.2843 * x ** 3 - 0.1015 * x ** 4)


def make_mesh(p, n_span=20, n_chord=40):
    """Indexed surface mesh; each vertex stores x,y,z,roll,twist,eta."""
    vertices, triangles = [], []
    u = (1 - np.cos(np.linspace(0, np.pi, n_chord + 1))) / 2
    section = [(v, float(thickness(v, p.thickness_ratio))) for v in u]
    section += [(v, -float(thickness(v, p.thickness_ratio))) for v in u[:0:-1]]
    count = len(section)
    for side in (-1, 1):
        offset = len(vertices)
        for eta in np.linspace(0, 1, n_span + 1):
            chord = p.c_root + (p.c_tip - p.c_root) * eta
            ph = side * float(_shape(p.roll_shape, np.asarray(eta)))
            ps = side * float(_shape(p.twist_shape, np.asarray(eta)))
            for x, z in section:
                vertices.append([(x - p.x_ea) * chord, side * eta * p.semispan,
                                 z * chord, ph, ps, eta])
        for j in range(n_span):
            for k in range(count):
                a = offset + j * count + k
                b = offset + j * count + (k + 1) % count
                triangles.extend([[a, b, a + count], [b, b + count, a + count]])
        for j in (0, n_span):
            start = offset + j * count
            for k in range(1, count - 1):
                triangles.append([start, start + k, start + k + 1])
    return {"vertices": vertices, "triangles": triangles}


def deform(vertices, q_h, q_twist, amplification=1.0):
    """Reference deformation used to verify the viewer transform."""
    v = np.asarray(vertices)
    theta = v[:, 4] * q_twist * amplification
    return np.column_stack((v[:, 0] * np.cos(theta) + v[:, 2] * np.sin(theta),
                            v[:, 1],
                            -v[:, 0] * np.sin(theta) + v[:, 2] * np.cos(theta)
                            + v[:, 3] * q_h * amplification))
