"""Export an offline WebGL viewer with actual model time histories embedded."""
from dataclasses import asdict
from pathlib import Path
import json
import numpy as np
from .analysis import sweep
from .geometry import make_mesh
from .excitation import Excitation, export_kernel, simulate_tail


def _number(value):
    return float(value) if np.isfinite(value) else None


def export_viewer(model, V_test, vmax, outpath, critical=None, duration=6.0):
    p = model.p
    vf, vd = critical if critical is not None else (
        model.flutter_speed(vmax), model.divergence_speed(vmax))
    forcing = Excitation()
    response = simulate_tail(model, V_test, duration=duration, excitation=forcing)
    velocities = np.linspace(1, vmax, 120)
    frequency, damping = sweep(model, velocities)
    nullable = lambda a: [[_number(v) for v in row] for row in a]
    data = {"params": asdict(p), "span": p.span, "area": p.area,
            "aspectRatio": p.aspect_ratio, "mesh": make_mesh(p), "response": response,
            "kernel": export_kernel(model), "forcing": asdict(forcing), "duration": duration,
            "testSpeed": V_test, "vmax": vmax, "flutterSpeed": _number(vf),
            "divergenceSpeed": _number(vd), "sweep": {"speed": velocities.tolist(),
            "frequency": nullable(frequency), "damping": nullable(100 * damping)}}
    assets = Path(__file__).with_name("viewer")
    html = (assets / "index.html").read_text(encoding="utf-8")
    css = (assets / "style.css").read_text(encoding="utf-8") + (assets / "controls.css").read_text(encoding="utf-8")
    html = html.replace("/* VIEWER_CSS */", css)
    html = html.replace("/* SIMULATOR_JS */", (assets / "simulator.js").read_text(encoding="utf-8"))
    html = html.replace("/* VIEWER_JS */", (assets / "viewer.js").read_text(encoding="utf-8"))
    html = html.replace("/* MODEL_DATA */", json.dumps(data, allow_nan=False).replace("</", "<\\/"))
    outpath = Path(outpath)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    outpath.write_text(html, encoding="utf-8")
    return outpath
