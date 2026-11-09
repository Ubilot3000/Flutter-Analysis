"""Command-line analysis, static plots, and interactive 3D export."""
import argparse
from pathlib import Path
import webbrowser
import numpy as np
from .parameters import Params
from .model import HStab


def positive(value):
    number = float(value)
    if not np.isfinite(number) or number <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return number


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--V", type=positive, default=15.0, help="test airspeed [m/s]")
    ap.add_argument("--vmax", type=positive, default=60.0, help="stability search upper speed [m/s]")
    ap.add_argument("--out", default="output", help="output directory (created if needed)")
    ap.add_argument("--duration", type=positive, default=6.0, help="3D response duration [s]")
    ap.add_argument("--no-plots", action="store_true", help="skip the static PNG report")
    ap.add_argument("--no-viewer", action="store_true", help="skip the interactive HTML report")
    ap.add_argument("--open", action="store_true", help="open the 3D viewer in your browser")
    ap.add_argument("--selftest", action="store_true", help="run physical validation checks")
    args = ap.parse_args(argv)
    if args.vmax < 1:
        ap.error("--vmax must be >= 1 m/s")
    if args.open and args.no_viewer:
        ap.error("--open requires the viewer; remove --no-viewer")
    p = Params()
    if args.selftest:
        from .validation import selftest
        return 0 if selftest(p) else 1
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    model = HStab(p)
    print(f"{p.airfoil}: full span {p.span:.5f} m, chord {p.c_root:.6f} m, AR {p.aspect_ratio:.2f}")
    print("Mass, modal frequencies, damping, CG and elastic axis are illustrative assumptions.")
    model.summary()
    vf = model.flutter_speed(args.vmax)
    vd = model.divergence_speed(args.vmax)
    for name, value in (("Flutter", vf), ("Divergence", vd)):
        print(f"{name} speed: {value:.2f} m/s" if np.isfinite(value)
              else f"No {name.lower()} crossing found in 1-{args.vmax:g} m/s")
    if not args.no_plots:
        from .plotting import plots
        plots(model, args.V, args.vmax, out, critical=(vf, vd))
        print(f"Static report: {(out / 'hstab_flutter.png').resolve()}")
    if not args.no_viewer:
        from .visualization import export_viewer
        viewer = export_viewer(model, args.V, args.vmax, out / "hstab_3d.html",
                               critical=(vf, vd), duration=args.duration)
        print(f"Interactive 3D: {viewer.resolve()}")
        if args.open:
            webbrowser.open(viewer.resolve().as_uri())
    return 0
