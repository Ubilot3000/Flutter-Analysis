"""Static stability sweep and response report."""
from pathlib import Path
import numpy as np
from .analysis import sweep


def plots(model, V_test, vmax, outdir=".", critical=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    speeds = np.arange(1, vmax + 0.001, 0.5)
    frequency, damping = sweep(model, speeds)
    vf, vd = critical if critical is not None else (
        model.flutter_speed(vmax), model.divergence_speed(vmax))
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    colors = ("#2a78d6", "#eb6834")
    names = ("Lower-frequency branch", "Higher-frequency branch")
    for k, (color, name) in enumerate(zip(colors, names)):
        axes[0, 0].plot(speeds, frequency[:, k], color=color, label=name)
        axes[0, 1].plot(speeds, 100 * damping[:, k], color=color)
    axes[0, 0].set(ylabel="Frequency [Hz]", title="Oscillatory branches")
    axes[0, 0].legend(frameon=False, fontsize=9)
    axes[0, 1].set(ylabel="Damping ratio [%]", title="Negative damping indicates oscillatory flutter")
    axes[0, 1].axhline(0, color="#666666", lw=1)
    for ax in axes[0]:
        ax.set_xlabel("Airspeed [m/s]")
        ax.axvline(V_test, color="#777777", ls=":", lw=1)
        for value, color, label in ((vf, "#bd3c32", "Flutter"), (vd, "#8053b5", "Divergence")):
            if np.isfinite(value):
                ax.axvline(value, color=color, ls="--", lw=1)
                ax.text(value, .97, f" {label} {value:.1f}", color=color,
                        transform=ax.get_xaxis_transform(), rotation=90, va="top", fontsize=8)
    t, state = model.simulate(V_test)
    axes[1, 0].plot(t, 1000 * state[0], color=colors[0])
    axes[1, 0].set(xlabel="Time [s]", ylabel="Right tip displacement [mm]", title=f"Response at {V_test:g} m/s")
    axes[1, 1].plot(t, np.rad2deg(state[1]), color=colors[1])
    axes[1, 1].set(xlabel="Time [s]", ylabel="Right tip twist [deg]", title="Coupled torsional response")
    for ax in axes.flat:
        ax.grid(color="#e4e3df", lw=.6)
        ax.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Horizontal-tail model - assumed structural properties", fontsize=13)
    fig.savefig(outdir / "hstab_flutter.png", dpi=150)
    plt.close(fig)
    return vf, vd
