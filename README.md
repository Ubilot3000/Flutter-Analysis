# Horizontal-stabilizer flutter analysis

The horizontal tail uses the supplied **NACA 0012**, **29.40 in full span** and
**10.11 in chord**. The planform is assumed rectangular and unswept because no
taper or sweep was specified. The viewer models two 14.70 in semispans, each with
its own bending and twist coordinates and a fixed center attachment. The default
bending shape has zero displacement and zero slope at the root. The vertical tail
is outside its scope.

## Run on Windows

From this directory, using the existing virtual environment:

```powershell
.\venv\Scripts\python.exe src/hstab_flutter.py --open
```

This generates `output/hstab_3d.html` and `output/hstab_flutter.png`, then opens the
interactive viewer. You can also open the HTML file directly; it works offline
with no CDN, JavaScript packages, server, or Python process running.

For a fresh environment:

```powershell
python -m venv venv
.\venv\Scripts\python.exe -m pip install -r requirements.txt
.\venv\Scripts\python.exe src/hstab_flutter.py --open
```

Useful options:

```powershell
# Change test speed, sweep range, animation duration, and destination.
.\venv\Scripts\python.exe src/hstab_flutter.py --V 25 --vmax 60 --duration 3 --out output
# Generate just the interactive viewer.
.\venv\Scripts\python.exe src/hstab_flutter.py --no-plots
# Physical validation and automated regression checks.
.\venv\Scripts\python.exe src/hstab_flutter.py --selftest
.\venv\Scripts\python.exe -m unittest discover -s tests -v
```

An optional browser integration check uses Node 22+ and an installed Chrome:
`node tests/check_viewer.mjs`. Set `CHROME_PATH` if Chrome is installed elsewhere.
It opens a headless browser on the local HTML file, verifies the interactive
controls and simulated values, and writes desktop/mobile previews under `output/`.
The Python suite also compares the browser integrator with Python when Node is
available; that comparison is skipped if Node is absent.

## Interactive 3D view

- Drag to orbit; scroll to zoom. The focused canvas also supports arrow keys and
  plus/minus. Choose perspective, top, front, or section views.
- Play, pause, scrub, restart, and loop the simulated time history. Playback starts
  slowed to 0.5 times real time so the vibration is visible. Reduced-motion
  preferences disable autoplay.
- Use the **mean airspeed slider** in 0.1 m/s increments. Each selected speed runs
  a fresh integration in a Web Worker inside the offline HTML, then restarts the
  record. It does not interpolate between precomputed speed cases. Mean airspeed
  remains constant within a record; the vertical gust speed varies with time.
- Choose continuous random gusts, periodic airflow, their combination, or release
  and decay. Adjust vertical gust RMS and periodic frequency/peak speed.
- Choose shared + local gusts (default), symmetric airflow (both tips together),
  or antisymmetric airflow (opposite tips). Mixed-mode forcing does not impose a
  permanent up/down direction on either tip.
- Show coupled motion, bending only, or twist only. These filters isolate coordinates
  from the coupled solution; they do not rerun uncoupled physics.
- Amplify the displayed motion from 1 to 20 times. Tip readouts retain physical
  values. The reference outline, elastic axis, and CG can be toggled.

The HTML embeds an initial Python-computed forced response, the half-tail strip
equations and a linear damping sweep. Its browser integrator uses RK4 with four
substeps per 600 Hz output interval; equation and trajectory parity are checked
against Python. It re-integrates both halves when speed or excitation changes.
Re-run Python after changing structural parameters. The PNG report retains the
unforced release response and linear stability sweep.

The default six-second record repeats when Loop is enabled. This is replay of a
finite gust realization; the loop boundary resets the state. A response that
reaches 10% semispan displacement or 10 degrees twist stops and disables automatic
looping. That boundary marks the small-motion model's range, not a physical
amplitude limit or a claim that the structure would survive larger motion.

### Representation of changing airflow

`excitation.py` generates a repeatable first-order colored random **vertical gust
velocity** with an Ornstein-Uhlenbeck process. Its ensemble RMS defaults to
0.25 m/s and its correlation time is `length_scale / airspeed`, using an
illustrative 0.5 m length scale. The mixed pattern combines a shared process and
independent local processes, giving left/right gust correlation 0.75. This makes
ongoing, irregular shaking possible well below flutter. Finite records can have
a nonzero sample mean; no artificial recentering is applied to the response.

This is a demonstration model, not a calibrated Dryden/von Karman atmosphere or
a measured propeller wake. Gust speed is added to the existing downwash and wake
lags; a separate Sears gust-admittance model is not implemented. Periodic airflow
is an optional sinusoidal vertical gust (default 20 Hz, 0.10 m/s peak), with no
assumed relationship to motor RPM. Actual turbulence spectra, wake harmonics and
mount acceleration should come from flight or bench measurements. The current
halves share airflow and a fixed root, with no flexible-boom/fuselage coupling.

## Inputs and interpretation

Edit `src/flutter_analysis/parameters.py`. All model dimensions are SI internally:
span 0.74676 m, semispan 0.37338 m, chord 0.256794 m. Inch conversion happens only in
the parameter defaults and display.

Only geometry and the airfoil designation came from the supplied table. Original
illustrative values are retained for **half-tail mass (50 g), dry uncoupled modal
frequencies (12 / 30 Hz), damping (1%), elastic axis (30% chord), CG (38% chord),
radius of gyration and lift polar**. They must be replaced with measured or
independently calculated properties. The mesh's 12% NACA thickness does not infer
material, spar construction, stiffness, or inertia.

For those assumptions, the default 1–60 m/s sweep finds no oscillatory flutter
crossing and finds static torsional divergence near **50.56 m/s**. This is a model
result, not a validated operating limit. A flutter crossing at a higher speed would
occur beyond an already statically unstable condition. The viewer labels both
instabilities separately. Missing oscillatory branches are left blank on the
damping plot, not treated as stable zero damping.

Each half's structural coordinates are tip displacement (m) and tip nose-up
twist (rad). The original `f_roll` input name is retained for compatibility and
now sets the uncoupled bending frequency. Set `roll_shape="rigid_roll"` explicitly
to study the former compliant-root assumption. Aerodynamics use spanwise strips,
Theodorsen apparent mass, and a two-lag Jones approximation for circulatory lift.
The `tanh` polar is illustrative; it is not a measured NACA 0012 low-Reynolds-number
polar, a 3D CFD result, or a validated large-deflection/limit-cycle model. The
default stability sweep linearizes about zero displacement and zero trim angle.

References: [NASA's two-degree-of-freedom flutter method](https://ntrs.nasa.gov/citations/19850014068)
[NASA airfoil geometry report](https://ntrs.nasa.gov/api/citations/19970008124/downloads/19970008124.pdf),
and [turbulence spectra and advection](https://www.mathworks.com/help/aeroblks/drydenwindturbulencemodelcontinuous.html).

## Code layout

| File under `src/flutter_analysis/` | Purpose |
| --- | --- |
| `parameters.py` | Dimensions, units, physical inputs, input validation |
| `aerodynamics.py` | Polar, Jones constants, Theodorsen functions |
| `shapes.py` | Tip-normalized roll/bending and twist shapes |
| `model.py` | Mass/stiffness/damping, state equations, stability search, integration |
| `excitation.py` | Correlated continuous gusts, optional wake tone, two-half integration |
| `analysis.py` | Eigenvalue branch tracking |
| `geometry.py` | NACA surface mesh and modal deformation |
| `plotting.py` | Static PNG report |
| `visualization.py` | Initial response, equation coefficients, standalone HTML export |
| `viewer/` | HTML, CSS, WebGL rendering, background response integrator |
| `validation.py` | Physical self-tests |
| `cli.py` | Command-line orchestration |

`src/hstab_flutter.py` remains the entry point and re-exports the original public
model functions for existing imports. `Params`, `HStab`, and `sweep` can also be
imported from `flutter_analysis` when `src` is on the Python path.
