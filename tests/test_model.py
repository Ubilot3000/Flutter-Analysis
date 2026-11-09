"""Numerical regression and geometry checks; run with unittest discovery."""
from dataclasses import asdict, replace
import json
from pathlib import Path
import re
import sys
import tempfile
import subprocess
import shutil
import unittest
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from flutter_analysis import Params, HStab
from flutter_analysis.aerodynamics import exact_C, jones_C
from flutter_analysis.geometry import thickness, make_mesh, deform
from flutter_analysis.visualization import export_viewer
from flutter_analysis.excitation import Excitation, gust_history, export_kernel, simulate_tail
from flutter_analysis.shapes import _shape


class ModelTests(unittest.TestCase):
    def test_supplied_geometry(self):
        p = Params()
        self.assertAlmostEqual(p.span, 0.74676)
        self.assertAlmostEqual(p.semispan, 0.37338)
        self.assertAlmostEqual(p.c_root, 0.256794)
        self.assertEqual(p.c_root, p.c_tip)
        self.assertAlmostEqual(p.area, p.span * p.c_root)

    def test_invalid_parameters(self):
        for kwargs in ({"mass_half": 0}, {"n_strips": 1}, {"semispan": float("nan")},
                       {"x_ea": 2}, {"twist_shape": "bad"}, {"zeta_roll": -1}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                Params(**kwargs)

    def test_mass_and_static_divergence(self):
        model = HStab(Params())
        self.assertTrue(np.all(np.linalg.eigvalsh(model.Ms) > 0))
        self.assertTrue(np.all(np.linalg.eigvalsh(model.Mtot) > 0))
        vd = model.divergence_speed(60)
        self.assertAlmostEqual(vd / model.divergence_speed_analytic(), 1, places=4)
        self.assertLess(model.divergence_margin(0.99 * vd), 0)
        self.assertGreater(model.divergence_margin(1.01 * vd), 0)

    def test_original_geometry_regression(self):
        # Original script's outputs for its placeholder geometry.
        model = HStab(replace(Params(), semispan=.30, c_root=.14, c_tip=.10, roll_shape="rigid_roll"))
        self.assertAlmostEqual(model.flutter_speed(40), 27.88, delta=.02)

    def test_theodorsen_approximation(self):
        k = np.linspace(.02, 3, 200)
        self.assertLess(np.max(np.abs(exact_C(k) - jones_C(k))), .02)

    def test_naca_and_antisymmetric_deformation(self):
        p = Params()
        x = np.linspace(0, 1, 1000)
        self.assertAlmostEqual(2 * thickness(x).max(), .12, delta=.0001)
        mesh = make_mesh(p, n_span=4, n_chord=12)
        v = np.asarray(mesh["vertices"])
        self.assertTrue(np.allclose(deform(v, 0, 0), v[:, :3]))
        shape = deform(v, .004, np.deg2rad(.5))
        count = len(v) // 2
        self.assertTrue(np.allclose(shape[:count, 1], -shape[count:, 1]))
        # Corresponding LE points move in opposite vertical directions.
        left_le = (4 * (2 * 12 + 1))
        right_le = left_le + count
        self.assertAlmostEqual(shape[left_le, 2], -shape[right_le, 2])
        self.assertGreater(shape[right_le, 2], .004)
        self.assertLess(np.max(mesh["triangles"]), len(v))

    def test_export_contains_actual_response_and_no_dependencies(self):
        model = HStab(Params())
        with tempfile.TemporaryDirectory() as directory:
            path = export_viewer(model, 15, 20, Path(directory) / "viewer.html",
                                 critical=(np.inf, np.inf), duration=.04)
            html = path.read_text(encoding="utf-8")
            payload = re.search(r'<script id="model-data" type="application/json">(.*?)</script>', html, re.S).group(1)
            data = json.loads(payload)
            expected = simulate_tail(model, 15, duration=.04)
            for name in ("leftHeave", "leftTwist", "rightHeave", "rightTwist"):
                np.testing.assert_allclose(data["response"][name], expected[name])
            self.assertNotIn("NaN", payload)
            self.assertNotRegex(html, r'<script[^>]*src=')
            self.assertNotIn("/* VIEWER_", html)

    def test_clamped_root_and_full_tail_forcing(self):
        p = Params()
        self.assertEqual(p.roll_shape, "cantilever")
        phi = _shape(p.roll_shape, np.array([0, 1e-6, 1]))
        self.assertEqual(phi[0], 0)
        self.assertLess(abs(phi[1] / 1e-6), 1e-5)
        self.assertAlmostEqual(phi[2], 1)
        model = HStab(p)
        for pattern, sign in (("symmetric", 1), ("antisymmetric", -1)):
            response = simulate_tail(model, 15, duration=.12, excitation=Excitation(pattern=pattern))
            np.testing.assert_allclose(response["leftHeave"], sign*np.array(response["rightHeave"]), atol=1e-13)
            self.assertGreater(np.std(response["rightHeave"]), 1e-5)
        response = simulate_tail(model, 15, duration=.12)
        self.assertGreater(np.max(np.abs(np.array(response["leftHeave"])-response["rightHeave"])), 1e-5)

    def test_gust_statistics_and_reproducibility(self):
        e = Excitation()
        _, gust = gust_history(15, 20, 1/2400, e)
        np.testing.assert_array_equal(gust[:100], gust_history(15, 99/2400, 1/2400, e)[1])
        np.testing.assert_allclose(gust.std(axis=0), e.gust_rms, rtol=.12)
        self.assertAlmostEqual(np.corrcoef(gust.T)[0,1], e.correlation, delta=.08)

    def test_rk4_matches_adaptive_unforced_solver(self):
        model = HStab(Params())
        r = simulate_tail(model, 15, duration=.2, excitation=Excitation(kind="release"))
        t, y = model.simulate(15, T=.2+1/1200, fs=600)
        np.testing.assert_allclose(r["rightHeave"], y[0], atol=2e-8)
        np.testing.assert_allclose(r["rightTwist"], y[1], atol=2e-8)

    @unittest.skipUnless(shutil.which("node"), "Node is needed for Python/browser solver comparison")
    def test_browser_solver_matches_python(self):
        model = HStab(Params())
        script = Path(__file__).with_name("compare_solver.mjs")
        options = [dict(speed=17.3, duration=.12, excitation=asdict(Excitation(pattern=pattern)))
                   for pattern in ("mixed", "symmetric", "antisymmetric")]
        options += [dict(speed=53.1, duration=.3, excitation=asdict(Excitation(kind="combined"))),
                    dict(speed=15, duration=.12, excitation=asdict(Excitation(kind="release")))]
        payload = {"kernel": export_kernel(model), "options": options,
                   "state": np.linspace(-.003,.003,model.n).tolist()}
        run = subprocess.run(["node", str(script)], input=json.dumps(payload), capture_output=True, text=True, check=True)
        data = json.loads(run.stdout)
        np.testing.assert_allclose(data["rhs"], model.rhs(0, np.asarray(payload["state"]), 17.3, .27), atol=1e-10)
        for opt, actual in zip(options, data["results"]):
            expected = simulate_tail(model, opt["speed"], duration=opt["duration"], excitation=Excitation(**opt["excitation"]))
            for name in ("leftHeave","leftTwist","rightHeave","rightTwist","gust"):
                np.testing.assert_allclose(actual[name], expected[name], atol=2e-11)

    def test_zero_input_stays_still_and_large_response_stops(self):
        model = HStab(Params())
        r = simulate_tail(model, 15, duration=.1, excitation=Excitation(gust_rms=0))
        self.assertTrue(np.all(np.asarray(r["rightHeave"]) == 0))
        r = simulate_tail(model, 60, duration=2, excitation=Excitation(kind="release"))
        self.assertTrue(r["stopped"])
        self.assertLess(r["time"][-1], 2)


if __name__ == "__main__":
    unittest.main()
