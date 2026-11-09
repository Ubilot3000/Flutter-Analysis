"""Compatible entry point. Implementation lives in flutter_analysis/.

Run: python src/hstab_flutter.py --open
     python src/hstab_flutter.py --selftest
"""
from flutter_analysis import Params, HStab, sweep
from flutter_analysis.aerodynamics import JONES_A, JONES_BETA, make_polar, exact_C, jones_C
from flutter_analysis.shapes import _cantilever, _shape
from flutter_analysis.plotting import plots
from flutter_analysis.validation import selftest
from flutter_analysis.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
