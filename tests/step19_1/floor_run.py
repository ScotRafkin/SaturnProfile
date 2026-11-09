"""SPEC_19 v0.5 section 5 item 5: `casspian-wind-fit` once, the three runs rebuilt from its wind, and what the
closing asks for: the passes and the group ratios (the command prints them) and delivered minus observed at
110 mbar, read from the rebuilt runs. Also the build for the clean rebuild at the acceptance commit.

Run from the repository root. About an hour alone.
"""

import subprocess
import sys
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
import accept_step19_1 as accept  # noqa: E402
from casspian.tools.lindal import radio_profiles as rp  # noqa: E402
from casspian.tools.wind.estimation import driver  # noqa: E402
from casspian.tools.wind.estimation.radio import interp_ln_p  # noqa: E402

#: Before the floor (REPORT_19_step1_v04 section 4), and the closing's bound on the change at 31.2 S.
BEFORE_31S_K = 3.34
BOUND_K = 1.0


def main():
    warnings.simplefilter("ignore")
    subprocess.run(["casspian-wind-fit", str(accept.CONTROL)], check=True)
    problem = driver.read_control(accept.CONTROL)
    wind_path = problem.output.relative_to(Path.cwd())
    at = np.array([accept.AT_110_Pa])
    moved = None
    for target in problem.radio.targets:
        name = accept.RUNS[target.name][2]
        product = accept.make_run(name, accept.RUNS[target.name][0], wind_path)
        root = driver.tree(product).to_dataset(inherit=False)
        T = interp_ln_p(root["pressure_Pa"].values, root["temperature_K"].values, at)[0]
        table = rp.read_csv(target.observed)
        value = T - interp_ln_p(table["pressure_mbar"] * 100.0, table["temperature_K"], at)[0]
        print(f"{target.name}: delivered minus observed at 110 mbar {value:+.3f} K; pressure identity largest "
              f"{float(np.nanmax(np.abs(root['pressure_identity_residual'].values))):.3e}; run {name}", flush=True)
        if target.name == "31.2 S":
            moved = value - BEFORE_31S_K
    print(f"31.2 S moved {moved:+.3f} K from {BEFORE_31S_K:+.2f} K; bound {BOUND_K} K: {'within' if abs(moved) <= BOUND_K else 'BEYOND, stop and report'}")
    return 0 if abs(moved) <= BOUND_K else 1


if __name__ == "__main__":
    raise SystemExit(main())
