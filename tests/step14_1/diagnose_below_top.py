"""Diagnosis for SPEC_14 Step 1 check 1: why the wind below 110 mbar is not the base's bit for bit.

The construction iterates until the wind changes nowhere, at any pressure, by more than
`TOLERANCE_MS`. Case (iii) changes the wind aloft, so the stopping test can end the two cases on
different iterations, and the shear at the three levels (which depends only on the wind at the
levels) is then taken from different iterates. This script measures that: the changes per
iteration in each case, the difference below 110 mbar at the shipped tolerance, and the same with
the tolerance lowered in this process only, to show the difference is the stopping rule's. No
code is changed. Writes `reports/step14_1/diagnose_below_top.txt`.

Run from the repository root.
"""

import tomllib
import warnings
from pathlib import Path

import numpy as np

from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.wind import shear as shear_tool

HERE = Path("reports/step14_1")
RUNS = ("lindal_iris_v_sl_ii_110", "lindal_iris_v_sl_iii_110")
TOP_Pa = 11000.0
warnings.simplefilter("ignore")


def build(name):
    run = Path("forward") / name
    with open(run / f"{name}_build.toml", "rb") as handle:
        section = tomllib.load(handle)["shear"]
    inputs = {key: shear_tool._read(run / section[key], kind) for key, kind in shear_tool.INPUT_KEYS.items()}
    source = shear_tool._read(run / section["source"], "wind")
    return lw.construct(source, inputs["temperatures"], inputs["composition"], inputs["gravity"],
                        inputs["rotation"], float(section["top_level_Pa"]), section["above_top"],
                        gradient=section.get("gradient", "slope"), window=section.get("window", "adaptive"),
                        hold_Pa=section.get("hold_Pa"), cross_Pa=section.get("cross_Pa"))


lines = []
shipped = lw.TOLERANCE_MS
for tolerance in (shipped, 1e-6):
    lw.TOLERANCE_MS = tolerance            # this process only; the module's value is unchanged on disk
    built = {name: build(name) for name in RUNS}
    a, b = (built[n] for n in RUNS)
    p = lw.pressure_grid()
    below = p >= TOP_Pa
    u_a = np.asarray(a.dataset["u_total_ms"].values)
    u_b = np.asarray(b.dataset["u_total_ms"].values)
    du = np.abs(u_b[:, below] - u_a[:, below])
    ds = np.abs(b.shear_ms - a.shear_ms)
    lines.append(f"tolerance {tolerance:g} m/s{' (shipped)' if tolerance == shipped else ' (diagnostic only)'}:")
    for name in RUNS:
        c = built[name]
        lines.append(f"  {name}: {c.iterations} iterations, changes (m/s) "
                     + ", ".join(f"{x:.3g}" for x in c.changes_ms))
    lines.append(f"  below 110 mbar, |u_total(iii) - u_total(ii)|: largest {du.max():.3e} m/s, "
                 f"bit-identical {bool(np.array_equal(u_a[:, below], u_b[:, below]))}")
    lines.append(f"  the shear at the three levels, |s(iii) - s(ii)|: largest {ds.max():.3e} m/s")
lw.TOLERANCE_MS = shipped
text = "\n".join(lines)
print(text)
(HERE / "diagnose_below_top.txt").write_text(text + "\n", encoding="utf-8", newline="\n")
