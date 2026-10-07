"""Acceptance checks for SPEC_10 v0.2 Step 1: the fit of the IRIS temperatures and their gradient.

Reads the committed `occul_data/lindal/iris_temperatures.csv` and writes only under
`reports/step10_1/`: one comparison figure per level and `comparison.json`. The author views the
comparison and chooses the degree and minimum that become `fit`'s defaults.

Every check prints its measured values. Run from the repository root.
"""

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from casspian.tools.lindal import iris_temperatures as iris

HERE = Path("reports/step10_1")
HERE.mkdir(parents=True, exist_ok=True)
TABLE = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
GRID = np.arange(90.0, -90.01, -0.5)
CANDIDATES = {"a": {"degree": 1, "min_points": 3, "label": "(a) local linear, at least 3"},
              "b": {"degree": 1, "min_points": 8, "label": "(b) local linear, at least 8"},
              "c": {"degree": 2, "min_points": 4, "label": "(c) local quadratic, at least 4"}}
COLORS = {"a": "C0", "b": "C3", "c": "C2"}
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def level(pressure):
    rows = TABLE["pressure_mbar"] == pressure
    return TABLE["latitude_planetographic_deg"][rows], TABLE["temperature_K"][rows]


# ---------------------------------------------------------------------------
# 0. The fit on fields it must reproduce (beyond the specification)
# ---------------------------------------------------------------------------
lat, _ = level(110)
probe = np.array([60.0, 20.0, 0.0, -40.0, -80.0])
lines, ok = [], True
for name, spec in CANDIDATES.items():
    for field, f, df, exact in (("linear", lambda p: 90.0 + 0.05 * p, lambda p: 0.05 + 0 * p, True),
                                ("quadratic", lambda p: 90.0 + 0.05 * p - 0.002 * p ** 2,
                                 lambda p: 0.05 - 0.004 * p, spec["degree"] == 2)):
        result = iris.fit(lat, f(lat), spec["degree"], spec["min_points"], grid=probe)
        value_error = float(np.max(np.abs(result.value_K - f(probe))))
        slope_error = float(np.max(np.abs(result.gradient_K_per_deg - df(probe))))
        if exact:
            good = value_error <= 1e-9 and slope_error <= 1e-11
            ok &= good
            lines.append(f"{spec['label']}, {field} field: |value error| {value_error:.1e} K, |slope "
                         f"error| {slope_error:.1e} K/deg (exact by construction; bound 1e-9, 1e-11)")
# The slope's standard error against the scatter of slopes over noisy realizations, report only.
rng = np.random.default_rng(10)
noise = 0.5
for name, spec in CANDIDATES.items():
    slopes, predicted = [], []
    for _ in range(300):
        result = iris.fit(lat, 85.0 + rng.normal(0.0, noise, lat.size), spec["degree"],
                          spec["min_points"], grid=probe)
        slopes.append(result.gradient_K_per_deg)
        predicted.append(result.gradient_error_K_per_deg)
    ratio = np.mean(predicted, axis=0) / np.std(slopes, axis=0)
    lines.append(f"{spec['label']}: mean predicted slope error over the scatter of slopes, 300 "
                 f"realizations of {noise} K noise on the 110 mbar latitudes, at {probe.tolist()} deg: "
                 f"{np.round(ratio, 2).tolist()} (reported, not bounded)")
record(0, "beyond the specification: the fit returns a linear field's value and slope exactly, and a "
          "local quadratic a quadratic's; the standard error against Monte Carlo scatter",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 1. The comparison
# ---------------------------------------------------------------------------
lines, summary = [], {}
for pressure in (110, 290, 730):
    lat, temp = level(pressure)
    fits = {name: iris.fit(lat, temp, s["degree"], s["min_points"], grid=GRID)
            for name, s in CANDIDATES.items()}
    at_points = {name: iris.fit(lat, temp, s["degree"], s["min_points"], grid=lat)
                 for name, s in CANDIDATES.items()}
    fig, (top, bottom) = plt.subplots(2, 1, figsize=(11, 8.5), sharex=True,
                                      gridspec_kw={"height_ratios": [1.2, 1.0]})
    top.plot(lat, temp, ".", color="0.4", markersize=3, label=f"{pressure} mbar, digitized")
    summary[pressure] = {}
    for name, s in CANDIDATES.items():
        f = fits[name]
        rms = float(np.sqrt(np.mean((temp - at_points[name].value_K) ** 2)))
        widened = np.isfinite(f.fwhm_deg) & (f.fwhm_deg > 4.0 + 1e-9)
        summary[pressure][name] = {
            "rms_K": rms,
            "widened_fraction": float(widened.sum() / np.isfinite(f.fwhm_deg).sum()),
            "widest_fwhm_deg": float(np.nanmax(f.fwhm_deg)),
            "largest_gradient_K_per_deg": float(np.nanmax(np.abs(f.gradient_K_per_deg))),
            "median_gradient_error_K_per_deg": float(np.nanmedian(f.gradient_error_K_per_deg))}
        top.plot(GRID, f.value_K, color=COLORS[name], linewidth=1.1,
                 label=f"{s['label']}; RMS {rms:.3f} K")
        bottom.plot(GRID, f.gradient_K_per_deg, color=COLORS[name], linewidth=1.1, label=s["label"])
        lines.append(f"{pressure} mbar, {s['label']}: RMS of the points about the fit {rms:.3f} K; "
                     f"window widened at {100 * summary[pressure][name]['widened_fraction']:.0f} percent "
                     f"of the grid, to at most {summary[pressure][name]['widest_fwhm_deg']:.1f} deg FWHM; "
                     f"largest |dT/dphi| {summary[pressure][name]['largest_gradient_K_per_deg']:.3f} K/deg; "
                     f"median standard error {summary[pressure][name]['median_gradient_error_K_per_deg']:.3f} K/deg")
    a = fits["a"]
    bottom.fill_between(GRID, a.gradient_K_per_deg - a.gradient_error_K_per_deg,
                        a.gradient_K_per_deg + a.gradient_error_K_per_deg, color=COLORS["a"],
                        alpha=0.2, linewidth=0, label="(a) +- 1 standard error")
    bottom.axhline(0.0, color="0.6", linewidth=0.8)
    top.set_ylabel("temperature (K)")
    bottom.set_ylabel("dT/dphi (K per deg)")
    bottom.set_xlabel("planetographic latitude (deg)")
    bottom.set_xlim(90, -90)
    top.legend(fontsize=8, loc="best")
    bottom.legend(fontsize=8, loc="best")
    for ax in (top, bottom):
        ax.grid(alpha=0.3)
    top.set_title(f"IRIS temperatures at {pressure} mbar: three fits, 4 deg FWHM, window widened to "
                  f"the minimum effective points")
    fig.tight_layout()
    fig.savefig(HERE / f"comparison_{pressure}mbar.png", dpi=130)
    plt.close(fig)
(HERE / "comparison.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
figures = sorted(str(p) for p in HERE.glob("comparison_*mbar.png"))
record(1, "the comparison at each level on a 0.5 deg grid: the points with fits (a), (b), (c), their "
          "gradients with (a)'s standard error, the RMS about each; for the author to choose",
       len(figures) == 3, "\n".join(lines) + "\nfigures " + ", ".join(figures))

failed = [r for r in results if not r[2]]
summary_line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary_line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary_line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
