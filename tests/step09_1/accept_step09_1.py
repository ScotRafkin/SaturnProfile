"""Acceptance checks for SPEC_09 v0.3 Step 1: the IRIS temperatures at three levels, digitized.

Runs the digitization of `casspian.tools.lindal.iris_temperatures` on
`docs/conrath_etal_IRIS_1983.pdf` and writes only under `reports/step09_1/`: the reproduced CSV,
the spread figure and the overlay. The committed `occul_data/lindal/iris_temperatures.csv` is
read, never written.

Every check prints its measured values. Run from the repository root.
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw

from casspian.tools.lindal import iris_temperatures as iris

HERE = Path("reports/step09_1")
HERE.mkdir(parents=True, exist_ok=True)
PDF = Path("docs/conrath_etal_IRIS_1983.pdf")
COMMITTED = Path("occul_data/lindal/iris_temperatures.csv")
FWHM_DEG = 4.0
COLORS = {730: (200, 0, 0), 290: (0, 140, 0), 110: (0, 0, 220)}
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


result = iris.digitize(PDF)
table = np.genfromtxt(COMMITTED, delimiter=",", names=True)

# ---------------------------------------------------------------------------
# 0. The committed CSV is what the tool writes
# ---------------------------------------------------------------------------
again = iris.write_csv(result, HERE / "iris_temperatures.csv")
same = again.read_bytes() == COMMITTED.read_bytes()
counts = {p: int((table["pressure_mbar"] == p).sum()) for p in (110, 290, 730)}
record(0, "running the tool reproduces the committed CSV byte for byte",
       same,
       f"{table.size} rows, by level {counts}; identical {same}\n"
       f"frame (top, bottom, left, right) {result.frame} px at {iris.DPI} dpi; single dot "
       f"{result.single_dot_area_px:g} px; {result.dropped_curve} remnant of the broken curve "
       f"dropped; {int((result.k > 1).sum())} dots split by k-means from components holding "
       f"2 to {int(result.k.max())} dots\n"
       f"temperature: {iris.TEMPERATURE_TICKS_K[0]:g} to {iris.TEMPERATURE_TICKS_K[-1]:g} K, "
       f"{1.0 / abs(result.temperature_fit[0]):.2f} px per K, tick residuals at most "
       f"{np.abs(result.temperature_residual_K).max():.3f} K\n"
       f"latitude: {iris.LATITUDE_TICKS_DEG[0]:g} to {iris.LATITUDE_TICKS_DEG[-1]:g} deg, "
       f"{1.0 / abs(result.latitude_fit[0]):.2f} px per deg, tick residuals at most "
       f"{np.abs(result.latitude_residual_deg).max():.3f} deg")

# ---------------------------------------------------------------------------
# 1. The spread about the running mean
# ---------------------------------------------------------------------------
lines, ok = [], True
fig, axes = plt.subplots(3, 1, figsize=(10, 11), sharex=True)
for ax, level in zip(axes, (110, 290, 730)):
    rows = table["pressure_mbar"] == level
    lat, temp = table["latitude_planetographic_deg"][rows], table["temperature_K"][rows]
    fit_at_dots = iris.running_mean(lat, temp, lat, FWHM_DEG)
    residual = temp - fit_at_dots
    rms = float(np.sqrt(np.mean(residual ** 2)))
    sigma_read = float(np.sqrt(np.mean((temp - iris.running_mean(
        lat, temp, lat, 4.0 * np.sqrt(8.0 * np.log(2.0)))) ** 2)))
    grid = np.arange(-90.0, 90.01, 0.25)
    fit = iris.running_mean(lat, temp, grid, FWHM_DEG)
    ok &= np.isfinite(rms) and np.isfinite(fit_at_dots).all()
    gaps = np.flatnonzero(np.diff(np.isfinite(fit).astype(int)) != 0)
    lines.append(f"{level} mbar: {lat.size} dots, {lat.max():.2f} to {lat.min():.2f} deg; RMS about the "
                 f"running mean (FWHM {FWHM_DEG:g} deg) {rms:.3f} K; with 4 deg read as the standard "
                 f"deviation instead, {sigma_read:.3f} K; largest |residual| {np.abs(residual).max():.2f} K; "
                 f"the fit undefined (no dot within {FWHM_DEG:g} deg) across "
                 f"{', '.join(f'{grid[a + 1]:g}' for a in gaps)} deg")
    ax.plot(lat, temp, ".", color="0.35", markersize=3, label=f"{level} mbar, digitized")
    ax.plot(grid, fit, color="C3", linewidth=1.2,
            label=f"Gaussian running mean, FWHM {FWHM_DEG:g} deg; RMS {rms:.2f} K")
    ax.set_ylabel("temperature (K)")
    ax.legend(loc="best", fontsize=8)
    ax.grid(alpha=0.3)
axes[0].set_title("Conrath and Pirraglia (1983) Fig. 1, digitized (top level written 110 mbar)")
axes[-1].set_xlabel("planetographic latitude (deg)")
axes[-1].set_xlim(90, -90)
fig.tight_layout()
spread = HERE / "spread.png"
fig.savefig(spread, dpi=130)
plt.close(fig)
record(1, "per level, the RMS of the dots about a Gaussian-weighted running mean in latitude, the fit "
          "drawn over the dots", ok, "\n".join(lines) + f"\nfigure {spread}")

# ---------------------------------------------------------------------------
# 2. The overlay
# ---------------------------------------------------------------------------
top, bottom, left, right = result.frame
image = Image.fromarray(np.where(result.ink, 0, 255).astype(np.uint8)).convert("RGB")
draw = ImageDraw.Draw(image)
for row, col, level in zip(result.row, result.col, result.pressure_mbar):
    draw.ellipse([col - 4, row - 4, col + 4, row + 4], fill=COLORS[int(level)])
for box in iris.LABEL_BOXES.values():
    draw.rectangle(box, outline=(255, 150, 0), width=3)
overlay = HERE / "overlay.png"
image.crop((left - 60, top - 60, right + 60, bottom + 60)).save(overlay)
inside = bool(np.all((result.row > top) & (result.row < bottom)
                     & (result.col > left) & (result.col < right)))
record(2, "the digitized dots over the scan, for the author to view", inside and overlay.exists(),
       f"{overlay}: each dot's center as a filled disk, red 730, green 290, blue 110 (printed 150) mbar; "
       f"the blanked label boxes in orange; every center inside the frame {inside}")

failed = [r for r in results if not r[2]]
summary = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(summary)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{summary}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
