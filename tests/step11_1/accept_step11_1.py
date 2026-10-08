"""Acceptance for SPEC_11 v0.3 Step 1: the Lindal-only wind, four cases, for the author.

Builds the four winds through the `lindal_iris` case of `casspian-wind-shear` (SPEC_11 v0.4 Step 2;
the tool's own `build` was removed then) from control files written under `reports/step11_1/`, on the `lindal_transfer` run's inputs: its cloud wind (the registered
`lindal_transfer_wind_source.nc`), composition, gravity and rotation, and the committed IRIS
temperatures. Nothing committed is written. The figures are the author's (SPEC_11 Step 1
acceptance 1 to 3); the checks below are what every wind the model runs on must satisfy.

Run from the repository root.
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from casspian.lib import io as cio
from casspian.lib.control import _admit_wind
from casspian.tools.lindal import lindal_wind as lw
from casspian.tools.wind import shear

HERE = Path("reports/step11_1")
HERE.mkdir(parents=True, exist_ok=True)
RUN = Path("forward/lindal_transfer/inputs")
INPUTS = {"cloud_wind": RUN / "lindal_transfer_wind_source.nc",
          "temperatures": Path("occul_data/lindal/iris_temperatures.csv"),
          "composition": RUN / "lindal_transfer_composition.nc",
          "gravity": RUN / "lindal_transfer_gravity.nc",
          "rotation": RUN / "lindal_transfer_rotation.nc"}
CASES = {"i_110": (11000.0, "held"), "ii_110": (11000.0, "relaxed"),
         "i_150": (15000.0, "held"), "ii_150": (15000.0, "relaxed")}
PDF = Path("docs/conrath_etal_IRIS_1983.pdf")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def load(path, kind):
    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def control(name, top, above):
    path = HERE / name / f"{name}.toml"
    path.parent.mkdir(parents=True, exist_ok=True)
    rel = {("source" if k == "cloud_wind" else k): Path("../../..") / v for k, v in INPUTS.items()}
    lines = ["[shear]"] + [f'{k} = "{v.as_posix()}"' for k, v in rel.items()]
    lines += [f'output = "{name}_wind.nc"', 'case = "lindal_iris"', f"top_level_Pa = {top!r}",
              f'above_top = "{above}"']
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


# ---------------------------------------------------------------------------
# The four winds, and the constructions behind them for the figures
# ---------------------------------------------------------------------------
cloud, comp = load(INPUTS["cloud_wind"], "wind"), load(INPUTS["composition"], "composition")
grav, rot = load(INPUTS["gravity"], "gravity"), load(INPUTS["rotation"], "rotation")
table = np.genfromtxt(INPUTS["temperatures"], delimiter=",", names=True)
latitude = np.asarray(cloud["latitude_planetocentric_deg"].values, dtype="float64")
latitude_g = np.asarray(cloud["latitude_planetographic_deg"].values, dtype="float64")

written, built = {}, {}
for name, (top, above) in CASES.items():
    written[name] = shear.build(control(name, top, above))
    built[name] = lw.construct(cloud, table, comp, grav, rot, top, above)

# ---------------------------------------------------------------------------
# 0. Each wind is one the model can run on, and equals its construction
# ---------------------------------------------------------------------------
lines, ok = [], True
for name, path in written.items():
    wind = load(path, "wind")
    pressure = np.asarray(wind["pressure_Pa"].values)
    message = None
    try:
        _admit_wind(wind, path.name)
    except Exception as exc:
        message = str(exc)
    same = np.array_equal(wind["u_total_ms"].values, built[name].dataset["u_total_ms"].values)
    reference = bool(np.any(pressure == float(wind["reference_level_pressure_Pa"])))
    parts = np.max(np.abs(wind["u_total_ms"].values - wind["u_reference_ms"].values[:, None]
                          - wind["u_shear_ms"].values))
    at_reference = float(np.max(np.abs(wind["u_shear_ms"].values[:, pressure == lw.REFERENCE_PRESSURE_Pa])))
    good = message is None and same and reference and at_reference == 0.0
    ok &= good
    lines.append(f"{name}: {pressure.size} pressures {pressure.min():g} to {pressure.max():g} Pa, "
                 f"reference {float(wind['reference_level_pressure_Pa']):g} Pa a node {reference}; "
                 f"admitted by the loaders' check (u_total, poles exactly zero) {message is None}; "
                 f"array-equal to the construction {same}; shear at the reference level "
                 f"{at_reference:.1e} m/s; |u_total - u_reference - u_shear| {parts:.1e} m/s; "
                 f"{built[name].iterations} iterations, last change {built[name].changes_ms[-1]:.3f} m/s "
                 f"(tolerance {lw.TOLERANCE_MS} m/s); u_total {float(wind['u_total_ms'].min()):.0f} to "
                 f"{float(wind['u_total_ms'].max()):.0f} m/s"
                 + ("" if message is None else f"\n    {message}"))
record(0, "each wind reads as kind W, passes the loaders' check, is array-equal to its construction, "
          "and carries zero shear at its reference level", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 1. The shear at the three levels
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
lines = []
for k, ax in enumerate(axes):
    for name, color in (("i_110", "C0"), ("i_150", "C3")):
        c = built[name]
        ax.plot(latitude, c.shear_ms[k], color=color, linewidth=1.0,
                label=f"top level at {c.level_Pa[0] / 100:g} mbar: shaped")
        ax.plot(latitude, c.shear_computed_ms[k], ".", color=color, markersize=1.5,
                label=f"top level at {c.level_Pa[0] / 100:g} mbar: computed")
    level = built["i_110"].level_Pa[k] / 100 if k else "110 or 150"
    ax.set_ylabel("du/dln p (m/s per H)")
    ax.set_title(f"the shear at the {'top' if k == 0 else f'{level:g} mbar'} level", fontsize=10)
    ax.axhline(0, color="0.6", linewidth=0.8)
    ax.axvspan(-lw.EQUATORIAL_BAND_DEG, lw.EQUATORIAL_BAND_DEG, color="0.9")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, loc="best")
    for name in ("i_110", "i_150"):
        s = built[name].shear_computed_ms[k]
        lines.append(f"{name}, level {k + 1} ({built[name].level_Pa[k] / 100:g} mbar): computed "
                     f"{np.nanmin(s):.1f} to {np.nanmax(s):.1f} m/s per H, at "
                     f"{latitude[np.nanargmin(s)]:.1f} and {latitude[np.nanargmax(s)]:.1f} deg")
axes[-1].set_xlabel("planetocentric latitude (deg)")
axes[-1].set_xlim(90, -90)
fig.tight_layout()
fig.savefig(HERE / "shear_levels.png", dpi=130)
plt.close(fig)
apart = max(float(np.max(np.abs(built[a].shear_ms - built[b].shear_ms)))
            for a, b in (("i_110", "ii_110"), ("i_150", "ii_150")))
record(1, "the shear at the three levels against latitude, for the author", True,
       "\n".join(lines) + f"\nthe held and relaxed cases' level shears differ by at most "
       f"{apart:.3f} m/s per H: the shape above the top level does not enter the balance at the "
       f"levels, and what remains is each iteration stopping within its {lw.TOLERANCE_MS} m/s "
       f"tolerance\nfigure {HERE / 'shear_levels.png'}")

# ---------------------------------------------------------------------------
# 2. The 150 mbar shear drawn as Conrath and Pirraglia Fig. 3
# ---------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 6.5), gridspec_kw={"width_ratios": [1, 1.3]})
c = built["i_150"]
axes[0].plot(c.shear_ms[0], latitude_g, color="k", linewidth=1.0, label="shaped")
axes[0].plot(c.shear_computed_ms[0], latitude_g, ".", color="C3", markersize=1.5, label="computed")
axes[0].axvline(0, color="0.5", linewidth=0.8)
axes[0].set_xlim(-40, 50)
axes[0].set_ylim(-90, 90)
axes[0].set_xlabel("-du/dz = du/dln p (m/s per H)")
axes[0].set_ylabel("planetographic latitude (deg)")
axes[0].set_title("this tool, top level at 150 mbar", fontsize=10)
axes[0].legend(fontsize=8)
axes[0].grid(alpha=0.3)
scan_note = "the scan of Fig. 3 is drawn beside it"
if PDF.exists():
    with tempfile.TemporaryDirectory() as directory:
        page = Path(directory) / "p4.png"
        from casspian.tools.lindal.iris_temperatures import find_ghostscript
        subprocess.run([find_ghostscript(), "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                        "-sDEVICE=pnggray", "-r300", "-dFirstPage=4", "-dLastPage=4",
                        f"-sOutputFile={page}", str(PDF)], check=True)
        with Image.open(page) as image:
            s = 300 / 72
            axes[1].imshow(np.asarray(image.crop((int(255 * s), int(435 * s), int(445 * s),
                                                   int(590 * s)))), cmap="gray")
else:
    scan_note = f"{PDF} is absent, so the scan is not drawn"
axes[1].axis("off")
axes[1].set_title("Conrath and Pirraglia (1983) Fig. 3", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "shear_150mbar_as_fig3.png", dpi=130)
plt.close(fig)
record(2, "the 150 mbar shear drawn as Conrath and Pirraglia's Fig. 3, for the author", True,
       f"-du/dz in m/s per scale height (z in scale heights is -ln p, so -du/dz is du/dln p) against "
       f"planetographic latitude, the figure's axes; {scan_note}\n"
       f"figure {HERE / 'shear_150mbar_as_fig3.png'}")

# ---------------------------------------------------------------------------
# 3. Each wind
# ---------------------------------------------------------------------------
profiles = (60.0, 30.0, 0.0, -30.0, -60.0)
for name, path in written.items():
    wind = load(path, "wind")
    p = np.asarray(wind["pressure_Pa"].values)
    u = np.asarray(wind["u_total_ms"].values)
    fig, (left, right) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.6, 1]})
    filled = left.contourf(latitude, p / 100.0, u.T, levels=np.arange(-600, 701, 50), cmap="RdBu_r",
                           extend="both")
    fig.colorbar(filled, ax=left, label="u_total (m/s)")
    left.set_yscale("log")
    left.invert_yaxis()
    left.set_xlim(90, -90)
    for level in built[name].level_Pa:
        left.axhline(level / 100.0, color="k", linewidth=0.6, linestyle=":")
    left.axhline(lw.REFERENCE_PRESSURE_Pa / 100.0, color="k", linewidth=0.9)
    left.set_xlabel("planetocentric latitude (deg)")
    left.set_ylabel("pressure (mbar)")
    left.set_title(f"u_total, case {name}", fontsize=10)
    for lat0 in profiles:
        right.plot([np.interp(lat0, latitude, u[:, j]) for j in range(p.size)], p / 100.0,
                   label=f"{lat0:+.0f} deg")
    right.set_yscale("log")
    right.invert_yaxis()
    right.axhline(lw.REFERENCE_PRESSURE_Pa / 100.0, color="k", linewidth=0.9)
    right.set_xlabel("u (m/s)")
    right.set_title("u(p) at five planetocentric latitudes", fontsize=10)
    right.legend(fontsize=8)
    right.grid(alpha=0.3)
    fig.suptitle(str(wind.attrs["vertical_structure"]), fontsize=9)
    fig.tight_layout()
    fig.savefig(HERE / f"wind_{name}.png", dpi=120)
    plt.close(fig)
record(3, "each wind's u_total(phi, p) and u(p) at five latitudes, for the author", True,
       "figures " + ", ".join(str(HERE / f"wind_{n}.png") for n in written))

# ---------------------------------------------------------------------------
# 4. The geometry's anchor (beyond the specification, reported)
# ---------------------------------------------------------------------------
R = float(grav["normalization_radius_m"])
lines = []
for delta in (-300e3, 300e3):
    other = lw.construct(cloud, table, comp, grav, rot, 11000.0, "held", equatorial_radius_m=R + delta)
    change = np.nanmax(np.abs(other.shear_computed_ms - built["i_110"].shear_computed_ms), axis=1)
    lines.append(f"equatorial radius {(R + delta) / 1e3:.0f} km against {R / 1e3:.0f} km: largest "
                 f"change of the computed shear at the three levels {np.round(change, 3).tolist()} "
                 f"m/s per H")
record(4, "beyond the specification: the shear's sensitivity to the radius the isobars are anchored "
          "at, reported", True, "\n".join(lines))

summary = {name: {"iterations": c.iterations, "changes_ms": c.changes_ms,
                  "level_Pa": c.level_Pa.tolist()} for name, c in built.items()}
(HERE / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
