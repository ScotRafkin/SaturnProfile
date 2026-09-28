"""Acceptance checks for SPEC_02 v0.7 Step 5, the standard diagnostics.

Every check reports its measured value. `casspian-refrac` with `figures = true` runs on a copy of
the Lindal manifest and its six inputs under this directory: the committed manifest carries no
`[diagnostics]` section, and adding one would change its hash, which the kind N product and the
earlier reports record. Nothing under `occul_data/` is written except by check 7, which draws
the committed product into its default, ignored `figures/` directory.
"""

import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
from matplotlib.image import imread

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib import reduction as red
from casspian.lib.gravity import g_eff_radial
from casspian.tools.plots import render, style

HERE = Path("reports/step02_5")
HERE.mkdir(parents=True, exist_ok=True)
PROFILE = Path("occul_data/lindal")
MANIFEST = PROFILE / "lindal_reduction.toml"
FIGURES_REPORT = Path("reports/figures")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def run(*args):
    done = subprocess.run(list(args), capture_output=True, text=True)
    return done.returncode, done.stdout.strip(), done.stderr.strip()


# ---------------------------------------------------------------------------
# 1. casspian-refrac with figures = true
# ---------------------------------------------------------------------------
case = HERE / "case"
if case.exists():
    shutil.rmtree(case)
case.mkdir(parents=True)
shutil.copy2(MANIFEST, case / MANIFEST.name)
for path in ctl.read_reduction_manifest(MANIFEST).inputs.values():
    shutil.copy2(path, case / path.name)
manifest_copy = case / MANIFEST.name
# From SPEC_02 v0.8 Step 6 the committed manifest carries [diagnostics] itself; append the
# section only to a manifest that lacks it.
if b"[diagnostics]" not in manifest_copy.read_bytes():
    manifest_copy.write_bytes(manifest_copy.read_bytes()
                              + b'\n[diagnostics]\nfigures = true\nformat  = "png"\ndpi     = 150\n')
code, out, err = run("casspian-refrac", str(manifest_copy))
figures_dir = case / "figures"
expected = [f"lindal_diag_F{k}_{n}.png" for k, n in
            ((1, "inputs"), (2, "gravity_profile"), (3, "gravity_latitude"), (4, "product"))]
present = {name: (figures_dir / name).exists() for name in expected}
pdf = figures_dir / "lindal_diag.pdf"
pages = len(re.findall(rb"/Type\s*/Page[^s]", pdf.read_bytes())) if pdf.exists() else 0
skips = [line for line in out.splitlines() if line.startswith("skipped")]
# SPEC_03 Step 3 retires the "F5 and F6 skipped" check: F5 and F6 belong to kind profile, and kind
# N renders F1 to F4 and reports nothing skipped (SPEC_02 v0.10 decision 10 superseded).
record(1, "casspian-refrac with figures = true writes F1 to F4 and the combined PDF, and reports "
       "nothing skipped (SPEC_03 Step 3)",
       code == 0 and all(present.values()) and pages == 4 and not skips,
       f"exit {code}{(': ' + err) if err else ''}\n{out}\n"
       + "\n".join(f"  {n:34s} {'present' if ok else 'MISSING'}" for n, ok in present.items())
       + f"\n  lindal_diag.pdf pages: {pages}")

# ---------------------------------------------------------------------------
# 2. casspian-plots by hand writes identical figures, apart from the generation time
# ---------------------------------------------------------------------------
by_hand = HERE / "by_hand"
if by_hand.exists():
    shutil.rmtree(by_hand)
code2, out2, err2 = run("casspian-plots", str(case / "lindal_refractivity.nc"), "--out", str(by_hand))
lines, identical = [], code2 == 0
for name in expected:
    a, b = imread(figures_dir / name), imread(by_hand / name)
    same_shape = a.shape == b.shape
    band = math.ceil(style.FOOTER_TIME_BAND * a.shape[0])
    body_same = same_shape and bool(np.array_equal(a[:-band], b[:-band]))
    band_differs = same_shape and not np.array_equal(a[-band:], b[-band:])
    bytes_same = (figures_dir / name).read_bytes() == (by_hand / name).read_bytes()
    identical = identical and body_same
    lines.append(f"  {name:34s} {a.shape[1]}x{a.shape[0]} px; identical above the "
                 f"{band} px time band: {body_same}; time band differs: {band_differs}; "
                 f"whole file byte identical: {bytes_same}")
record(2, "casspian-plots on the same file by hand writes identical figures, the generation "
       "time masked", identical,
       f"exit {code2}{(': ' + err2) if err2 else ''}\n" + "\n".join(lines)
       + "\nThe PNGs are compared as images with the bottom band, which holds only the "
         "generation time, masked. A time identical to the second gives byte identical files.")

# ---------------------------------------------------------------------------
# 3. Every figure carries the footer
# ---------------------------------------------------------------------------
product = case / "lindal_refractivity.nc"
result = render(product, HERE / "in_process")
tree = cio.read(product, "refractivity")
commit = str(tree.attrs["casspian_git_commit"])
sha12 = cio.sha256(product)[:12]
footer_ok = all(product.name in f and commit in f and sha12 in f for f in result.footers.values())
inked = []
for path in result.written:
    image = imread(path)[..., :3]
    height = image.shape[0]
    top = height - math.ceil(style.FOOTER_TOP * height)
    inked.append((path.name, bool(np.any(image[top:] < 0.9)),
                  bool(np.any(image[-math.ceil(style.FOOTER_TIME_BAND * height):] < 0.9))))
record(3, "every figure carries the footer: file name, casspian_git_commit, SHA-256 prefix and "
       "the generation time",
       footer_ok and all(a and b for _, a, b in inked) and bool(result.generated_at),
       "\n".join(f"  {key}: {text}" for key, text in result.footers.items())
       + f"\n  generated {result.generated_at}\n"
       + "\n".join(f"  {n}: footer ink {a}, time band ink {b}" for n, a, b in inked))

# ---------------------------------------------------------------------------
# 4. The wind view is written and the gravity file is refused by kind
# ---------------------------------------------------------------------------
wind_out = HERE / "wind_view"
if wind_out.exists():
    shutil.rmtree(wind_out)
code_w, out_w, err_w = run("casspian-plots", str(PROFILE / "lindal_wind.nc"), "--out", str(wind_out))
code_g, out_g, err_g = run("casspian-plots", str(PROFILE / "lindal_gravity.nc"), "--out",
                           str(HERE / "gravity_view"))
wind_files = sorted(p.name for p in wind_out.glob("*")) if wind_out.exists() else []
record(4, "casspian-plots lindal_wind.nc writes the single figure wind view, and "
       "lindal_gravity.nc is refused by kind",
       code_w == 0 and wind_files == ["lindal_diag_wind.png"] and code_g != 0
       and "gravity" in err_g and not (HERE / "gravity_view").exists(),
       f"wind: exit {code_w}, {out_w}; files {wind_files}\n"
       f"gravity: exit {code_g}, {err_g}")

# ---------------------------------------------------------------------------
# 5. F4's recovered temperature is within 1e-12 at every level
# ---------------------------------------------------------------------------
fraction = np.asarray(result.data["F4"]["recovered_temperature_fraction"])
root = tree.to_dataset(inherit=False)
thermo = tree["inputs/thermo"].to_dataset(inherit=False)
independent = red.temperature_from_refractivity(
    thermo["pressure_Pa"].values, root["mean_refractivity_m3"].values,
    root["refractivity"].values) / thermo["temperature_K"].values - 1.0
record(5, "F4's recovered temperature panel shows a fractional difference below 1e-12 at every "
       "level",
       fraction.size == 66 and float(np.max(np.abs(fraction))) < 1e-12
       and np.array_equal(fraction, independent),
       f"{fraction.size} levels; max |T_recovered / T - 1| = {float(np.max(np.abs(fraction))):.2e}; "
       f"the panel's array equals an independent evaluation: {np.array_equal(fraction, independent)}")

# ---------------------------------------------------------------------------
# 6. F2's g_eff at 1 bar agrees with lib.gravity called directly
# ---------------------------------------------------------------------------
manifest = ctl.read_reduction_manifest(MANIFEST)
inputs = ctl.load_reduction_inputs(manifest)
k = int(np.flatnonzero(thermo["pressure_Pa"].values == 1.0e5)[0])
F2 = result.data["F2"]
gravity, rotation, wind = inputs.gravity, inputs.rotation, inputs.wind
lat = wind["latitude_planetocentric_deg"].values
column = int(np.flatnonzero(wind["pressure_Pa"].values == 1.0e5)[0])
phi_c = float(root["latitude_planetocentric_deg"].values)
u_direct = float(np.interp(phi_c, lat, wind["u_total_ms"].values[:, column]))
g_direct = float(g_eff_radial(u_direct, float(root["radius_m"].values[k]), math.radians(phi_c),
                              float(rotation["angular_rate_rad_s"]), float(gravity["GM_m3s2"]),
                              gravity["J"].values, gravity["degree"].values,
                              float(gravity["normalization_radius_m"])))
g_panel = float(F2["g_eff_ms2"][k])
record(6, "F2's g_eff at the 1 bar level agrees with lib.gravity called directly to round-off",
       abs(g_panel - g_direct) <= 4 * np.finfo(float).eps * abs(g_direct),
       f"panel {g_panel!r} m/s2, direct {g_direct!r} m/s2 (from the files on disk, not the "
       f"embedded copies), difference {g_panel - g_direct:.1e}\n"
       f"u at phi_c on the 1 bar level: panel {float(F2['u_ms'][k])!r}, direct {u_direct!r}")

# ---------------------------------------------------------------------------
# 7. Beyond the specification: default location, the T and C views, and F5, F6 when present
# ---------------------------------------------------------------------------
code_d, out_d, err_d = run("casspian-plots", str(PROFILE / "lindal_refractivity.nc"))
default_dir = PROFILE / "figures"
ignored = run("git", "check-ignore", str(default_dir / "lindal_diag.pdf"))[0] == 0
code_t, out_t, err_t = run("casspian-plots", str(PROFILE / "lindal_thermo.nc"), "--out", str(HERE / "views"))
code_c, out_c, err_c = run("casspian-plots", str(PROFILE / "lindal_composition.nc"), "--out", str(HERE / "views"))

# SPEC_03 Step 3 retires the synthetic F5 and F6 on kind N: those figures belong to kind profile
# and are tested in the Step 03_3 acceptance. Kind N renders F1 to F4 and reports nothing skipped.
full = render(PROFILE / "lindal_refractivity.nc", HERE / "product_figures")
tree.close()
record(7, "beyond the specification: the default figures/ directory is ignored by git, the T and "
       "C views render, and kind N renders F1 to F4 with nothing skipped",
       code_d == 0 and ignored and code_t == 0 and code_c == 0
       and [p.name for p in full.written] == expected and not full.skipped,
       f"default: exit {code_d}; {len(out_d.splitlines())} lines; figures/ ignored by git: {ignored}\n"
       f"thermo view: exit {code_t}, {out_t}\ncomposition view: exit {code_c}, {out_c}\n"
       f"kind N rendered: written {[p.name for p in full.written]}, skipped {full.skipped}")

# ---------------------------------------------------------------------------
# Figures attached to the report
# ---------------------------------------------------------------------------
FIGURES_REPORT.mkdir(parents=True, exist_ok=True)
attached = []
# F1 to F4 are attached from the committed product's own figures/ directory (check 7), not from
# the case copy, whose product was written from the working tree and carries -dirty.
for source in [default_dir / n for n in expected] + [wind_out / "lindal_diag_wind.png",
                                                     HERE / "views" / "lindal_diag_thermo.png",
                                                     HERE / "views" / "lindal_diag_composition.png"]:
    target = FIGURES_REPORT / f"step02_5_{source.name}"
    shutil.copy2(source, target)
    attached.append(target.as_posix())
print("attached to the report: " + ", ".join(attached))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
