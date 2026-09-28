"""Acceptance checks for SPEC_02 v0.8 Step 6, anchoring the geoid on the equatorial radius.

A comparison step. Every check reports its measured value; the comparison table and the wind
scaling runs are reported in full beside the expected values of SPEC_02 v0.8, and a value outside
an expected range is a finding, not a failure. The author decides the default rule from the table.
"""

import dataclasses
import math
import re
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

import numpy as np

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.control import ControlFileError
from casspian.refrac.anchor import freeze_anchor, geoid_setup
from casspian.refrac.reduce import reduce_profile

HERE = Path("reports/step02_6")
HERE.mkdir(parents=True, exist_ok=True)
PROFILE = Path("occul_data/lindal")
BUILD = PROFILE / "lindal_build.toml"
MANIFEST = PROFILE / "lindal_reduction.toml"
PRODUCT = PROFILE / "lindal_refractivity.nc"
STEP7_ASYMMETRY_M = 28745.0447
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}", flush=True)
    for line in str(detail).splitlines():
        print(f"        {line}", flush=True)


def run(*args):
    done = subprocess.run(list(args), capture_output=True, text=True)
    return done.returncode, done.stdout.strip(), done.stderr.strip()


def recorded_hashes():
    # The registered hash of each file is the last one recorded. SPEC_01 v0.19 (REPORT_01_step8),
    # SPEC_02 Step 6 (REPORT_02_step6) and SPEC_03 Step 0 (REPORT_03_step0, the whole chain
    # rebuilt) supersede REPORT_01_step9 for the files they name; a later row overrides an earlier.
    recorded = {}
    # SPEC_04 Step 0 rebuilt the whole chain once more, for kind W in three parts and kind C
    # as one structure; REPORT_04_step0 section 6 records the swept hashes.
    for name in ("reports/REPORT_01_step9.md", "reports/REPORT_01_step8.md",
                 "reports/REPORT_02_step6.md", "reports/REPORT_03_step0.md",
                 "reports/REPORT_03_step3.md", "reports/REPORT_04_step0.md"):
        if Path(name).exists():
            report = Path(name).read_text(encoding="utf-8")
            recorded.update(re.findall(
                r"\| `([\w/]+\.(?:nc|toml))` \| \w+ \| `([0-9a-f]{64})` \|", report))
    return recorded


# ---------------------------------------------------------------------------
# 1. casspian-lindal-inputs writes the new manifest and keeps kinds T and D
# ---------------------------------------------------------------------------
before = recorded_hashes()
bad_control = PROFILE / "lindal_build_mismatch_acceptance.toml"
text = BUILD.read_text(encoding="utf-8")
bad_control.write_bytes(text.replace('geoid_anchor_quantity     = "radius_equatorial_m"',
                                     'geoid_anchor_quantity     = "radius_polar_m"').encode("utf-8"))
try:
    code_bad, out_bad, err_bad = run("casspian-lindal-inputs", str(bad_control))
finally:
    bad_control.unlink()
code, out, err = run("casspian-lindal-inputs", str(BUILD))
with open(MANIFEST, "rb") as handle:
    doc = tomllib.load(handle)
kept = {name: cio.sha256(PROFILE / name) == before[name] for name in ("lindal_thermo.nc",
                                                                       "lindal_geodesy.nc")}
commits = {}
for name, kind in (("lindal_thermo.nc", "thermo"), ("lindal_geodesy.nc", "geodesy")):
    handle = cio.read(PROFILE / name, kind)
    commits[name] = str(handle.attrs["casspian_git_commit"])
    handle.close()
record(1, "casspian-lindal-inputs writes [geoid] equatorial_radius on radius_equatorial_m and "
       "[diagnostics], keeps kinds T and D, and refuses a mismatched rule and quantity",
       code == 0 and doc["geoid"]["anchor_rule"] == "equatorial_radius"
       and doc["geoid"]["anchor_quantity"] == "radius_equatorial_m"
       and doc.get("diagnostics") == {"figures": True, "format": "png", "dpi": 150}
       and all(kept.values()) and code_bad != 0 and "do not belong together" in err_bad,
       f"exit {code}\n{out}\n"
       f"[geoid] {doc['geoid']}\n[diagnostics] {doc.get('diagnostics')}\n"
       f"manifest SHA-256 {cio.sha256(MANIFEST)} (was {before['lindal_reduction.toml'][:16]}...)\n"
       + "\n".join(f"  {n}: hash as REPORT_01_step9 records {k}, commit {commits[n][:12]}"
                   for n, k in kept.items())
       + f"\nmismatched control copy: exit {code_bad}\n  {err_bad.splitlines()[-1] if err_bad else ''}")

# ---------------------------------------------------------------------------
# 2. lib.control accepts the pair and refuses the mismatches
# ---------------------------------------------------------------------------
case = HERE / "case"
if case.exists():
    shutil.rmtree(case)
case.mkdir(parents=True)
manifest = ctl.read_reduction_manifest(MANIFEST)
shutil.copy2(MANIFEST, case / MANIFEST.name)
for path in manifest.inputs.values():
    shutil.copy2(path, case / path.name)
outcomes = []
base_text = (case / MANIFEST.name).read_text(encoding="utf-8")
for rule, quantity in (("equatorial_radius", "radius_equatorial_m"),
                       ("equatorial_radius", "radius_polar_m"),
                       ("mean_polar_radius", "radius_equatorial_m"),
                       ("mean_polar_radius", "radius_polar_m")):
    edited = re.sub(r'^(anchor_rule\s*=\s*)"[^"]+"', rf'\g<1>"{rule}"', base_text, flags=re.M)
    edited = re.sub(r'^(anchor_quantity\s*=\s*)"[^"]+"', rf'\g<1>"{quantity}"', edited, flags=re.M)
    (case / MANIFEST.name).write_bytes(edited.encode("utf-8"))
    try:
        ctl.read_reduction_manifest(case / MANIFEST.name)
        outcomes.append((rule, quantity, "accepted", ""))
    except ControlFileError as exc:
        outcomes.append((rule, quantity, "refused", str(exc).split(": ", 1)[-1]))
expected = ["accepted", "refused", "refused", "accepted"]
record(2, "lib.control accepts equatorial_radius with radius_equatorial_m and refuses either rule "
       "with the other's quantity",
       [o[2] for o in outcomes] == expected,
       "\n".join(f"  {r:18s} {q:20s} {v}{('  ' + m) if m else ''}" for r, q, v, m in outcomes))

# ---------------------------------------------------------------------------
# 3. The equator is a march node and the anchoring residual is below 1e-3 m
# ---------------------------------------------------------------------------
inputs = ctl.load_reduction_inputs(manifest)
setup = geoid_setup(inputs, manifest)
at_equator = setup.march(np.array([0.0]))
record(3, "under equatorial_radius the anchor is read at the exact equator node, residual below "
       "1e-3 m",
       at_equator.anchor_node_latitude_rad == 0.0 and abs(at_equator.anchor_residual_m) < 1e-3
       and at_equator.radius[0] == at_equator.equator_radius_m,
       f"anchor node latitude {at_equator.anchor_node_latitude_rad!r} rad (exactly 0.0: "
       f"{at_equator.anchor_node_latitude_rad == 0.0})\n"
       f"anchor radius {setup.r_anchor_m!r} m; marched radius at the node "
       f"{at_equator.equator_radius_m!r} m; residual {at_equator.anchor_residual_m:.3e} m\n"
       f"the radius asked for at 0.0 equals the node value exactly, so nothing is interpolated: "
       f"{at_equator.radius[0] == at_equator.equator_radius_m}")

# ---------------------------------------------------------------------------
# 4. The product: fixed point, r0 from the march, figures without being asked, listing
# ---------------------------------------------------------------------------
figures = PROFILE / "figures"
if figures.exists():
    shutil.rmtree(figures)
code_r, out_r, err_r = run("casspian-refrac", str(MANIFEST))
tree = cio.read(PRODUCT, "refractivity")
root = tree.to_dataset(inherit=False)
rec = dict(tree["reduction_record"].attrs)
phi_c = math.radians(float(root["latitude_planetocentric_deg"].values))
r0_file = float(root["anchor_isobar_radius_m"].values)
r0_march = float(setup.march(np.array([phi_c])).radius[0])
drawn = sorted(p.name for p in figures.glob("*.png")) if figures.exists() else []
listing = ["raw/lindal_table1.csv", "raw/lindal_scalars.toml", "raw/lindal_raw.nc",
           "raw/notes.md", "lindal_build.toml", "lindal_thermo.nc", "lindal_composition.nc",
           "lindal_geodesy.nc", "lindal_gravity.nc", "lindal_rotation.nc", "lindal_wind.nc",
           "lindal_reduction.toml", "lindal_refractivity.nc"]
missing = [n for n in listing if not (PROFILE / n).exists()]
record(4, "casspian-refrac under equatorial_radius: at most eight iterations, r0 equals the march "
       "at phi_c, F1 to F4 rendered unasked, the section 2.2 listing complete",
       code_r == 0 and int(rec["fixed_point_iteration_count"]) <= 8 and r0_file == r0_march
       and len([d for d in drawn if "_F" in d]) == 4 and not missing
       and str(rec["anchor_rule"]) == "equatorial_radius",
       f"exit {code_r}{(': ' + err_r) if err_r else ''}\n{out_r}\n"
       f"anchor_rule {rec['anchor_rule']}; iterations {int(rec['fixed_point_iteration_count'])}; "
       f"iterates {np.round(np.asarray(rec['fixed_point_iterates_deg']), 7).tolist()}\n"
       f"anchor_isobar_radius_m {r0_file!r}; an independent march at phi_c {r0_march!r}; equal "
       f"{r0_file == r0_march}\nfigures rendered: {drawn}\nlisting missing: {missing or 'none'}\n"
       f"product SHA-256 {cio.sha256(PRODUCT)}, commit {tree.attrs['casspian_git_commit'][:12]}")
tree.close()

# ---------------------------------------------------------------------------
# 5. The asymmetry against Step 7
# ---------------------------------------------------------------------------
eq_anchor = freeze_anchor(inputs, manifest)
mp_manifest = dataclasses.replace(manifest, anchor_rule="mean_polar_radius",
                                  anchor_quantity="radius_polar_m")
mp_anchor = freeze_anchor(inputs, mp_manifest)
mean_eq = 0.5 * (eq_anchor.polar_north_m + eq_anchor.polar_south_m)
mean_mp = 0.5 * (mp_anchor.polar_north_m + mp_anchor.polar_south_m)
scaling = (eq_anchor.polar_asymmetry_m / mp_anchor.polar_asymmetry_m - 1.0) / (mean_eq / mean_mp - 1.0)
# SPEC_02 v0.9 Step 6: the asymmetry scales as r cubed with the size of the surface, because the
# Eq. B3 slope carries accelerations in r over a gravity in 1/r^2. The expected value is the
# Step 7 asymmetry scaled by the cube of the ratio of the mean polar radii.
expected_asymmetry = STEP7_ASYMMETRY_M * (mean_eq / mean_mp) ** 3
record(5, "the polar asymmetry under equatorial_radius equals the Step 7 value scaled by the cube "
       "of the ratio of mean polar radii, to 1 m (SPEC_02 v0.9)",
       abs(eq_anchor.polar_asymmetry_m - expected_asymmetry) <= 1.0,
       f"expected {STEP7_ASYMMETRY_M:.3f} m x ({mean_eq / 1e3:.3f} / {mean_mp / 1e3:.3f})^3 = "
       f"{expected_asymmetry:.3f} m (the review quotes 28,740.31)\n"
       f"measured under equatorial_radius {eq_anchor.polar_asymmetry_m:.3f} m; difference "
       f"{eq_anchor.polar_asymmetry_m - expected_asymmetry:+.3f} m\n"
       f"unscaled, against Step 7 ({mp_anchor.polar_asymmetry_m:.3f} m under mean_polar_radius): "
       f"{eq_anchor.polar_asymmetry_m - mp_anchor.polar_asymmetry_m:+.3f} m\n"
       f"relative change of the mean polar radius {mean_eq / mean_mp - 1.0:+.3e}; of the asymmetry "
       f"{eq_anchor.polar_asymmetry_m / mp_anchor.polar_asymmetry_m - 1.0:+.3e}; ratio {scaling:.3f}")

# ---------------------------------------------------------------------------
# 6. The comparison table: the four rules
# ---------------------------------------------------------------------------
rows = {}
for rule, quantity in ctl.ANCHOR_QUANTITY_FOR_RULE.items():
    m_rule = dataclasses.replace(manifest, anchor_rule=rule, anchor_quantity=quantity)
    a_rule = {"equatorial_radius": eq_anchor, "mean_polar_radius": mp_anchor}.get(rule) \
        or freeze_anchor(inputs, m_rule)
    row = {"r0_km": a_rule.r0_m / 1e3, "phi_c_deg": math.degrees(a_rule.phi_c_rad),
           "north_km": a_rule.polar_north_m / 1e3, "south_km": a_rule.polar_south_m / 1e3,
           "mean_polar_km": 0.5 * (a_rule.polar_north_m + a_rule.polar_south_m) / 1e3,
           "asymmetry_km": a_rule.polar_asymmetry_m / 1e3,
           "equatorial_km": a_rule.equator_radius_m / 1e3,
           "iterations": a_rule.iteration_count, "residual_m": a_rule.anchor_residual_m}
    if rule in ("equatorial_radius", "mean_polar_radius"):
        red = reduce_profile(inputs, m_rule, a_rule)
        row.update({
            "dr0_dr_anchor": red.partials["dr0_dr_anchor"],
            "anchor_term_km": red.terms["r0_anchor_radius_m"] / 1e3,
            "label_term_km": red.terms["r0_label_latitude_m"] / 1e3,
            "radius_uncertainty_km": float(red.companions["radius_uncertainty_m"].value[0]) / 1e3,
            "phi_c_uncertainty_deg": math.degrees(
                float(red.companions["latitude_planetocentric_uncertainty_rad"].value)),
        })
    rows[rule] = row

ranges = {
    "r0_km": (58516.0, 58516.5), "north_km": (54420.5, 54420.7), "south_km": (54449.2, 54449.4),
    "mean_polar_km": (54434.7, 54435.0), "asymmetry_km": (28.744, 28.746),
    "equatorial_km": (60366.999, 60367.001), "dr0_dr_anchor": (0.85, 0.97),
    "anchor_term_km": (3.4, 3.9), "label_term_km": (10.72, 10.74),
    "radius_uncertainty_km": (11.2, 11.4), "phi_c_uncertainty_deg": (0.1085, 0.1095),
}
eq = rows["equatorial_radius"]
lines = [f"{'quantity':24s}" + "".join(f"{r:>20s}" for r in rows)]
for key in ("r0_km", "phi_c_deg", "north_km", "south_km", "mean_polar_km", "asymmetry_km",
            "equatorial_km", "iterations", "residual_m", "dr0_dr_anchor", "anchor_term_km",
            "label_term_km", "radius_uncertainty_km", "phi_c_uncertainty_deg"):
    cells = []
    for rule in rows:
        value = rows[rule].get(key)
        cells.append(f"{'':>20s}" if value is None else f"{value:>20.6f}" if isinstance(value, float)
                     else f"{value:>20}")
    lines.append(f"{key:24s}" + "".join(cells))
findings = []
for key, (lo, hi) in ranges.items():
    value = eq.get(key)
    inside = value is not None and lo <= value <= hi
    findings.append(f"  {key:24s} measured {value:.6f}, expected {lo} to {hi}: "
                    f"{'inside' if inside else 'OUTSIDE (a finding)'}")
phi_shift = eq["phi_c_deg"] - rows["mean_polar_radius"]["phi_c_deg"]
findings.append(f"  phi_c shift from mean_polar_radius {phi_shift:+.6f} deg, expected under 0.002: "
                f"{'inside' if abs(phi_shift) < 0.002 else 'OUTSIDE (a finding)'}")
record(6, "the comparison table of the four rules, reported completely with the expected values",
       all(np.isfinite(v) for r in rows.values() for v in r.values()),
       "\n".join(lines) + "\nequatorial_radius against the SPEC_02 v0.8 expected values:\n"
       + "\n".join(findings))

# ---------------------------------------------------------------------------
# 7. The wind scaling diagnostic, report only
# ---------------------------------------------------------------------------
def scaled_inputs(scale):
    wind = inputs.wind.copy(deep=True)
    # SPEC_04 Step 0 deliverable 1: kind W is in three parts, u_total = u_reference + u_shear,
    # and u_cylindrical_ms is retired. Scaling all three keeps the sum identity, as scaling the
    # old three did.
    for name in ("u_total_ms", "u_reference_ms", "u_shear_ms"):
        wind[name].values[:] = wind[name].values * scale
    return dataclasses.replace(inputs, wind=wind)


scaling_rows, lines = [], []
for rule, m_rule, base in (("equatorial_radius", manifest, eq_anchor),
                           ("mean_polar_radius", mp_manifest, mp_anchor)):
    for scale in (0.95, 1.05):
        a_s = freeze_anchor(scaled_inputs(scale), m_rule)
        scaling_rows.append((rule, scale, a_s, base))
        lines.append(
            f"  {rule:18s} u x {scale:4.2f}: r0 {a_s.r0_m / 1e3:.4f} km ({(a_s.r0_m - base.r0_m) / 1e3:+.4f}), "
            f"phi_c {math.degrees(a_s.phi_c_rad):.6f} ({math.degrees(a_s.phi_c_rad - base.phi_c_rad):+.6f}), "
            f"north {a_s.polar_north_m / 1e3:.4f}, south {a_s.polar_south_m / 1e3:.4f}, "
            f"equatorial {a_s.equator_radius_m / 1e3:.4f} km")
shift = {}
for rule in ("equatorial_radius", "mean_polar_radius"):
    pair = [a for r, _, a, _ in scaling_rows if r == rule]
    shift[rule] = abs(pair[1].r0_m - pair[0].r0_m) / 2e3
lines.append(f"half range of r0 over u x 0.95 to 1.05: equatorial_radius {shift['equatorial_radius']:.4f} km "
             f"(expected a few km), mean_polar_radius {shift['mean_polar_radius']:.4f} km (expected under 1 km)")
record(7, "the wind scaling diagnostic, four runs, reported completely", len(scaling_rows) == 4,
       "\n".join(lines))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
