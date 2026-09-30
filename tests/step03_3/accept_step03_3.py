"""Acceptance checks for SPEC_03 v0.10 Step 3: the season identifier and relative paths across the
chain, the run directory and its inputs, the run namelist and its refusals, kind profile with F5
and F6, and the input_hashes warning for every derived kind.

Every check prints its measured value; checks beyond the specification are labeled so. The chain
was rebuilt in the working tree and the candidate kind N built in memory by `rebuild.sh` and
`build_candidate.py` (SPEC_03 v0.6 section 0); the run's inputs were written by
`casspian-run-inputs`. The `-dirty` refusal is relaxed in this process by `relaxed.py` and named in
the output; the one case that tests the refusal itself turns the relaxation off. Nothing under
`occul_data/` or `forward/` is written by this script.
"""

import math
import os
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import relaxed  # noqa: E402

import netCDF4  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

from casspian.lib import control as ctl  # noqa: E402
from casspian.lib import geopotential as gp  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.lib.constants import BOLTZMANN_CONSTANT  # noqa: E402
from casspian.lib.schema import CasspianSchemaError, uncertainty_companion  # noqa: E402
from casspian.refrac.anchor import wind_of_latitude  # noqa: E402
from casspian.tools.gravity import build_gravity  # noqa: E402
from casspian.tools.plots import render  # noqa: E402
from casspian.tools.run import run_inputs  # noqa: E402

D = Path("occul_data/lindal")
RUN = Path("forward/lindal_closure")
HERE = Path("reports/step03_3")
FIXTURES = Path(__file__).resolve().parent / "fixtures"
BEFORE = FIXTURES / "before"
CANDIDATE = FIXTURES / "candidate" / "lindal_refractivity.nc"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def fresh(path):
    path = Path(path)
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


REDUCTION = [("raw/lindal_raw.nc", "raw"), ("lindal_thermo.nc", "thermo"),
             ("lindal_geodesy.nc", "geodesy"), ("lindal_gravity.nc", "gravity"),
             ("lindal_rotation.nc", "rotation"), ("lindal_wind.nc", "wind"),
             ("lindal_composition.nc", "composition")]
RUN_INPUTS = [("inputs/lindal_closure_composition.nc", "composition"),
              ("inputs/lindal_closure_gravity.nc", "gravity"),
              ("inputs/lindal_closure_rotation.nc", "rotation"),
              ("inputs/lindal_closure_wind.nc", "wind")]
FILES = ([(D / n, k) for n, k in REDUCTION] + [(CANDIDATE, "refractivity")]
         + [(RUN / n, k) for n, k in RUN_INPUTS])
EXPECT_SEASON = {"raw": 18.2, "thermo": 18.2, "geodesy": 18.2, "wind": 18.2, "refractivity": 18.2,
                 "gravity": "uniform", "rotation": "uniform", "composition": "uniform"}
print(f"relaxation: -dirty refusal relaxed in this process ({relaxed.MARK})\n")

# ---------------------------------------------------------------------------
# 1. Every product and every run input validates with the season attributes
# ---------------------------------------------------------------------------
lines, ok = [], True
for path, kind in FILES:
    try:
        handle = cio.read(path, kind)
        attrs = dict(handle.attrs)
        handle.close()
    except CasspianSchemaError as exc:
        ok = False
        lines.append(f"{path}: REFUSED {exc}")
        continue
    season = attrs.get("solar_longitude_deg", attrs.get("season_absent_meaning"))
    want = EXPECT_SEASON[kind]
    good = (float(season) == want) if isinstance(want, float) else (season == want)
    # SPEC_00 v0.19, SPEC_01 v0.28 (REVIEW_03_step3 finding 1): one date for the whole chain.
    good = good and attrs.get("epoch") == "1981-08-26"
    ok = ok and good
    lines.append(f"{path.as_posix():58s} {kind:12s} epoch {attrs.get('epoch')!r:34s} season "
                 f"{season!r}{'' if good else f'  <-- expected {want!r}'}"
                 + (f"; epoch_note {attrs['epoch_note'][:60]!r}..." if "epoch_note" in attrs else ""))
record(1, "every product of the rebuilt chain, the candidate kind N and every run input validates "
          "with the season attributes, every one dated 1981-08-26 (SPEC_03 v0.12)", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. Neither and both are refused
# ---------------------------------------------------------------------------
cases = fresh(HERE / "season_cases")
neither, both = cases / "lindal_rotation_neither.nc", cases / "lindal_rotation_both.nc"
shutil.copy2(D / "lindal_rotation.nc", neither)
shutil.copy2(D / "lindal_rotation.nc", both)
with netCDF4.Dataset(neither, "a") as h:
    h.delncattr("season_absent_meaning")
with netCDF4.Dataset(both, "a") as h:
    h.setncattr("solar_longitude_deg", 18.2)
    h.setncattr("solar_longitude_source", "acceptance edit")
prose = cases / "lindal_rotation_prose_epoch.nc"
shutil.copy2(D / "lindal_rotation.nc", prose)
with netCDF4.Dataset(prose, "a") as h:
    h.setncattr("epoch", "Voyager 1 and 2, 1980 to 1981")
lines, ok = [], True
for label, path in (("neither", neither), ("both", both),
                    ("an epoch that is not an ISO date (SPEC_03 v0.11)", prose)):
    try:
        cio.read(path, "rotation").close()
        ok = False
        lines.append(f"{label}: not refused")
    except CasspianSchemaError as exc:
        lines.append(f"{label}: {exc}")
record(2, "a copy with neither solar_longitude_deg nor season_absent_meaning, a copy with both, and "
          "a copy whose epoch is not an ISO date are each refused", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. Kind N's season equals kind T's; 4. the transcribed pair satisfies sin(d) = sin(e) sin(Ls)
# ---------------------------------------------------------------------------
t = cio.read(D / "lindal_thermo.nc", "thermo")
n = cio.read(CANDIDATE, "refractivity")
keys = ("epoch", "solar_longitude_deg", "solar_longitude_source")
same = {k: (t.attrs[k], n.attrs[k]) for k in keys}
equal = all(np.all(a == b) for a, b in same.values())
record(3, "kind N's season equals kind T's", equal,
       "\n".join(f"{k}: T {a!r}\n   N {b!r}" for k, (a, b) in same.items()))
t.close()
n.close()
raw = cio.read(D / "raw/lindal_raw.nc", "raw")
ls, delta = float(raw.attrs["solar_longitude_deg"]), float(raw.attrs["subsolar_latitude_deg"])
source = raw["scalars/source"].attrs
raw.close()
EPSILON = 26.73
implied = math.degrees(math.asin(math.sin(math.radians(EPSILON)) * math.sin(math.radians(ls))))
record(4, "the transcribed pair satisfies sin(delta_s) = sin(26.73 deg) sin(Ls) to 0.05 deg",
       abs(implied - delta) <= 0.05,
       f"Ls {ls} deg, sub-solar latitude {delta} deg (raw bundle root and scalars/source: "
       f"{float(source['solar_longitude_deg'])}, {float(source['subsolar_latitude_deg'])})\n"
       f"asin(sin 26.73 sin {ls}) = {implied:.4f} deg; difference {implied - delta:+.4f} deg\n"
       f"epoch {source['observation_date']!r}; source {str(source['solar_longitude_source'])[:120]}...")

# ---------------------------------------------------------------------------
# 5. role from the build file: reduction products, forward inputs, a section without role refused
# ---------------------------------------------------------------------------
lines, ok = [], True
for path, kind in FILES:
    handle = cio.read(path, kind)
    role, crole = handle.attrs.get("role"), handle.attrs.get("composition_role")
    handle.close()
    want = "forward" if "forward" in path.as_posix() else "reduction"
    good = role == want and (crole is None or crole == want)
    ok = ok and good
    lines.append(f"{path.name:36s} role {role!r}" + (f", composition_role {crole!r}" if crole else "")
                 + ("" if good else f"  <-- expected {want!r}"))
no_role = fresh(HERE / "role_case") / "lindal_build.toml"
text = (D / "lindal_build.toml").read_text(encoding="utf-8")
assert text.count('role   = "reduction"   # SPEC_01 v0.26: required in every section\n') == 1
no_role.write_text(text.replace('role   = "reduction"   # SPEC_01 v0.26: required in every section\n', ""),
                   encoding="utf-8")
try:
    build_gravity.build(no_role, "gravity")
    ok = False
    lines.append("a [gravity] section without role: not refused")
except ctl.ControlFileError as exc:
    ok = ok and "role" in str(exc)
    lines.append(f"a [gravity] section without role: {exc}")
record(5, "beyond the specification: role written from the build file, and a section without role "
          "refused (section 8 ruling 1)", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 6. Every path-bearing attribute is relative: no drive letter, no leading slash
# ---------------------------------------------------------------------------
PATH_ATTRS = ("input_hashes", "control_file", "raw_bundle", "raw_sources", "master_table",
              "master_table_hash", "latitude_conversion_inputs", "decomposition_geometry")
ABSOLUTE = re.compile(r"^(?:[A-Za-z]:[\\/]|/)")


def scan(path):
    found, bad = 0, []

    def walk(group, where):
        nonlocal found
        for name in group.ncattrs():
            value = str(group.getncattr(name))
            if name not in PATH_ATTRS and "sha256:" not in value:
                continue
            for line in value.splitlines():
                token = line.rsplit(" sha256:", 1)[0].strip()
                if name == "decomposition_geometry" and " sha256:" not in line:
                    continue
                found += 1
                if ABSOLUTE.match(token) or re.search(r"[A-Za-z]:[\\/]", token):
                    bad.append(f"{where}{name}: {token}")
        for child_name, child in group.groups.items():
            walk(child, f"{where}{child_name}/")

    with netCDF4.Dataset(path, "r") as root:
        walk(root, "/")
    return found, bad


lines, ok = [], True
for path, _ in FILES:
    found, bad = scan(path)
    ok = ok and not bad and found > 0
    lines.append(f"{path.as_posix():58s} {found:3d} recorded paths, absolute: {bad or 'none'}")
_, bad = scan(D / "lindal_thermo.nc")
with netCDF4.Dataset(D / "lindal_thermo.nc") as h:
    lines.append("kind T input_hashes, as written:\n  " + h.getncattr("input_hashes").replace("\n", "\n  "))
with netCDF4.Dataset(RUN / "inputs/lindal_closure_wind.nc") as h:
    lines.append("run wind latitude_conversion_inputs, as written:\n  "
                 + h.getncattr("latitude_conversion_inputs").replace("\n", "\n  "))
record(6, "every product of the rebuilt chain and every run input carries only relative paths "
          "(a scan of every path-bearing attribute in every group)", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 7. The input_hashes warning, resolved against the file's own directory (deliverable 4)
# ---------------------------------------------------------------------------
beside = fresh(HERE / "warning_beside")
(beside / "raw").mkdir()
shutil.copy2(D / "raw/lindal_raw.nc", beside / "raw/lindal_raw.nc")
edited = beside / "lindal_thermo.nc"
shutil.copy2(D / "lindal_thermo.nc", edited)
with netCDF4.Dataset(edited, "a") as h:
    entries = h.getncattr("input_hashes").splitlines()
    changed = []
    for line in entries:
        if line.startswith("raw/lindal_raw.nc sha256:"):
            line = line[:-1] + ("0" if line[-1] != "0" else "1")
        changed.append(line)
    h.setncattr("input_hashes", "\n".join(changed))
alone = fresh(HERE / "warning_alone") / "lindal_thermo.nc"
shutil.copy2(edited, alone)


def warnings_on_read(path):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        cio.read(path, "thermo").close()
    return [str(w.message) for w in caught if "input_hashes" in str(w.message)]


w_beside, w_alone = warnings_on_read(edited), warnings_on_read(alone)
w_clean = warnings_on_read(D / "lindal_thermo.nc")
record(7, "a kind T copy whose recorded raw bundle hash is edited warns with the bundle beside it at "
          "the recorded relative path, and reads silently alone",
       len(w_beside) == 1 and not w_alone and not w_clean,
       f"beside ({beside.as_posix()}, raw/ holding the bundle): {w_beside}\n"
       f"alone ({alone.parent.as_posix()}): {w_alone or 'silent'}\n"
       f"beyond the specification, the unedited kind T in place: {w_clean or 'silent'}")

# ---------------------------------------------------------------------------
# 8. The run's four inputs; 9. content-identical to the anchor's embedded copies
# ---------------------------------------------------------------------------
lines, ok = [], True
loaded = {}
for name, kind in RUN_INPUTS:
    path = RUN / name
    handle = cio.read(path, kind)
    attrs = dict(handle.attrs)
    loaded[name.split("_")[-1][:-3]] = handle.load()
    handle.close()
    reduction = D / name.replace("inputs/", "").replace("lindal_closure_", "lindal_")
    h_run, h_red = cio.sha256(path), cio.sha256(reduction)
    good = (attrs["profile_or_run"] == "lindal_closure" and attrs["role"] == "forward"
            and path.name.startswith("lindal_closure_") and h_run != h_red
            and (kind != "composition" or attrs["composition_role"] == "forward"))
    ok = ok and good
    lines.append(f"{path.name:36s} kind {kind:11s} profile_or_run {attrs['profile_or_run']!r}, role "
                 f"{attrs['role']!r}; sha256 {h_run[:16]}... against the reduction's {reduction.name} "
                 f"{h_red[:16]}... differ {h_run != h_red}")
record(8, "casspian-run-inputs wrote the four inputs, each validating under its kind with the run "
          "prefix and role forward, their hashes differing from the reduction's", ok,
       "\n".join(lines) + "\n" + "\n".join(l for l in (FIXTURES / "rebuild.txt").read_text(
           encoding="utf-8").splitlines() if "lindal_closure" in l))

anchor = cio.read(CANDIDATE, "refractivity")
anchor_mem = anchor.load()
anchor.close()
comparison = ctl.check_closure_inputs(loaded, anchor_mem)
lines, ok = [], True
for c in comparison:
    ok = ok and c.identical
    lines.append(f"{c.kind}: identical {c.identical}" + (f"; differences {list(c.differences)}"
                                                         if c.differences else ""))
    for group, names in c.dropped.items():
        lines.append(f"  dropped from {group}: {list(names)}")
record(9, "each run input is content-identical to the candidate kind N's embedded copy under the "
          "comparison of deliverable 2 (root and every group), attributes dropped listed per file",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 10. Beyond the specification: the driver's own refusals
# ---------------------------------------------------------------------------
build_text = (RUN / "lindal_closure_build.toml").read_text(encoding="utf-8")
lines, ok = [], True
for label, mutate in (
        ("a [stage_two] section", lambda s: s + '\n[stage_two]\nrole = "forward"\n'),
        ("a section with role = reduction", lambda s: s.replace('prefix = "lindal_closure"\nrole   = "forward"',
                                                             'prefix = "lindal_closure"\nrole   = "reduction"', 1)),
        ("an output outside inputs/", lambda s: s.replace('output = "inputs/lindal_closure_gravity.nc"',
                                                        'output = "lindal_closure_gravity.nc"'))):
    case = fresh(HERE / "driver_case") / "lindal_closure_build.toml"
    mutated = mutate(build_text)
    assert mutated != build_text, label
    case.write_text(mutated, encoding="utf-8")
    try:
        run_inputs.build(case)
        ok = False
        lines.append(f"{label}: not refused")
    except ctl.ControlFileError as exc:
        lines.append(f"{label}: {exc}")
record(10, "beyond the specification: casspian-run-inputs refuses a build file that is not a run's",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 11. The namelist loads and resolves; the run's inputs load against the candidate anchor
# ---------------------------------------------------------------------------
committed = ctl.read_run_namelist(RUN / "lindal_closure.toml")
base_text = (RUN / "lindal_closure.toml").read_text(encoding="utf-8")
ANCHOR_LINE = 'path   = "../../occul_data/lindal/lindal_refractivity.nc"'
assert base_text.count(ANCHOR_LINE) == 1


def run_case(name, text, composition_edit=None, extra_inputs=()):
    """A run directory named lindal_closure under cases/<name>, inputs copied, anchor at the candidate."""
    directory = fresh(HERE / "cases" / name / "lindal_closure")
    (directory / "inputs").mkdir()
    for input_name, _ in RUN_INPUTS:
        shutil.copy2(RUN / input_name, directory / input_name)
    for source, target in extra_inputs:
        shutil.copy2(source, directory / target)
    if composition_edit:
        composition_edit(directory / "inputs/lindal_closure_composition.nc")
    relative = Path(os.path.relpath(CANDIDATE, directory.resolve()))
    namelist = directory / "lindal_closure.toml"
    namelist.write_text(text.replace(ANCHOR_LINE, f'path   = "{relative.as_posix()}"'), encoding="utf-8")
    return namelist


good_case = run_case("resolves", base_text)
parsed = ctl.read_run_namelist(good_case)
inputs = ctl.load_run_inputs(parsed)
record(11, "the Lindal closure namelist loads and resolves, and the run's inputs load against the "
           "candidate anchor",
       committed.mode == "closure" and inputs.gauge_level_index == 29
       and all(c.identical for c in inputs.closure),
       f"committed namelist: run {committed.name!r}, mode {committed.mode!r}, Ls "
       f"{committed.solar_longitude_deg}, p_b_rule {committed.p_b_rule!r}, gauge "
       f"{committed.gauge_isobar_Pa} Pa, anchor {committed.anchors[0].path}\n"
       f"  inputs {[p.name for p in committed.inputs.values()]}; product {committed.product}\n"
       f"case with the anchor at the candidate: gauge level {inputs.gauge_level_index}; closure "
       f"{[(c.kind, c.identical) for c in inputs.closure]}\n"
       f"commits {dict(inputs.commits)}")

# ---------------------------------------------------------------------------
# 12. The thirteen refusal cases
# ---------------------------------------------------------------------------


def edit_composition(path):
    with netCDF4.Dataset(path, "a") as h:
        he, h2 = h.variables["x_He"], h.variables["x_H2"]
        he[:] = he[:] + 0.01
        h2[:] = h2[:] - 0.01


def refused(label, needle, action):
    try:
        action()
        return False, f"{label}: NOT REFUSED"
    except (ctl.ControlFileError, CasspianSchemaError) as exc:
        message = str(exc)
        return needle.lower() in message.lower(), f"{label}: {message}"


def load(namelist):
    return ctl.load_run_inputs(ctl.read_run_namelist(namelist))


def loading_without_relaxation(namelist):
    with relaxed.disabled():
        return load(namelist)


CASES = [
    ("inputs_missing", "[inputs] missing", "[inputs] is missing",
     lambda s: re.sub(r"\[inputs\].*?\n\n", "", s, count=1, flags=re.S), load),
    ("target", "[target] present", "[target] is not accepted in closure mode",
     lambda s: s + "\n[target]\nlatitude_planetocentric_deg = 10.0\n", load),
    ("grid", "[grid] present", "not implemented in this specification",
     lambda s: s + "\n[grid]\ngeopotential_spacing_m2s2 = 5.0e3\n", load),
    ("unknown_key", "an unknown key", "unknown key",
     lambda s: s.replace('mode        = "closure"', 'mode        = "closure"\nflavor      = "vanilla"'), load),
    ("second_anchor", "a second [[anchors]] entry", "exactly one [[anchors]] entry",
     lambda s: s.replace("[inputs]", '[[anchors]]\nslug = "lindal"\npath = "x.nc"\nweight = 0.0\n'
                                     'measurement_uncertainty_scale = 1.0\n\n[inputs]', 1), load),
    ("slug", "a slug that does not match", "profile_or_run",
     lambda s: s.replace('slug   = "lindal"', 'slug   = "cassini"'), load),
    ("p_b_both", "p_b_Pa alongside p_b_rule", "exactly one",
     lambda s: s.replace('p_b_location = "top_of_anchor_profile"',
                         'p_b_Pa       = 20.0\np_b_location = "top_of_anchor_profile"'), load),
    ("gauge_9000", "a gauge isobar of 9,000 Pa", "not a tabulated level",
     lambda s: s.replace("gauge_isobar_Pa = 1.0e4", "gauge_isobar_Pa = 9.0e3"), load),
    ("dirty_anchor", "a -dirty anchor", "-dirty",
     lambda s: s, loading_without_relaxation),
    ("unprefixed_input", "an input without the run prefix", "run prefix",
     lambda s: s.replace('wind        = "inputs/lindal_closure_wind.nc"',
                         'wind        = "inputs/lindal_wind.nc"'), load),
    ("geodesy", "a namelist that names a geodesy file", "geodesy",
     lambda s: s.replace('wind        = "inputs/lindal_closure_wind.nc"',
                         'wind        = "inputs/lindal_closure_wind.nc"\n'
                         'geodesy     = "inputs/lindal_closure_geodesy.nc"'), load),
    ("composition_x_He", "a composition copy with x_He raised and x_H2 lowered by 0.01", "x_He",
     lambda s: s, load),
    ("season", "a namelist whose solar_longitude_deg differs from the anchor's", "not a closure",
     lambda s: s.replace("solar_longitude_deg = 18.2", "solar_longitude_deg = 18.3"), load),
]
lines, ok = [], True
for name, label, needle, mutate, action in CASES:
    mutated = mutate(base_text)
    if name not in ("dirty_anchor", "composition_x_He"):
        assert mutated != base_text, name
    extra = ((RUN / "inputs/lindal_closure_wind.nc", "inputs/lindal_wind.nc"),) if name == "unprefixed_input" else ()
    namelist = run_case(name, mutated, edit_composition if name == "composition_x_He" else None, extra)
    good, text_line = refused(label, needle, lambda: action(namelist))
    ok = ok and good
    lines.append(text_line)
record(12, "thirteen refusal cases, each refused with the right message (the composition case needed "
           "only the x_He and x_H2 edits: the kind C reader holds no closure-declaration check)",
       ok and len(CASES) == 13, "\n".join(lines))

# ---------------------------------------------------------------------------
# 13. A synthetic kind profile validates and reads back; its three refusals
# ---------------------------------------------------------------------------
synthetic_dir = fresh(HERE / "synthetic")
root = anchor_mem.to_dataset(inherit=False)
thermo = anchor_mem["inputs/thermo"].to_dataset(inherit=False)
p_tab = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
T_tab = np.asarray(thermo["temperature_K"].values, dtype="float64")
r = np.asarray(root["radius_m"].values, dtype="float64")
h = np.asarray(root["height_above_anchor_isobar_m"].values, dtype="float64")
N = np.asarray(root["refractivity"].values, dtype="float64")
R_bar = np.asarray(root["mean_refractivity_m3"].values, dtype="float64")
M_bar = np.asarray(root["mean_molar_mass_kg_mol"].values, dtype="float64")
phi_c = math.radians(float(root["latitude_planetocentric_deg"].values))
gravity, rotation = inputs.gravity, inputs.rotation
constants = (float(rotation["angular_rate_rad_s"]), float(gravity["GM_m3s2"]),
             np.asarray(gravity["J"].values), np.asarray(gravity["degree"].values),
             float(gravity["normalization_radius_m"]))
u = float(wind_of_latitude(inputs.wind)(np.array([phi_c]))[0])
gauge = inputs.gauge_level_index
Phi = gp.geopotential_along_profile(u, r, h, phi_c, gauge, *constants).geopotential_m2s2
# A synthetic residual, so that F6 has something to draw; not the Step 4 production.
p_syn = p_tab * (1.0 + 5.0e-4 * np.sin(np.arange(p_tab.size)))
p_syn[0] = p_tab[0]
n_syn = N / R_bar
T_syn = p_syn * R_bar / (BOLTZMANN_CONSTANT * N)
UNSTATED = ("composition: not propagated through the integral; anchor_radius: not propagated; "
            "label_latitude: not propagated. A first order companion carrying only the terms that "
            "pass through the integral would be read as a total; the Monte Carlo wrapper supplies "
            "the product's uncertainty (SPEC_03 decision 5)")


def a(units, long_name, provenance, **extra):
    return {"units": units, "long_name": long_name, "provenance": provenance, **extra}


def companion(units, name):
    return (("level",), np.full(p_tab.size, np.nan),
            a(units, f"uncertainty on {name}", "modeled", uncertainty_kind="1sigma",
              uncertainty_method="not propagated", uncertainty_terms_unstated=UNSTATED))


def profile_dataset():
    ds = xr.Dataset(
        {
            "radius_m": (("level",), r, a("m", "radius, copied from the anchor", "derived")),
            "height_above_anchor_isobar_m": (("level",), h, a("m", "height above the anchor isobar, copied", "measured")),
            "refractivity": (("level",), N, a("1", "refractivity, copied", "derived")),
            "mean_refractivity_m3": (("level",), R_bar, a("m3", "mean refractivity per molecule", "modeled")),
            "mean_refractivity_uncertainty_m3": companion("m3", "mean_refractivity_m3"),
            "mean_molar_mass_kg_mol": (("level",), M_bar, a("kg mol-1", "mean molar mass", "modeled")),
            "mean_molar_mass_uncertainty_kg_mol": companion("kg mol-1", "mean_molar_mass_kg_mol"),
            "number_density_m3": (("level",), n_syn, a("m-3", "number density, synthetic", "modeled")),
            "number_density_uncertainty_m3": companion("m-3", "number_density_m3"),
            "pressure_Pa": (("level",), p_syn, a("Pa", "pressure, synthetic acceptance values", "modeled",
                                                 boundary_pressure_Pa=float(p_tab[0]))),
            "pressure_uncertainty_Pa": companion("Pa", "pressure_Pa"),
            "temperature_K": (("level",), T_syn, a("K", "temperature, synthetic", "modeled")),
            "temperature_uncertainty_K": companion("K", "temperature_K"),
            "pressure_tabulated_Pa": (("level",), p_tab, a("Pa", "tabulated pressure of the anchor", "derived")),
            "temperature_tabulated_K": (("level",), T_tab, a("K", "tabulated temperature of the anchor", "derived")),
            "geopotential_uncertainty_m2s2": companion("m2 s-2", "geopotential_m2s2"),
            "latitude_planetocentric_deg": ((), float(root["latitude_planetocentric_deg"].values),
                                            a("degrees_north", "latitude of the anchor", "derived")),
            "psi_deg": ((), float(root["psi_deg"].values), a("degrees", "tilt, copied", "derived")),
            "gauge_isobar_Pa": ((), 1.0e4, a("Pa", "gauge isobar", "index")),
            "gauge_level_index": ((), np.int32(gauge), a("1", "gauge level", "index")),
            "boundary_pressure_Pa": ((), float(p_tab[0]), a("Pa", "boundary pressure", "index")),
            "boundary_level_index": ((), np.int32(0), a("1", "boundary level", "index")),
        },
        coords={"geopotential_m2s2": (("level",), Phi, a("m2 s-2", "geopotential, gauge at 1e4 Pa", "modeled",
                                                       positive="up"))},
    )
    ds.attrs = {
        "title": "synthetic kind profile, SPEC_03 Step 3 acceptance only",
        "profile_or_run": "lindal_closure", "role": "forward",
        "source": "SPEC_03 Step 3 acceptance; synthetic pressure, not a production",
        "epoch": "none: season declared as solar longitude only",
        "solar_longitude_deg": 18.2, "solar_longitude_source": "declared in the namelist",
        "anchor_solar_longitudes_deg": np.array([18.2]),
        "boundary_pressure_Pa": float(p_tab[0]),
        "latitude_planetocentric_absent_meaning": "point",
    }
    return ds


synthetic = synthetic_dir / "lindal_closure_profile.nc"
groups = {}
for node in anchor_mem.subtree:
    groups["anchors/lindal" + (node.path if node.path != "/" else "")] = node.to_dataset(inherit=False)
for key in ctl.RUN_INPUT_KINDS:
    obj = getattr(inputs, key)
    if isinstance(obj, xr.DataTree):
        for node in obj.subtree:
            groups[f"inputs/{key}" + (node.path if node.path != "/" else "")] = node.to_dataset(inherit=False)
    else:
        groups[f"inputs/{key}"] = obj
groups["namelist"] = xr.Dataset(attrs={"text": parsed.text, "sha256": parsed.sha256})
groups["production_record"] = xr.Dataset(attrs={"mode": "closure", "note": "synthetic, acceptance only"})
ds = profile_dataset()
ds.attrs["input_hashes"] = cio.input_hashes(
    [CANDIDATE] + [good_case.parent / n for n, _ in RUN_INPUTS] + [good_case], synthetic)
cio.write(synthetic, ds, "profile", groups=groups, created_by="accept_step03_3")
back = cio.read(synthetic, "profile")
is_tree = isinstance(back, xr.DataTree)
back_groups = sorted(cio.group_paths(synthetic))
back.close()

lines, ok = [], is_tree
variants = []
v = profile_dataset(); v = v.drop_vars("temperature_tabulated_K"); variants.append(("pressure_tabulated_Pa without temperature_tabulated_K", v, "temperature_tabulated_K"))
v = profile_dataset(); values = v["geopotential_m2s2"].values.copy(); values[5], values[6] = values[6], values[5]
v = v.assign_coords(geopotential_m2s2=(("level",), values, v["geopotential_m2s2"].attrs)); variants.append(("a coordinate that is not monotonic", v, "monotonic"))
v = profile_dataset(); v = v.drop_vars("pressure_uncertainty_Pa"); variants.append(("a modeled variable without its NaN companion", v, "companion"))
for label, variant, needle in variants:
    variant.attrs["input_hashes"] = ds.attrs["input_hashes"]
    try:
        cio.write(synthetic_dir / "must_not_exist.nc", variant, "profile", groups=groups,
                  created_by="accept_step03_3")
        ok = False
        lines.append(f"{label}: not refused")
    except CasspianSchemaError as exc:
        ok = ok and needle in str(exc)
        lines.append(f"{label}: {exc}")
record(13, "a synthetic kind profile validates and reads back as a DataTree, and is refused without "
           "temperature_tabulated_K, with a non-monotonic coordinate, or without a NaN companion",
       ok, f"written {synthetic.as_posix()}; read back as {'DataTree' if is_tree else 'NOT a DataTree'}; "
           f"{len(back_groups)} groups, first {back_groups[:3]}\n" + "\n".join(lines))

# ---------------------------------------------------------------------------
# 14. casspian-plots: F5 and F6 on the profile, F1 to F4 on kind N with nothing skipped
# ---------------------------------------------------------------------------
cli_profile = subprocess.run(["casspian-plots", str(synthetic), "--out", str(synthetic_dir / "figures")],
                             capture_output=True, text=True)
cli_n = subprocess.run(["casspian-plots", str(CANDIDATE), "--out", str(HERE / "n_figures")],
                       capture_output=True, text=True)
rp = render(synthetic, synthetic_dir / "figures_render")
rn = render(CANDIDATE, HERE / "n_figures_render")
F6 = rp.data.get("F6", {})
band = np.asarray(F6.get("envelope", []))
p_term = np.asarray(F6.get("envelope_pressure_term", []))
h_term = np.asarray(F6.get("envelope_height_term", []))
names_p = [p.name for p in rp.written]
names_n = [p.name for p in rn.written]
ok = (cli_profile.returncode == 0 and cli_n.returncode == 0
      and names_p == ["lindal_closure_diag_F5_geopotential.png", "lindal_closure_diag_F6_hydrostatic.png"]
      and band.size and np.all(np.isfinite(band)) and np.all(band > 0)
      and np.all(p_term[:-1] == 0.0) and p_term[-1] > 0.0
      and names_n == [f"lindal_diag_F{k}_{m}.png" for k, m in
                      ((1, "inputs"), (2, "gravity_profile"), (3, "gravity_latitude"), (4, "product"))]
      and not rn.skipped and "skipped" not in cli_n.stdout)
record(14, "casspian-plots renders F5 and F6 with the envelope band on the synthetic profile, and F1 to "
           "F4 with nothing skipped on kind N",
       ok,
       f"profile, CLI exit {cli_profile.returncode}: {cli_profile.stdout.strip()}\n"
       f"profile, render: {names_p}; skipped {rp.skipped}\n"
       f"F6 envelope: pressure term zero at {int(np.sum(p_term == 0))} of {p_term.size} levels, "
       f"{p_term[-1]:.3e} at the bottom row; height term {h_term.min():.3e} to {h_term.max():.3e}; "
       f"band {band.min():.3e} to {band.max():.3e}; scale height "
       f"{np.asarray(F6.get('scale_height_m')).min() / 1e3:.1f} to {np.asarray(F6.get('scale_height_m')).max() / 1e3:.1f} km\n"
       f"kind N, CLI exit {cli_n.returncode}: {cli_n.stdout.strip()}\n"
       f"kind N, render: {names_n}; skipped {rn.skipped or 'nothing'}")

# ---------------------------------------------------------------------------
# 15. Beyond the specification: the rebuild changed no number
# ---------------------------------------------------------------------------
old = cio.read(BEFORE / "lindal_refractivity.nc", "refractivity") if False else None
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    old_tree = xr.open_datatree(BEFORE / "lindal_refractivity.nc", engine="netcdf4").load()
new_tree = anchor_mem
old_root, new_root = old_tree.to_dataset(inherit=False), new_tree.to_dataset(inherit=False)
lines, ok = [], True
differ = [name for name in old_root.variables
          if not np.array_equal(np.asarray(old_root[name].values), np.asarray(new_root[name].values), equal_nan=True)]
old_rec, new_rec = dict(old_tree["reduction_record"].attrs), dict(new_tree["reduction_record"].attrs)
skip = {"casspian_git_commit"} | {k for k in old_rec if k.startswith("input_sha256_")}
rec_differ = [k for k in old_rec if k not in skip and not np.array_equal(np.asarray(old_rec[k]), np.asarray(new_rec.get(k)))]
ok = not differ and not rec_differ
lines.append(f"kind N variables compared {len(old_root.variables)}; differing {differ or 'none'}")
lines.append(f"reduction_record attributes compared {len(old_rec) - len(skip)}; differing {rec_differ or 'none'}")
for name, kind in REDUCTION:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        o = xr.open_datatree(BEFORE / name, engine="netcdf4").load()
        nw = xr.open_datatree(D / name, engine="netcdf4").load()
    bad = []

    def same(x, y):
        x, y = np.asarray(x), np.asarray(y)
        if x.dtype.kind in "fc" and y.dtype.kind in "fc":
            return np.array_equal(x, y, equal_nan=True)
        return x.shape == y.shape and bool(np.all(x == y))

    for node in o.subtree:
        a_ds = node.to_dataset(inherit=False)
        b_ds = nw[node.path].to_dataset(inherit=False) if node.path != "/" else nw.to_dataset(inherit=False)
        for var in a_ds.variables:
            if var not in b_ds.variables or not same(a_ds[var].values, b_ds[var].values):
                bad.append(f"{node.path}:{var}")
    ok = ok and not bad
    added = sorted(set(nw.attrs) - set(o.attrs))
    removed = sorted(set(o.attrs) - set(nw.attrs))
    lines.append(f"{name:26s} values differing {bad or 'none'}; attributes added {added}, removed {removed or 'none'}")
record(15, "beyond the specification: the rebuild changed no value in any product; kind N and its "
           "reduction_record are unchanged bit for bit", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 16. Beyond the specification: commit and SHA-256 of every file, provisional until the sweep
# ---------------------------------------------------------------------------
rows = []
for path, kind in FILES:
    with netCDF4.Dataset(path) as h_:
        commit = h_.getncattr("casspian_git_commit")
    rows.append(f"{path.as_posix():58s} {kind:12s} {cio.sha256(path)}  {commit}")
rows.append(f"{(D / 'lindal_reduction.toml').as_posix():58s} {'manifest':12s} {cio.sha256(D / 'lindal_reduction.toml')} "
            f"(before {cio.sha256(BEFORE / 'lindal_reduction.toml')})")
record(16, "beyond the specification: commit and SHA-256 of every file, provisional until the sweep",
       True, "\n".join(rows))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
