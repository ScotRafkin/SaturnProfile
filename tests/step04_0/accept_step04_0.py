"""Acceptance checks for SPEC_04 v0.9 Step 0.

Kind W in three parts, kind C on a latitude grid, `forward/lindal_transfer/`, the transfer
namelist for M anchors, `load_run_inputs` in transfer mode, the anchor as it arrives, and the
propagation hook. Every check prints its measured value; a check beyond the specification is
labeled so.

**The relaxation.** This step changes input files, so the products under `occul_data/` and
`forward/` were rebuilt on this working tree by `rebuild_chain.py` and carry
`casspian_git_commit` ending in `-dirty`. `forward` refuses a `-dirty` input, so this script
relaxes the one refusal point, `casspian.lib.control._refuse_dirty_commit`, for its own
duration and nowhere else, and check 0 shows the refusal firing before the relaxation is put
in place. Nothing else is patched, no commit string is altered on disk, and every other
refusal below is exercised in full. The registered products are rebuilt from the clean tree at
the sweep and their hashes recorded in the report then.

**The swept copies.** The pre-step products, the form of every file before this step, were
rebuilt by `build_swept.py` in a git worktree at `5f23a8b` whose tree is clean, so they carry
that commit and not `-dirty`. Nothing in the reduction chain changed between the registered
sweep at `6ae113e` and `5f23a8b`, so their values are the swept values. `SWEPT` below points at
their committed copies, `tests/step04_0/fixtures/swept/`.

Run from the repository root.
"""

import os
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

from casspian.forward import propagate as fp
from casspian.lib import composition as ctl_comp
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError

#: The Lindal anchor's own planetocentric latitude, the column the reduction reads (decision L).
PHI_C_DEG = 30.80556842739218

ROOT = Path(".").resolve()
HERE = Path("reports/step04_0")
LINDAL = Path("occul_data/lindal")
CLOSURE = Path("forward/lindal_closure")
TRANSFER = Path("forward/lindal_transfer")
SWEPT = Path(os.environ.get(
    "CASSPIAN_SWEPT",
    Path(__file__).resolve().parent / "fixtures" / "swept"))

results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def nodes(tree):
    """`{group path: Dataset}` for a DataTree, the root as '/'."""
    out = {}
    base = tree.path.rstrip("/")
    for node in tree.subtree:
        relative = node.path[len(base):] if base else node.path
        out[relative or "/"] = node.to_dataset(inherit=False)
    return out


def open_any(path):
    """Read a netCDF file into memory as a DataTree, groups and all, with no schema check.

    The values are pulled in and the file released, as `lib.control._load_into_memory` does:
    this script opens every product twice over and netCDF4 does not survive that many live
    handles on Windows.
    """
    handle = xr.open_datatree(path, engine="netcdf4")
    try:
        return handle.load()
    finally:
        handle.close()


def refusal(callable_, *args, **kwargs):
    """Run and return the refusal message, or None when nothing was refused."""
    try:
        callable_(*args, **kwargs)
    except Exception as exc:                      # the refusals are of several types
        return f"{type(exc).__name__}: {exc}"
    return None


HERE.mkdir(parents=True, exist_ok=True)
WORK = HERE / "cases"
if WORK.exists():
    shutil.rmtree(WORK)
WORK.mkdir()

# ---------------------------------------------------------------------------
# 0. The -dirty refusal fires on this working tree, and is then relaxed, named
# ---------------------------------------------------------------------------
NAMELIST = TRANSFER / "lindal_transfer.toml"

# The check's subject is a `-dirty` copy this script writes, not whatever state the working tree
# happens to be in (SPEC_04 section 12 ruling 4, on finding 4 of REPORT_04_step1). At Step 0 the
# tree's own products were all `-dirty`, so asking the tree was the same question; the Step 0 sweep
# then rebuilt them clean, and the check had nothing left to fire on. With its own subject it means
# "the refusal fires on a `-dirty` input" whatever the tree holds.
DIRTY = WORK / "dirty"
DIRTY.mkdir()
dirty_input = DIRTY / "lindal_transfer_wind.nc"
shutil.copy(TRANSFER / "inputs/lindal_transfer_wind.nc", dirty_input)
with netCDF4.Dataset(dirty_input, "a") as handle:
    clean_commit = str(handle.getncattr("casspian_git_commit"))
    handle.setncattr("casspian_git_commit", f"{clean_commit}-dirty")
with netCDF4.Dataset(dirty_input) as handle:
    dirty_attrs = {name: handle.getncattr(name) for name in handle.ncattrs()}
on_dirty_copy = refusal(ctl._refuse_dirty_commit, dirty_input, dirty_attrs, "casspian-forward")
on_clean_file = refusal(ctl._refuse_dirty_commit, TRANSFER / "inputs/lindal_transfer_wind.nc",
                        {"casspian_git_commit": clean_commit}, "casspian-forward")
on_this_tree = refusal(lambda: ctl.load_run_inputs(ctl.read_run_namelist(NAMELIST)))

_original_refuse = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(
    attrs.get("casspian_git_commit", ""))
record(0, "the -dirty refusal fires on a -dirty input and not on a clean one, and is then relaxed "
          "for this script",
       on_dirty_copy is not None and "-dirty" in str(on_dirty_copy) and on_clean_file is None,
       f"subject: {dirty_input}, a copy of the registered transfer wind with its commit rewritten "
       f"{clean_commit[:12]} -> {clean_commit[:12]}-dirty\n"
       f"refused on the -dirty copy: {on_dirty_copy}\n"
       f"refused on the same file clean: {on_clean_file}\n"
       f"for information, not asserted: loading the run from this working tree gave "
       f"{on_this_tree}; the tree's products are clean after the Step 0 sweep, so the refusal has "
       f"nothing to fire on there (SPEC_04 section 12 ruling 4)\n"
       "RELAXATION: casspian.lib.control._refuse_dirty_commit is replaced for this script only "
       "(SPEC_03 section 0, in-memory candidate rule). Nothing else is patched.")

WIND_FILES = {
    "reduction": (LINDAL / "lindal_wind.nc", SWEPT / "occul_data/lindal/lindal_wind.nc"),
    "closure": (CLOSURE / "inputs/lindal_closure_wind.nc",
                SWEPT / "forward/lindal_closure/inputs/lindal_closure_wind.nc"),
    "transfer": (TRANSFER / "inputs/lindal_transfer_wind.nc", None),
}

# ---------------------------------------------------------------------------
# 1. Kind W in three parts, the retired items gone
# ---------------------------------------------------------------------------
lines, ok = [], True
for name, (path, _) in WIND_FILES.items():
    wind = cio.read(path, "wind")
    try:
        present = set(wind.variables)
        required = {"u_reference_ms", "u_total_ms", "u_shear_ms", "u_total_uncertainty_ms"}
        retired_vars = {"u_cylindrical_ms"} & present
        retired_attrs = {"decomposition", "decomposition_geometry"} & set(wind.attrs)
        dims = tuple(wind["u_reference_ms"].dims)
        good = (required <= present and not retired_vars and not retired_attrs
                and dims == ("latitude_planetocentric",))
        ok &= good
        lines.append(f"{name}: {path.name} u_reference_ms{dims} shape "
                     f"{wind['u_reference_ms'].shape}, u_total_ms"
                     f"{tuple(wind['u_total_ms'].dims)}, u_shear_ms"
                     f"{tuple(wind['u_shear_ms'].dims)}; retired variables "
                     f"{sorted(retired_vars) or 'none'}, retired attributes "
                     f"{sorted(retired_attrs) or 'none'}")
    finally:
        wind.close()
record(1, "kind W is data in three parts and the retired items are gone, all three files",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. u_total unchanged to the bit against the swept copies
# ---------------------------------------------------------------------------
lines, ok = [], True
for name, (path, swept) in WIND_FILES.items():
    if swept is None:
        continue
    new, old = open_any(path).to_dataset(inherit=False), open_any(swept).to_dataset(inherit=False)
    same = {}
    for var in ("u_total_ms", "u_total_uncertainty_ms", "value_provenance",
                "latitude_planetocentric_deg", "pressure_Pa", "reference_level_pressure_Pa"):
        same[var] = bool(np.array_equal(np.asarray(new[var].values),
                                        np.asarray(old[var].values), equal_nan=True))
    reference_column = int(np.flatnonzero(
        np.asarray(new["pressure_Pa"].values) == float(new["reference_level_pressure_Pa"]))[0])
    from_old = bool(np.array_equal(
        np.asarray(new["u_reference_ms"].values),
        np.asarray(old["u_cylindrical_ms"].values)[:, reference_column]))
    ok &= all(same.values()) and from_old
    lines.append(f"{name}: " + ", ".join(f"{k} identical {v}" for k, v in same.items())
                 + f"; u_reference_ms equals the swept u_cylindrical_ms at the reference "
                   f"column {reference_column}: {from_old}")
# the transfer wind carries the same field as the closure wind, which is the point of the state
new = open_any(WIND_FILES["transfer"][0]).to_dataset(inherit=False)
old = open_any(WIND_FILES["closure"][0]).to_dataset(inherit=False)
same_field = bool(np.array_equal(np.asarray(new["u_total_ms"].values),
                                 np.asarray(old["u_total_ms"].values)))
ok &= same_field
lines.append(f"transfer: u_total_ms equals the closure run's, bit for bit: {same_field}")
record(2, "u_total and its uncertainty are unchanged to the bit against the swept copies",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 3. The sum identity and the shear at the reference level
# ---------------------------------------------------------------------------
lines, ok = [], True
for name, (path, _) in WIND_FILES.items():
    d = open_any(path).to_dataset(inherit=False)
    total = np.asarray(d["u_total_ms"].values)
    reference = np.asarray(d["u_reference_ms"].values)
    shear = np.asarray(d["u_shear_ms"].values)
    worst = float(np.nanmax(np.abs(total - (reference[:, None] + shear))))
    column = int(np.flatnonzero(
        np.asarray(d["pressure_Pa"].values) == float(d["reference_level_pressure_Pa"]))[0])
    at_reference = float(np.max(np.abs(shear[:, column])))
    ok &= (worst == 0.0 and at_reference == 0.0)
    lines.append(f"{name}: worst |u_total - (u_reference + u_shear)| = {worst:.3e} m/s; "
                 f"max |u_shear| at the reference level = {at_reference:.3e} m/s")
record(3, "the sum identity holds exactly and the shear vanishes at the reference level",
       ok, "\n".join(lines))

# 4. Retired by SPEC_08 decision 3: kind W no longer checks the sum identity.

# ---------------------------------------------------------------------------
# 5. Kind C on a latitude grid, every column the point file's
# ---------------------------------------------------------------------------
grid = open_any(TRANSFER / "inputs/lindal_transfer_composition.nc").to_dataset(inherit=False)
point = open_any(SWEPT / "occul_data/lindal/lindal_composition.nc").to_dataset(inherit=False)
lines, ok = [], True
lat = np.asarray(grid["latitude_planetocentric_deg"].values)
lines.append(f"latitude nodes {lat.size} from {lat[0]:g} to {lat[-1]:g} deg, step "
             f"{np.diff(lat)[0]:g}; unique steps {np.unique(np.round(np.diff(lat), 12)).tolist()}")
ok &= (lat.size == 181 and lat[0] == -90.0 and lat[-1] == 90.0)
for var in ("x_H2", "x_He", "x_NH3", "x_H2_uncertainty", "nh3_provenance"):
    values = np.asarray(grid[var].values)
    equal_columns = bool(np.all(values == values[:, :1]))
    equals_point = bool(np.array_equal(values[:, 0], np.asarray(point[var].values),
                                       equal_nan=True))
    ok &= equal_columns and equals_point
    lines.append(f"{var}{tuple(grid[var].dims)}: every column equal {equal_columns}; "
                 f"column equals the point file's bit for bit {equals_point}")
pressure_1d = tuple(grid["pressure_Pa"].dims) == ("level",)
same_p = bool(np.array_equal(np.asarray(grid["pressure_Pa"].values),
                             np.asarray(point["pressure_Pa"].values)))
absent = "latitude_planetocentric_absent_meaning" in grid.attrs
ok &= pressure_1d and same_p and not absent
lines.append(f"pressure_Pa stays one dimensional {pressure_1d} and unchanged {same_p}; "
             f"latitude_planetocentric_absent_meaning present {absent}")
lines.append("latitude_grid_rule: " + str(grid.attrs.get("latitude_grid_rule", "absent")))
# Decision O: kind C is one structure for every use, so the reduction's and the closure's
# compositions are fields too, and the reduction reads the column at phi_c by decision L.
for name, path in (("reduction", LINDAL / "lindal_composition.nc"),
                   ("closure", CLOSURE / "inputs/lindal_closure_composition.nc")):
    field = open_any(path)
    root = field.to_dataset(inherit=False)
    gridded = ("latitude_planetocentric" in root.sizes
               and "latitude_planetocentric_absent_meaning" not in root.attrs)
    column = ctl_comp.column_at(field, PHI_C_DEG).to_dataset(inherit=False)
    same = all(bool(np.array_equal(np.asarray(column[v].values),
                                   np.asarray(point[v].values), equal_nan=True))
               for v in ("x_H2", "x_He", "x_NH3", "x_H2_uncertainty", "nh3_provenance"))
    ok &= gridded and same
    lines.append(f"the {name} composition is a field on {dict(root.sizes)} with no absence "
                 f"marker: {gridded}; its column at phi_c = {PHI_C_DEG:.6f} deg equals the "
                 f"pre-step point file bit for bit in every species variable: {same}")
outside = refusal(ctl_comp.column_at, open_any(LINDAL / "lindal_composition.nc"), 95.0)
ok &= outside is not None
lines.append(f"a latitude outside the grid is refused by the interpolant: {outside}")
record(5, "kind C is one structure, a field on (level, latitude), carrying the source's column "
          "at every node", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 6. The cascade: every product rebuilt, no value changed
# ---------------------------------------------------------------------------
PAIRS = [
    (LINDAL / "raw/lindal_raw.nc", SWEPT / "occul_data/lindal/raw/lindal_raw.nc"),
    (LINDAL / "lindal_gravity.nc", SWEPT / "occul_data/lindal/lindal_gravity.nc"),
    (LINDAL / "lindal_rotation.nc", SWEPT / "occul_data/lindal/lindal_rotation.nc"),
    (LINDAL / "lindal_wind.nc", SWEPT / "occul_data/lindal/lindal_wind.nc"),
    (LINDAL / "lindal_composition.nc", SWEPT / "occul_data/lindal/lindal_composition.nc"),
    (LINDAL / "lindal_thermo.nc", SWEPT / "occul_data/lindal/lindal_thermo.nc"),
    (LINDAL / "lindal_geodesy.nc", SWEPT / "occul_data/lindal/lindal_geodesy.nc"),
    (LINDAL / "lindal_refractivity.nc",
     SWEPT / "occul_data/lindal/lindal_refractivity.nc"),
    (CLOSURE / "inputs/lindal_closure_gravity.nc",
     SWEPT / "forward/lindal_closure/inputs/lindal_closure_gravity.nc"),
    (CLOSURE / "inputs/lindal_closure_rotation.nc",
     SWEPT / "forward/lindal_closure/inputs/lindal_closure_rotation.nc"),
    (CLOSURE / "inputs/lindal_closure_wind.nc",
     SWEPT / "forward/lindal_closure/inputs/lindal_closure_wind.nc"),
    (CLOSURE / "inputs/lindal_closure_composition.nc",
     SWEPT / "forward/lindal_closure/inputs/lindal_closure_composition.nc"),
    (CLOSURE / "output/lindal_closure_profile.nc",
     SWEPT / "forward/lindal_closure/output/lindal_closure_profile.nc"),
]
#: The differences this step is allowed to make: the three parts of kind W replacing the
#: cylindrical split, and the provenance every rebuild rewrites. Anything else is a changed
#: value and fails the check.
#:
#: `u_column_ms` and its NaN companion are SPEC_04 Step 1 deliverable 3's addition to kind
#: `profile`, made after this step was accepted. They are named here for the same reason the
#: first two are: the check is for values that changed, and a variable the specification
#: required a later step to add is not that. Nothing else is admitted, so a value that does
#: change still fails. Added under REPORT_04_step1 finding 5.
ALLOWED_ADDED = {"u_reference_ms", "latitude_planetocentric_deg",
                 "u_column_ms", "u_column_uncertainty_ms"}
ALLOWED_REMOVED = {"u_cylindrical_ms"}
ALLOWED_ATTRS_REMOVED = {"decomposition", "decomposition_geometry",
                         "latitude_planetocentric_absent_meaning"}
#: Decision O gives every kind C variable a latitude axis it did not have. That is a change of
#: shape, not of value: the check below requires every column of the new variable to equal the
#: old one bit for bit, which is the statement that the field is uniform and unchanged.
LATITUDE_DIM = "latitude_planetocentric"


def same_values(new, old) -> bool:
    """Whether `new` carries the same values as `old`, allowing a gained latitude axis."""
    a, b = np.asarray(new.values), np.asarray(old.values)
    equal_nan = a.dtype.kind == "f"
    if LATITUDE_DIM in new.dims and LATITUDE_DIM not in old.dims:
        axis = new.dims.index(LATITUDE_DIM)
        moved = np.moveaxis(a, axis, -1)
        return bool(np.all([np.array_equal(moved[..., j], b, equal_nan=equal_nan)
                            for j in range(moved.shape[-1])]))
    return bool(np.array_equal(a, b, equal_nan=equal_nan))


#: `production_record` records which attributes the closure comparison dropped for each input.
#: Retiring `decomposition_geometry` from kind W leaves one fewer to drop, so this attribute
#: changes by the retirement and not by any value. It is printed below rather than waved past.
CLOSURE_DROPPED_RECORD = "closure_dropped_wind"


def provenance_attr(name, value):
    return (name in ctl.CLOSURE_DROPPED_ATTRIBUTES or "sha256:" in str(value)
            or name.startswith("input_sha256_") or name == CLOSURE_DROPPED_RECORD
            # Decision O's prose about the grid, and the history line that names it.
            or name in ("latitude_grid_rule", "closure_dropped_composition"))


lines, ok = [], True
for new_path, old_path in PAIRS:
    new_tree, old_tree = nodes(open_any(new_path)), nodes(open_any(old_path))
    unexpected = []
    if set(new_tree) != set(old_tree):
        unexpected.append(f"groups {sorted(set(new_tree) ^ set(old_tree))}")
    added_all, removed_all, attrs_removed_all = set(), set(), set()
    for group in sorted(set(new_tree) & set(old_tree)):
        a, b = new_tree[group], old_tree[group]
        added, removed = set(a.variables) - set(b.variables), set(b.variables) - set(a.variables)
        added_all |= added
        removed_all |= removed
        if added - ALLOWED_ADDED:
            unexpected.append(f"{group}: unexpected new variables {sorted(added - ALLOWED_ADDED)}")
        if removed - ALLOWED_REMOVED:
            unexpected.append(f"{group}: variables gone {sorted(removed - ALLOWED_REMOVED)}")
        for var in sorted(set(a.variables) & set(b.variables)):
            if not same_values(a[var], b[var]):
                unexpected.append(f"{group}: {var} changed value")
        gone = set(b.attrs) - set(a.attrs)
        attrs_removed_all |= gone
        if gone - ALLOWED_ATTRS_REMOVED:
            unexpected.append(f"{group}: attributes gone {sorted(gone - ALLOWED_ATTRS_REMOVED)}")
        for key in sorted(set(a.attrs) & set(b.attrs)):
            if provenance_attr(key, b.attrs[key]):
                continue
            if not ctl._same_value(a.attrs[key], b.attrs[key]):
                unexpected.append(f"{group}: attribute {key} changed")
    ok &= not unexpected
    lines.append(f"{new_path.as_posix()}: variables added {sorted(added_all) or 'none'}, "
                 f"removed {sorted(removed_all) or 'none'}, attributes removed "
                 f"{sorted(attrs_removed_all) or 'none'}; unexpected differences "
                 f"{unexpected or 'none'}")
    for group in sorted(set(new_tree) & set(old_tree)):
        a, b = new_tree[group], old_tree[group]
        if CLOSURE_DROPPED_RECORD in a.attrs and not ctl._same_value(
                a.attrs[CLOSURE_DROPPED_RECORD], b.attrs.get(CLOSURE_DROPPED_RECORD)):
            lines.append(f"    {group}: {CLOSURE_DROPPED_RECORD} was "
                         f"{b.attrs.get(CLOSURE_DROPPED_RECORD)!r}, now "
                         f"{a.attrs[CLOSURE_DROPPED_RECORD]!r}; the retired attribute is no "
                         "longer there to drop")
record(6, "the cascade rebuilt every product and changed no value, group by group",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 7. casspian-run-inputs writes the transfer run's four inputs
# ---------------------------------------------------------------------------
lines, ok = [], True
for key, kind in ctl.RUN_INPUT_KINDS.items():
    path = TRANSFER / f"inputs/lindal_transfer_{key}.nc"
    exists = path.exists()
    ok &= exists
    if exists:
        handle = cio.read(path, kind)
        try:
            lines.append(f"{path.name}: kind {kind}, role {handle.attrs.get('role')!r}, "
                         f"profile_or_run {handle.attrs.get('profile_or_run')!r}, "
                         f"sha256 {cio.sha256(path)[:16]}...")
        finally:
            handle.close()
record(7, "casspian-run-inputs wrote the transfer run's four inputs, each valid under its kind",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 8. The transfer namelist loads and resolves
# ---------------------------------------------------------------------------
namelist = ctl.read_run_namelist(NAMELIST)
detail = "\n".join([
    f"mode {namelist.mode!r}, season {namelist.solar_longitude_deg}, date {namelist.date!r}",
    f"anchors {[(a.slug, a.weight, a.measurement_uncertainty_scale) for a in namelist.anchors]}",
    f"target {namelist.target_latitude_deg} deg; gauge {namelist.gauge_isobar_Pa} Pa; "
    f"datum {namelist.datum_isobar_Pa} Pa",
    f"grid {dict(namelist.grid)}",
    f"numerics {{" + ", ".join(
        f"{k}: {v if not isinstance(v, type(namelist.grid)) else dict(v)}"
        for k, v in namelist.numerics.items()) + "}",
    f"estimation {dict(namelist.estimation)}",
    f"product {namelist.product.relative_to(ROOT).as_posix()}",
])
record(8, "the transfer namelist loads and resolves", namelist.mode == "transfer", detail)

# ---------------------------------------------------------------------------
# 9. The anchor as it arrives: the two uncertainty columns
# ---------------------------------------------------------------------------
inputs = ctl.load_run_inputs(namelist)
lines, ok = [], True
ok &= len(inputs.anchors) == len(namelist.anchors)
for anchor in inputs.anchors:
    tree = open_any(anchor.path).to_dataset(inherit=False)
    N = np.asarray(tree["refractivity"].values)
    sigma_N = np.asarray(tree["refractivity_uncertainty"].values)
    expected = anchor.measurement_uncertainty_scale * sigma_N / N
    measurement_ok = bool(np.array_equal(anchor.sigma_ln_N_measurement, expected))
    season_zero = bool(np.all(anchor.sigma_ln_N_season == 0.0))
    lnN_ok = bool(np.allclose(np.exp(anchor.ln_N), N, rtol=0, atol=0)
                  or np.array_equal(anchor.ln_N, np.log(N)))
    ok &= measurement_ok and season_zero and anchor.season_term == "absent" and lnN_ok
    lines.append(
        f"{anchor.slug}: phi_c {anchor.latitude_planetocentric_deg:.6f} deg, "
        f"{anchor.ln_N.size} levels, gauge level {anchor.gauge_level_index} at "
        f"{anchor.label_pressure_Pa[anchor.gauge_level_index]:g} Pa, weight {anchor.weight} "
        f"(construction {anchor.is_construction})")
    lines.append(
        f"    sigma_ln_N_measurement {anchor.sigma_ln_N_measurement.min():.6e} to "
        f"{anchor.sigma_ln_N_measurement.max():.6e}, equal to scale * sigma_N / N: "
        f"{measurement_ok}")
    lines.append(
        f"    sigma_ln_N_season all zero {season_zero}, season_term {anchor.season_term!r}; "
        f"anchor season {anchor.anchor_season_deg}, run season {anchor.run_season_deg}, "
        f"matches {anchor.season_matches_run}")
record(9, "the loader returns every anchor as one object with its two uncertainty columns",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 10. The propagation hook
# ---------------------------------------------------------------------------
lines, ok = [], True
for anchor in inputs.anchors:
    moved = fp.propagate(anchor, namelist.solar_longitude_deg)
    unchanged = bool(np.array_equal(moved.ln_N, anchor.ln_N))
    season_zero = bool(np.all(moved.sigma_ln_N_season == 0.0))
    has_record = set(moved.propagation) >= {
        "anchor_solar_longitude_deg", "run_solar_longitude_deg", "propagation", "season_term"}
    ok &= unchanged and season_zero and has_record
    ok &= moved.propagation["propagation"] == fp.PROPAGATION_ABSENT
    lines.append(f"{anchor.slug}: ln N unchanged bit for bit {unchanged}; season term still "
                 f"zero {season_zero}")
    lines.append(f"    record {dict(moved.propagation)}")
# a run at another season passes through the same path, and is recorded, not refused
elsewhere = fp.propagate(inputs.anchors[0], 200.0)
lines.append(f"the same anchor moved to season 200.0: {dict(elsewhere.propagation)}")
ok &= bool(np.array_equal(elsewhere.ln_N, inputs.anchors[0].ln_N))
record(10, "the hook is called for every anchor, is the identity, and records its absence",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 11. The refusals, each message quoted
# ---------------------------------------------------------------------------
CASE = WORK / "run"
(CASE / "inputs").mkdir(parents=True)
(CASE / "output").mkdir()
for key in ctl.RUN_INPUT_KINDS:
    shutil.copy(TRANSFER / f"inputs/lindal_transfer_{key}.nc",
                CASE / f"inputs/lindal_transfer_{key}.nc")
ANCHOR_REL = os.path.relpath(LINDAL / "lindal_refractivity.nc", CASE).replace("\\", "/")
BASE = NAMELIST.read_text(encoding="utf-8").replace(
    '"../../occul_data/lindal/lindal_refractivity.nc"', f'"{ANCHOR_REL}"')


def case(text, name="lindal_transfer"):
    """Write a namelist variant into the case directory and return its path."""
    path = CASE / f"{name}.toml"
    path.write_text(text, encoding="utf-8")
    return path


def load(text):
    """Read and load a namelist variant, returning the refusal message or None."""
    return refusal(lambda: ctl.load_run_inputs(ctl.read_run_namelist(case(text))))


def edit(old, new, text=None):
    """One substitution in the namelist, refusing to pass silently if it matched nothing.

    The namelist's `[inputs]` keys are column aligned, and a replacement written with single
    spaces matched nothing and left the case testing the unedited namelist, which reported no
    refusal for the right reason and the wrong one.
    """
    text = BASE if text is None else text
    count = text.count(old)
    if count != 1:
        raise AssertionError(f"the replacement {old!r} matches {count} times, not once")
    return text.replace(old, new)


# SPEC_04 v0.9 acceptance, decision N: these refusals and no others. The pipeline assumes its
# inputs came from the prior steps, so what is left is the schema, namelist typo guards,
# extrapolation (the interpolants', exercised at check 5), and named numerical failures (the
# later steps'). The coverage, point-composition, season and retrieval refusals of v0.7 are
# gone; the season is recorded and warned below.
cases = []
cases.append(("a missing [target]", load(
    edit("[target]\nlatitude_planetocentric_deg = 10.0\n", ""))))
cases.append(("an unknown [numerics] table", load(
    edit("[numerics.outer_loop]",
         "[numerics.vorticity]\nscheme = \"rk4\"\n[numerics.outer_loop]"))))
cases.append(("a scheme not implemented", load(
    edit("[numerics.isobar_tracing]\nscheme = \"rk4\"",
         "[numerics.isobar_tracing]\nscheme = \"euler\""))))
cases.append(("an anchor weight that is neither 1 nor 0", load(
    edit("weight = 1.0", "weight = 0.5"))))

ok = all(message is not None for _, message in cases)
record(11, f"{len(cases)} refusal cases, each message quoted (decision N: these and no others)",
       ok, "\n".join(f"{name}:\n    {message}" for name, message in cases))

# ---------------------------------------------------------------------------
# 12. A season that differs from the run's is recorded and warned, never refused
# ---------------------------------------------------------------------------
lines, ok = [], True
for key, name in (("wind", "lindal_transfer_wind"),
                  ("composition", "lindal_transfer_composition")):
    edited = CASE / f"inputs/{name}_season.nc"
    shutil.copy(CASE / f"inputs/{name}.nc", edited)
    with netCDF4.Dataset(edited, "a") as handle:
        if "season_absent_meaning" in handle.ncattrs():
            handle.delncattr("season_absent_meaning")
        handle.solar_longitude_deg = 200.0
        handle.solar_longitude_source = (
            "edited by the Step 0 acceptance to exercise the recorded season")
    variant = case(edit(f'{key:11s} = "inputs/{name}.nc"',
                        f'{key:11s} = "inputs/{name}_season.nc"'))
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        message = refusal(lambda: ctl.load_run_inputs(ctl.read_run_namelist(variant)))
        recorded = ({} if message is not None
                    else dict(ctl.load_run_inputs(ctl.read_run_namelist(variant)).seasons))
    warned = [str(w.message) for w in caught if "solar_longitude_deg" in str(w.message)]
    ok &= message is None and bool(warned) and recorded.get(key) == 200.0
    lines.append(f"the run's {key} at solar_longitude_deg = 200.0: refused {message}; recorded "
                 f"as {recorded.get(key)!r} beside the run's {recorded.get('run')!r}; warnings "
                 f"raised {len(warned)}")
    if warned:
        lines.append(f"    {warned[0]}")
lines.append(f"the unedited run's seasons, all equal: {dict(inputs.seasons)}")
ok &= dict(inputs.seasons) == {"run": 18.2, "wind": 18.2, "composition": "uniform",
                               "anchor:lindal": 18.2}
record(12, "a season that differs from the run's is recorded and warned, never refused",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 13. M > 1 and a validation anchor load
# ---------------------------------------------------------------------------
lines, ok = [], True
two = edit(
    'weight = 1.0\nmeasurement_uncertainty_scale = 1.0',
    'weight = 1.0\nmeasurement_uncertainty_scale = 1.0\n\n[[anchors]]\nslug   = "lindal"\n'
    f'path   = "{ANCHOR_REL}"\nweight = 0.0\nmeasurement_uncertainty_scale = 2.0')
loaded = ctl.load_run_inputs(ctl.read_run_namelist(case(two)))
ok &= len(loaded.anchors) == 2
roles = [(a.slug, a.weight, a.is_construction, a.measurement_uncertainty_scale) for a in
         loaded.anchors]
lines.append(f"two [[anchors]] entries pointing at the same file load: {len(loaded.anchors)} "
             f"anchors, (slug, weight, construction, scale) {roles}")
scaled = bool(np.allclose(loaded.anchors[1].sigma_ln_N_measurement,
                          2.0 * loaded.anchors[0].sigma_ln_N_measurement, rtol=0, atol=0))
ok &= scaled
lines.append(f"measurement_uncertainty_scale = 2.0 doubles sigma_ln_N_measurement: {scaled}")
lines.append("beyond the specification: the constant relative uncertainty of the Lindal "
             f"anchor is sigma_N / N = {loaded.anchors[0].sigma_ln_N_measurement[0]:.7f} at "
             "every level")

# A composition with one latitude node edited is refused nowhere: it is data, and whether the
# field varies with latitude is a property of its values (SPEC_04 Step 0 deliverable 2). The
# composition term it then makes nonzero is Step 3's.
varying = CASE / "inputs/lindal_transfer_composition_varying.nc"
shutil.copy(CASE / "inputs/lindal_transfer_composition.nc", varying)
with netCDF4.Dataset(varying, "a") as handle:
    node = int(np.flatnonzero(handle["latitude_planetocentric_deg"][:] == 45.0)[0])
    moved = float(handle["x_He"][0, node]) + 1.0e-3
    handle["x_He"][:, node] = np.asarray(handle["x_He"][:, node]) + 1.0e-3
    handle["x_H2"][:, node] = np.asarray(handle["x_H2"][:, node]) - 1.0e-3
message = refusal(lambda: ctl.load_run_inputs(ctl.read_run_namelist(case(edit(
    'composition = "inputs/lindal_transfer_composition.nc"',
    'composition = "inputs/lindal_transfer_composition_varying.nc"')))))
ok &= message is None
lines.append(f"a composition with the node at 45 deg edited (x_He there raised by 1.0e-3 to "
             f"{moved:.6f}, x_H2 lowered to keep the sum) loads: refused {message}")
record(13, "M > 1 loads, and a validation anchor loads and is marked as one", ok,
       "\n".join(lines))

# ---------------------------------------------------------------------------
# 14. The closure namelist still loads, and the SPEC_03 Step 3 suite still passes
# ---------------------------------------------------------------------------
closure_namelist = ctl.read_run_namelist(CLOSURE / "lindal_closure.toml")
closure_inputs = ctl.load_run_inputs(closure_namelist)
lines = [f"the closure namelist loads: mode {closure_namelist.mode!r}, "
         f"{len(closure_namelist.anchors)} anchor, gauge level "
         f"{closure_inputs.gauge_level_index}, closure comparisons "
         f"{[(c.kind, c.identical) for c in closure_inputs.closure]}",
         f"transfer-only fields are empty in closure mode: target "
         f"{closure_namelist.target_latitude_deg}, grid {dict(closure_namelist.grid)}, "
         f"numerics {dict(closure_namelist.numerics)}, estimation "
         f"{dict(closure_namelist.estimation)}"]
ok = (closure_namelist.mode == "closure"
      and all(c.identical for c in closure_inputs.closure)
      and closure_namelist.target_latitude_deg is None)
completed = subprocess.run(
    [sys.executable, "tests/step03_3/accept_step03_3.py"],
    capture_output=True, text=True, encoding="utf-8", errors="replace")
tail = [line for line in completed.stdout.splitlines() if line.startswith(("[PASS", "[FAIL"))]
passed = sum(1 for line in tail if line.startswith("[PASS"))
lines.append(f"tests/step03_3/accept_step03_3.py: {passed} of {len(tail)} pass "
             f"(exit {completed.returncode}); the thirteen closure refusals are its check 12")
ok &= (passed == len(tail) and len(tail) > 0)
(HERE / "step03_3_rerun.txt").write_text(completed.stdout + completed.stderr, encoding="utf-8")
record(14, "the closure namelist and the SPEC_03 Step 3 acceptance are unchanged", ok,
       "\n".join(lines))

# ---------------------------------------------------------------------------
print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass.")
for number, description, _, _ in failed:
    print(f"  FAIL {number}. {description}")
sys.exit(1 if failed else 0)
