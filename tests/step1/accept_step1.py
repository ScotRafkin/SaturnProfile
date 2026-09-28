"""Acceptance checks for SPEC_01 Step 1. Prints the actual result of each check.

The five checks are the ones listed under Acceptance in SPEC_01 v0.3 Step 1. Check 6 is
beyond the specification and is run because Step 2 writes a group bearing raw bundle.
"""

import shutil
import subprocess
import sys
from pathlib import Path

import netCDF4
import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError

OUT = Path("reports/step1/products")
OUT.mkdir(parents=True, exist_ok=True)

results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    mark = "PASS" if passed else "FAIL"
    print(f"[{mark}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def author_globals(title):
    return {
        "title": title,
        "profile_or_run": "step1_acceptance",
        "role": "reduction",
        "source": "SPEC_01 Step 1 acceptance check, not a physical dataset",
        # SPEC_00 v0.17 section 5 (SPEC_03 Step 3): every file carries its epoch and season.
        "epoch": "2026-09-10",
        "season_absent_meaning": "uniform",
    }


# ---------------------------------------------------------------------------
# 1. Write a two level raw file with one variable and read it back
# ---------------------------------------------------------------------------
path = OUT / "acceptance_raw.nc"
ds = xr.Dataset(
    {"example_value_Pa": ("level", np.array([1.0, 2.0], dtype="float64"))},
    coords={"level": ("level", np.array([0, 1], dtype="int32"))},
    attrs=author_globals("Two level raw file, Step 1 acceptance"),
)
cio.write(path, ds, "raw", created_by="accept_step1")
# `raw` is a grouped kind, so `read` returns a DataTree by the rule of SPEC_00 section 3.1.
tree_back = cio.read(path, "raw")
back = tree_back.dataset

expected_types = {
    "Conventions": str,
    "casspian_kind": str,
    "casspian_schema_version": (int, np.integer),
    "title": str,
    "profile_or_run": str,
    "role": str,
    "source": str,
    "created_by": str,
    "created_at": str,
    "casspian_git_commit": str,
    "codata_release": str,
    "history": str,
}
missing, wrong = [], []
for name, want in expected_types.items():
    if name not in back.attrs:
        missing.append(name)
    elif not isinstance(back.attrs[name], want):
        wrong.append(f"{name} is {type(back.attrs[name]).__name__}")
values_ok = (
    back.sizes["level"] == 2
    and list(back.data_vars) == ["example_value_Pa"]
    and np.array_equal(back["example_value_Pa"].values, np.array([1.0, 2.0]))
)
detail = "\n".join(
    [f"{k} = {back.attrs[k]!r}" for k in expected_types if k in back.attrs]
    + [f"level size = {back.sizes['level']}, values = {back['example_value_Pa'].values}"]
)
record(1, "raw file written and read back, every section 5 global present and typed",
       not missing and not wrong and values_ok,
       detail if not (missing or wrong) else f"missing={missing} wrong={wrong}")
tree_back.close()

# ---------------------------------------------------------------------------
# 2. read refuses it when asked for kind thermo
# ---------------------------------------------------------------------------
try:
    cio.read(path, "thermo")
    record(2, "read refuses the raw file when asked for kind thermo", False, "it did not refuse")
except CasspianSchemaError as exc:
    record(2, "read refuses the raw file when asked for kind thermo", True, str(exc))

# ---------------------------------------------------------------------------
# 3. read refuses a file whose casspian_kind was edited to wind
# ---------------------------------------------------------------------------
edited = OUT / "acceptance_raw_edited_to_wind.nc"
shutil.copy(path, edited)
with netCDF4.Dataset(edited, "a") as handle:
    handle.setncattr("casspian_kind", "wind")
try:
    cio.read(edited, "wind")
    record(3, "read refuses a file relabeled wind that lacks the wind variables", False,
           "it did not refuse")
except CasspianSchemaError as exc:
    record(3, "read refuses a file relabeled wind that lacks the wind variables", True, str(exc))

# ---------------------------------------------------------------------------
# 4. sha256 matches sha256sum on the command line
# ---------------------------------------------------------------------------
ours = cio.sha256(path)
proc = subprocess.run(["sha256sum", str(path)], capture_output=True, text=True)
# GNU coreutils escapes a line whose filename contains a backslash by prefixing it with one,
# which every Windows path does. The prefix is not part of the digest.
theirs = proc.stdout.split()[0].lstrip("\\") if proc.returncode == 0 else "sha256sum unavailable"
record(4, "sha256 matches sha256sum on the command line", ours == theirs,
       f"lib.io.sha256 = {ours}\nsha256sum    = {theirs}")

# ---------------------------------------------------------------------------
# 5. a dataset with a variable named latitude is refused by write
# ---------------------------------------------------------------------------
bad = xr.Dataset(
    {"latitude": ("level", np.array([1.0, 2.0], dtype="float64"))},
    attrs=author_globals("Bare latitude, must be refused"),
)
try:
    cio.write(OUT / "must_not_exist.nc", bad, "raw", created_by="accept_step1")
    record(5, "write refuses a dataset with a variable named latitude", False, "it did not refuse")
except CasspianSchemaError as exc:
    record(5, "write refuses a dataset with a variable named latitude", True, str(exc))

# ---------------------------------------------------------------------------
# 6. Beyond the specification: groups survive a write and read round trip
# ---------------------------------------------------------------------------
grouped = OUT / "acceptance_raw_grouped.nc"
root = xr.Dataset(attrs=author_globals("Group bearing raw file, beyond Step 1 acceptance"))
groups = {
    "table1": xr.Dataset(
        {"pressure_Pa": ("level", np.array([20.0, 100000.0], dtype="float64"))},
        coords={"level": ("level", np.array([0, 1], dtype="int32"))},
    ),
    "scalars": xr.Dataset(attrs={"note": "a scalars group"}),
    "scalars/latitude": xr.Dataset(attrs={"planetographic_deg": 36.3}),
}
cio.write(grouped, root, "raw", groups=groups, created_by="accept_step1")
tree = cio.read(grouped, "raw")
found = cio.group_paths(grouped)
nested = float(tree["scalars/latitude"].attrs["planetographic_deg"])
record(6, "groups round trip, including a nested group (beyond the specification)",
       sorted(found) == ["scalars", "scalars/latitude", "table1"] and nested == 36.3,
       f"groups in file = {sorted(found)}\nscalars/latitude planetographic_deg = {nested}")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
