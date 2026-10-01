"""Acceptance checks for SPEC_02 v0.2 Step 1, `lib.control` manifest reading and input loading.

Every check reports its measured value. The refusal cases run on copies of the manifest and its
six inputs under this directory, because SPEC_00 section 8 refuses a manifest that points outside
its own directory, so a tampered manifest has to sit beside the files it names. Nothing under
`occul_data/` is written.
"""

import importlib
import re
import shutil
import sys
from pathlib import Path

import netCDF4

from casspian.lib import control as ctl
from casspian.lib.control import ControlFileError

HERE = Path("reports/step02_1")
HERE.mkdir(parents=True, exist_ok=True)
SOURCE = Path("occul_data/lindal")
MANIFEST = SOURCE / "lindal_reduction.toml"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def copy_case(name):
    """A fresh directory holding the manifest and the six inputs it names."""
    case = HERE / name
    if case.exists():
        shutil.rmtree(case)
    case.mkdir(parents=True)
    shutil.copy2(MANIFEST, case / MANIFEST.name)
    for path in ctl.read_reduction_manifest(MANIFEST).inputs.values():
        shutil.copy2(path, case / path.name)
    return case


def edit_manifest(case, pattern, replacement):
    target = case / MANIFEST.name
    text = target.read_text(encoding="utf-8")
    edited, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"the acceptance edit {pattern!r} matched {count} times")
    target.write_bytes(edited.encode("utf-8"))
    return target


def refusal(action):
    try:
        action()
    except ControlFileError as exc:
        return True, str(exc)
    return False, "not refused"


# ---------------------------------------------------------------------------
# 1. The real manifest loads, with the hashes REPORT_01_step9 section 6 recorded
# ---------------------------------------------------------------------------
report = Path("reports/REPORT_01_step9.md").read_text(encoding="utf-8")
recorded = dict(re.findall(r"\| `([\w/]+\.(?:nc|toml))` \| \w+ \| `([0-9a-f]{64})` \|", report))
# SPEC_01 v0.19 rebuilt the composition; the addendum to REPORT_01_step8 records its new hash,
# which supersedes the Step 9 table for that file. A later row overrides an earlier one.
# SPEC_03 v0.5 Step 0 rebuilt the whole chain; REPORT_03_step0 records the new hashes.
# SPEC_03 Step 3 rebuilt the whole chain again; REPORT_03_step3 records those hashes.
# SPEC_04 Step 0 rebuilt it once more, for kind W in three parts; REPORT_04_step0 records those,
# provisionally until its sweep record replaces them (handoff of 17 September section 7).
for line in [text for name in ("reports/REPORT_01_step8.md", "reports/REPORT_02_step6.md",
                               "reports/REPORT_03_step0.md", "reports/REPORT_03_step3.md",
                               "reports/REPORT_04_step0.md",
                               "reports/REPORT_05_step0.md")
             if Path(name).exists()
             for text in Path(name).read_text(encoding="utf-8").splitlines()]:
    if line.startswith("| `lindal_") and line.count("`") == 4:
        cells = [c.strip().strip("`") for c in line.strip().strip("|").split("|")]
        if len(cells) == 3 and len(cells[2]) == 64:
            recorded[cells[0]] = cells[2]

manifest = ctl.read_reduction_manifest(MANIFEST)
inputs = ctl.load_reduction_inputs(manifest)
lines, all_match = [], True
for key in ctl.INPUT_KINDS:
    name = manifest.inputs[key].name
    match = recorded.get(name) == inputs.sha256[key]
    all_match = all_match and match
    lines.append(f"  {key:12s} {name:24s} sha256 {inputs.sha256[key][:16]}...  "
                 f"{'matches' if match else 'DIFFERS from'} REPORT_01_step9  "
                 f"commit {inputs.commits[key][:12]}")
manifest_match = recorded.get(MANIFEST.name) == manifest.sha256
record(1, "lindal_reduction.toml loads with all six inputs and the recorded hashes",
       all_match and manifest_match and len(inputs.sha256) == 6,
       "\n".join(lines)
       + f"\n  manifest     {MANIFEST.name:24s} sha256 {manifest.sha256[:16]}...  "
         f"{'matches' if manifest_match else 'DIFFERS from'} REPORT_01_step9\n"
       f"slug {manifest.slug!r}; anchor isobar {manifest.anchor_isobar_Pa:g} Pa; "
       f"anchor rule {manifest.anchor_rule!r} on the {manifest.anchor_surface_Pa:g} Pa surface "
       f"by {manifest.anchor_quantity!r}\n"
       f"fixed_point_tolerance_deg = {manifest.fixed_point_tolerance_deg!r}, returned in "
       f"degrees as written; the record offers no radian form "
       f"({not any('rad' in f for f in manifest.__dataclass_fields__)}), because SPEC_02 "
       f"section 0 converts once, at the refrac boundary\n"
       f"product {manifest.product.name!r} (to be written by refrac); diagnostics "
       f"{dict(manifest.diagnostics)}, the section 7.1 default, since this manifest carries "
       f"no optional section (v0.12 removed [sensitivity])\n"
       f"composition loaded as {type(inputs.composition).__name__}, the other five as "
       f"{sorted({type(getattr(inputs, k)).__name__ for k in ctl.INPUT_KINDS if k != 'composition'})}")

# ---------------------------------------------------------------------------
# 2. A manifest naming a -dirty product is refused, with the file named
# ---------------------------------------------------------------------------
case = copy_case("case_dirty")
with netCDF4.Dataset(case / "lindal_wind.nc", "a") as handle:
    handle.setncattr("casspian_git_commit", handle.getncattr("casspian_git_commit") + "-dirty")
refused, message = refusal(lambda: ctl.load_reduction_inputs(
    ctl.read_reduction_manifest(case / MANIFEST.name)))
record(2, "a manifest naming a -dirty product is refused with the file named",
       refused and "lindal_wind.nc" in message and "-dirty" in message,
       f"lindal_wind.nc in the copy carries casspian_git_commit ending in -dirty\n"
       f"refused: {refused}\n  {message}")

# ---------------------------------------------------------------------------
# 3. pressure_Pa = 9000 is refused on both counts
# ---------------------------------------------------------------------------
case = copy_case("case_9000")
edit_manifest(case, r"^(\[anchor_isobar\]\s*\npressure_Pa\s*=\s*)\S+", r"\g<1>9000")
refused, message = refusal(lambda: ctl.load_reduction_inputs(
    ctl.read_reduction_manifest(case / MANIFEST.name)))
both = "not a surface of the geodesy file" in message and "not a tabulated level" in message
record(3, "[anchor_isobar] pressure_Pa = 9000 is refused as neither a geodesy surface nor a "
       "tabulated level", refused and both,
       f"refused: {refused}; the message names both failures: {both}\n  {message}\n"
       "The manifest itself parses: 9000 is a positive pressure. The refusal comes from "
       "checking it against the inputs, which is why it is a load time rule and not a parse "
       "time one.")

# ---------------------------------------------------------------------------
# 4. An extra polar_radius_m under [geoid] is refused as a physical value
# ---------------------------------------------------------------------------
case = copy_case("case_physical")
edit_manifest(case, r"^\[geoid\]\s*$", "[geoid]\npolar_radius_m = 54438000")
refused, message = refusal(lambda: ctl.read_reduction_manifest(case / MANIFEST.name))
record(4, "an extra polar_radius_m = 54438000 under [geoid] is refused as a physical value",
       refused and "physical value" in message,
       f"refused at parse time: {refused}\n  {message}")

# ---------------------------------------------------------------------------
# 5. Beyond the specification: the other parse time refusals
# ---------------------------------------------------------------------------
case = copy_case("case_unknown")
edit_manifest(case, r"^\[geoid\]\s*$", '[geoid]\nnote_to_self = "remember this"')
unknown_refused, unknown_message = refusal(
    lambda: ctl.read_reduction_manifest(case / MANIFEST.name))
case = copy_case("case_outside")
edit_manifest(case, r'^(wind\s*=\s*)"lindal_wind\.nc"', r'\g<1>"../lindal_wind.nc"')
outside_refused, outside_message = refusal(
    lambda: ctl.read_reduction_manifest(case / MANIFEST.name))
case = copy_case("case_rule")
edit_manifest(case, r'^(anchor_rule\s*=\s*)"[^"]+"', r'\g<1>"latitude"')
rule_refused, rule_message = refusal(lambda: ctl.read_reduction_manifest(case / MANIFEST.name))
record(5, "beyond the specification: an unknown key, a path outside the directory, and an "
       "inexpressible anchor rule are each refused by name",
       unknown_refused and "unknown key" in unknown_message and "physical" not in unknown_message
       and outside_refused and "outside" in outside_message
       and rule_refused and "latitude" in rule_message,
       f"unknown text key: {unknown_refused}, and named as unknown rather than physical\n"
       f"  {unknown_message}\n"
       f"wind = \"../lindal_wind.nc\": {outside_refused}\n  {outside_message}\n"
       f"anchor_rule = \"latitude\": {rule_refused}\n  {rule_message}")

# ---------------------------------------------------------------------------
# 6. Beyond the specification: the old module is absorbed, not shimmed
# ---------------------------------------------------------------------------
try:
    importlib.import_module("casspian.tools.gravity.control")
    gone = False
except ModuleNotFoundError:
    gone = True
# The gravity tool's own key table, so the check follows the section vocabulary as it grows
# (SPEC_01 v0.26 role, v0.28 epoch) rather than a copy of it.
from casspian.tools.gravity.build_gravity import SECTION_KEYS as GRAVITY_KEYS  # noqa: E402

section = ctl.load_section(SOURCE / "lindal_build.toml", "gravity", GRAVITY_KEYS,
                           path_keys=("source", "output"))
record(6, "beyond the specification: tools/gravity/control.py is absorbed into lib.control, "
       "with its named-key path resolution kept",
       gone and section["source"].is_absolute() and section["source"].exists(),
       f"casspian.tools.gravity.control is no longer importable: {gone}\n"
       f"load_section on lindal_build.toml [gravity] resolves source to an absolute path that "
       f"exists: {section['source'].is_absolute() and section['source'].exists()}\n"
       f"  {section['source'].as_posix()}\n"
       "The five tools that imported the old module now import lib.control; their suites are "
       "rerun in the regression.")

print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
