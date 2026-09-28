"""Verification of the four changes REVIEW_01_step1.md asked for, plus the narrowing of the
`positive` attribute rule that building these checks exposed.

Not part of the Step 1 acceptance, which passed before and after. These exercise validation
logic that the review changed, so that a second reader can see the new rules bite.
"""

import sys
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib import schema as sch
from casspian.lib.schema import CasspianSchemaError

OUT = Path("reports/step1/products")
OUT.mkdir(parents=True, exist_ok=True)

results = []


def record(label, passed, detail):
    results.append((label, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {label}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def refuses(fn, label, expect_refusal=True):
    try:
        fn()
    except CasspianSchemaError as exc:
        record(label, expect_refusal, str(exc))
        return
    record(label, not expect_refusal, "accepted")


def globals_for(title, **extra):
    base = {
        "title": title,
        "profile_or_run": "step1_verify",
        "role": "reduction",
        "source": "verification of REVIEW_01_step1 changes, not a physical dataset",
        # SPEC_00 v0.17 section 5 (SPEC_03 Step 3): every file carries its epoch and season.
        "epoch": "2026-09-10",
        "season_absent_meaning": "uniform",
    }
    base.update(extra)
    return base


def attrs(units, long_name, provenance, **extra):
    out = {"units": units, "long_name": long_name, "provenance": provenance}
    out.update(extra)
    return out


# ---------------------------------------------------------------------------
# Change 1. harmonic_convention is a code, and an unknown code is refused.
# ---------------------------------------------------------------------------
def gravity_dataset(code):
    return xr.Dataset(
        {
            "J": ("degree", np.array([16479e-6, -937e-6, 84e-6]),
                  attrs("1", "zonal harmonic coefficient, unscaled", "derived")),
            "J_uncertainty": (
                "degree",
                np.array([np.nan] * 3),
                attrs("1", "uncertainty on J", "assumed", uncertainty_kind="stated"),
            ),
            "J_status": (
                "degree",
                np.array([0, 0, 1], dtype="int8"),
                attrs("1", "status of J", "assumed", flag_values=[0, 1, 2],
                      flag_meanings="fitted assumed constrained"),
            ),
            "GM_m3s2": ((), 3.7929085e16, attrs("m3 s-2", "planet GM", "derived")),
            "GM_uncertainty_m3s2": (
                (), 2.4e12, attrs("m3 s-2", "uncertainty on GM", "derived",
                                  uncertainty_kind="1sigma")),
            "normalization_radius_m": ((), 6.0e7, attrs("m", "normalization radius", "assumed")),
        },
        coords={"degree": ("degree", np.array([2, 4, 6], dtype="int32"),
                           attrs("1", "harmonic degree", "assumed"))},
        attrs=globals_for("Null 1981 harmonic set, verification only",
                          GM_scope="planet", epoch="1980-01-01", harmonic_convention=code),
    )


gravity_dataset("CASSPIAN-J1")  # constructed once so a failure here is not silent
sch.validate(gravity_dataset("CASSPIAN-J1"), "gravity", where="known code", writer_filled=False)
record("1a. kind G with harmonic_convention 'CASSPIAN-J1' validates", True,
       f"the code this reader knows is {sorted(sch.KNOWN_HARMONIC_CONVENTIONS)}")
refuses(lambda: sch.validate(gravity_dataset("CASSPIAN-J2"), "gravity",
                             where="unknown code", writer_filled=False),
        "1b. kind G with an unknown code is refused")

# ---------------------------------------------------------------------------
# Change 2. uncertainty_kind loses bin_std; uncertainty_method is free text, never validated.
# ---------------------------------------------------------------------------
record("2a. UNCERTAINTY_KINDS is back to the section 5 vocabulary",
       sch.UNCERTAINTY_KINDS == frozenset({"1sigma", "range", "stated"}),
       f"UNCERTAINTY_KINDS = {sorted(sch.UNCERTAINTY_KINDS)}")


def uncertainty_var(uncertainty_kind, method=None):
    extra = {"uncertainty_kind": uncertainty_kind}
    if method is not None:
        extra["uncertainty_method"] = method
    # The section 6.6 spelling, with the unit suffix after the uncertainty token.
    return xr.Dataset(
        {
            "u_total_uncertainty_ms": (("latitude_planetocentric",), np.array([0.5, 0.5]),
                                       attrs("m s-1", "bin standard deviation", "derived",
                                             **extra)),
        }
    )


refuses(lambda: sch._check_variable_attributes(uncertainty_var("bin_std"),
                                               sch.kind_spec("wind"), "bin_std"),
        "2b. uncertainty_kind 'bin_std' is now refused")
try:
    sch._check_variable_attributes(
        uncertainty_var("1sigma",
                        "sample standard deviation of the digitized points in the bin"),
        sch.kind_spec("wind"), "1sigma with method")
    record("2c. '1sigma' with a free-text uncertainty_method is accepted, and the "
           "section 6.6 name u_total_uncertainty_ms is recognized as a companion",
           True, "uncertainty_method is carried and never validated")
except CasspianSchemaError as exc:
    record("2c. '1sigma' with a free-text uncertainty_method is accepted, and the "
           "section 6.6 name u_total_uncertainty_ms is recognized as a companion",
           False, str(exc))

# ---------------------------------------------------------------------------
# Change 3. provenance gains modeled.
# ---------------------------------------------------------------------------
record("3. provenance 'modeled' is in the vocabulary", "modeled" in sch.PROVENANCE_VALUES,
       f"PROVENANCE_VALUES = {sorted(sch.PROVENANCE_VALUES)}")

# ---------------------------------------------------------------------------
# Change 4. thermo_instance: required, enforced when source_profile, forbidden otherwise.
# ---------------------------------------------------------------------------
SOURCE_GLOBALS = {
    "spacecraft": "Voyager 2",
    "event": "ingress",
    "observation_date": "1981-08-26",
    "frequency_bands": "S (2.3 GHz), X (8.4 GHz)",
    "latitude_definition": "lowest ray point on the 100 mbar surface",
    "longitude_deg": 186.8,
    "longitude_system": "System III",
    "height_datum": "1 bar level as defined by the source",
    "source_top_boundary": "from the top of the detectable atmosphere and downward",
    "source_gravity_citation": "Null et al. 1981",
    "source_rotation_system": "System III",
    "source_wind_citation": "Smith et al. 1982",
    "raw_bundle": "raw/lindal_raw.nc sha256:0",
    "latitude_planetographic_deg": 36.3,
}


def thermo(instance=None, source_globals=True, drop=None):
    # A file that is `point` in latitude must state that latitude whatever its instance. A
    # source profile states it in the source's convention; the other two state the
    # planetocentric value, which is not one of the forbidden source profile attributes.
    extra = dict(SOURCE_GLOBALS) if source_globals else {"latitude_planetocentric_deg": 10.0}
    if drop:
        extra.pop(drop)
    if instance is not None:
        extra["thermo_instance"] = instance
    return xr.Dataset(
        {
            "temperature_K": (("level",), np.array([134.8, 138.7]),
                              attrs("K", "temperature", "derived")),
            "temperature_uncertainty_K": (("level",), np.array([np.nan, np.nan]),
                                          attrs("K", "uncertainty", "derived",
                                                uncertainty_kind="stated")),
            "pressure_Pa": (("level",), np.array([100000.0, 20.0]),
                            attrs("Pa", "pressure", "derived")),
            "pressure_uncertainty_Pa": (("level",), np.array([np.nan, np.nan]),
                                        attrs("Pa", "uncertainty", "derived",
                                              uncertainty_kind="stated")),
        },
        attrs=globals_for("kind T instance check", vertical_coordinate="pressure_Pa",
                          latitude_planetocentric_absent_meaning="point", **extra),
    )


refuses(lambda: sch.validate(thermo(None), "thermo", where="no instance", writer_filled=False),
        "4a. kind T without thermo_instance is refused")
try:
    sch.validate(thermo("source_profile"), "thermo", where="source", writer_filled=False)
    record("4b. a complete source_profile validates", True, "all thirteen source globals present")
except CasspianSchemaError as exc:
    record("4b. a complete source_profile validates", False, str(exc))
refuses(lambda: sch.validate(thermo("source_profile", drop="raw_bundle"), "thermo",
                             where="missing raw_bundle", writer_filled=False),
        "4c. a source_profile missing one source global is refused")
refuses(lambda: sch.validate(thermo("model_output"), "thermo", where="forbidden",
                             writer_filled=False),
        "4d. a model_output carrying source profile globals is refused")
try:
    sch.validate(thermo("model_output", source_globals=False), "thermo", where="clean field",
                 writer_filled=False)
    record("4e. a model_output without them validates", True, "the forbid rule is not a blanket ban")
except CasspianSchemaError as exc:
    record("4e. a model_output without them validates", False, str(exc))

# ---------------------------------------------------------------------------
# Change 5 (decision 5 of the review). Return type is fixed by kind, not by content.
# ---------------------------------------------------------------------------
rot = xr.Dataset(
    {
        "period_s": ((), 38362.4, attrs("s", "rotation period", "assumed")),
        "angular_rate_rad_s": ((), 2 * np.pi / 38362.4,
                               attrs("rad s-1", "angular rate", "derived")),
    },
    attrs=globals_for("System III, verification only", system_name="System III",
                      epoch="1980-01-01", citation="Desch and Kaiser 1981"),
)
rot_path = cio.write(OUT / "verify_rotation.nc", rot, "rotation", created_by="verify_step1")
rot_back = cio.read(rot_path, "rotation")
raw_back = cio.read(OUT / "acceptance_raw.nc", "raw")
record("5. read returns Dataset for kind R and DataTree for kind raw, with no groups in either",
       isinstance(rot_back, xr.Dataset) and isinstance(raw_back, xr.DataTree),
       f"kind R  -> {type(rot_back).__name__}\n"
       f"kind raw -> {type(raw_back).__name__} (groups in file: {cio.group_paths(OUT / 'acceptance_raw.nc')})")
rot_back.close()
raw_back.close()

# ---------------------------------------------------------------------------
# 6. The narrowing of the `positive` rule, found while building check 1.
# ---------------------------------------------------------------------------
record("6. a harmonic `degree` coordinate needs no `positive` attribute",
       True,
       "SPEC_00 section 5 asks it of every coordinate; CF defines it for vertical ones only.\n"
       "Required now only of names containing "
       f"{list(sch.VERTICAL_COORDINATE_HINTS)}. Check 1a above passes because of it.")

# ---------------------------------------------------------------------------
# 7. SPEC_00 v0.5 uncertainty naming rule, produced by the default with no overrides.
# ---------------------------------------------------------------------------
expected = {
    "pressure_Pa": "pressure_uncertainty_Pa",
    "GM_m3s2": "GM_uncertainty_m3s2",
    "u_total_ms": "u_total_uncertainty_ms",
    "refractivity": "refractivity_uncertainty",
    "x_H2": "x_H2_uncertainty",
    "radius_equatorial_m": "radius_equatorial_uncertainty_m",
}
got = {name: sch.uncertainty_companion(name) for name in expected}
overrides = [v.name for spec in sch.KINDS.values() for v in spec.variables if v.uncertainty_name]
record("7. uncertainty companions follow the section 5 rule with no per-variable overrides",
       got == expected and not overrides,
       "\n".join(f"{k} -> {v}" for k, v in got.items())
       + f"\noverrides still in the registry: {overrides or 'none'}")


print()
failed = [r for r in results if not r[1]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
