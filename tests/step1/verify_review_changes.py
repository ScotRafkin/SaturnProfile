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
sch.validate(gravity_dataset("CASSPIAN-J1"), "gravity", where="known code")
record("1a. kind G with harmonic_convention 'CASSPIAN-J1' validates", True,
       f"the code this reader knows is {sorted(sch.KNOWN_HARMONIC_CONVENTIONS)}")
refuses(lambda: sch.validate(gravity_dataset("CASSPIAN-J2"), "gravity",
                             where="unknown code"),
        "1b. kind G with an unknown code is refused")

# Changes 2 to 4 and the `positive` narrowing (old checks 2a to 4e and 6) are retired by
# SPEC_R1: uncertainty_kind, provenance, thermo_instance, the source profile globals and the
# `positive` attribute are no longer required or checked.

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
