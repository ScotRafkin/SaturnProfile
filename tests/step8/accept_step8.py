"""Acceptance checks for SPEC_01 Step 8, the Lindal composition and kind C.

Every check reports its measured value.
"""

import math
import shutil
import sys
from pathlib import Path

import netCDF4
import numpy as np

from casspian.lib import composition as comp
from casspian.lib import io as cio
from casspian.lib.constants import (
    ATOMIC_MASS_CONSTANT,
    AVOGADRO_CONSTANT,
    LOSCHMIDT_CONSTANT,
)
from casspian.lib.schema import CasspianSchemaError

HERE = Path("reports/step8")
HERE.mkdir(parents=True, exist_ok=True)
PRODUCT = Path("occul_data/lindal/lindal_composition.nc")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


tree = cio.read(PRODUCT, "composition")
field = tree.dataset
# SPEC_04 decision O: kind C is one structure for every use, a field on (level, latitude), and
# the one-latitude form is retired. This suite checks the source's column, so it takes one, by
# the rule of decision L that every reader uses. The Lindal composition is uniform in latitude
# by construction, which is asserted here so that taking a single column hides nothing: the
# reduction takes the column at phi_c and gets these same values.
_varying = [name for name, v in field.data_vars.items()
            if comp.LATITUDE_DIM in v.dims and not np.all(v.values == v.values[:, :1])]
assert not _varying, f"the composition varies with latitude in {_varying}; this suite reads "\
                     "the source's column and assumes the Lindal field is uniform"
root = comp.column_at(tree, float(field[comp.LATITUDE_NAME].values[0])).dataset
species_group = tree["species"].dataset

pressure = root["pressure_Pa"].values
mbar = pressure / 100.0
x_h2 = root["x_H2"].values
x_he = root["x_He"].values
x_nh3 = root["x_NH3"].values
names = [str(s) for s in species_group["species_name"].values]
molar = species_group["molar_mass_kg_mol"].values
script_R = species_group["refractivity_per_molecule_m3"].values

AMU = ATOMIC_MASS_CONSTANT * AVOGADRO_CONSTANT


def mean_molar_mass(i):
    return (x_h2[i] * molar[0] + x_he[i] * molar[1] + x_nh3[i] * molar[2]) / AMU


# ---------------------------------------------------------------------------
# 1. Mole fractions sum to one
# ---------------------------------------------------------------------------
total = x_h2 + x_he + x_nh3
worst = float(np.max(np.abs(total - 1.0)))
record(1, "the mole fractions sum to one at every level to 1e-12", worst <= 1e-12,
       f"{pressure.size} levels; max |sum - 1| = {worst:.3e}\n"
       "The sum is exact by construction rather than by normalization: ammonia is assigned "
       "first and the remainder is split, so the three fractions are x, s(1-x) and (1-s)(1-x), "
       "which sum to one identically. What is measured above is floating point round-off.")

# ---------------------------------------------------------------------------
# 2. The ammonia clamp
# ---------------------------------------------------------------------------
flags = root["nh3_provenance"].values
meanings = root["nh3_provenance"].attrs["flag_meanings"].split()
tabulated = np.flatnonzero(flags == 0)
lowest_tabulated = float(mbar[tabulated].min())
aloft = mbar < lowest_tabulated
first_aloft = float(mbar[aloft].max())
at_79433 = float(x_nh3[int(np.argmin(np.abs(mbar - 794.33)))])
flag_79433 = meanings[int(flags[int(np.argmin(np.abs(mbar - 794.33)))])]
at_83176 = float(x_nh3[int(np.argmin(np.abs(mbar - 831.76)))]) * 1e6
flag_83176 = meanings[int(flags[int(np.argmin(np.abs(mbar - 831.76)))])]
all_zero = bool(np.all(x_nh3[aloft] == 0.0))
all_assumed = bool(np.all(flags[aloft] == 2))
record(2, "ammonia is exactly zero above the highest tabulated level, flagged assumed",
       all_zero and all_assumed and at_79433 == 0.0 and abs(at_83176 - 2.6) < 1e-9,
       f"highest tabulated level {lowest_tabulated:.2f} mbar; the first grid level above it "
       f"is {first_aloft:.2f} mbar\n"
       f"at 794.33 mbar: {at_79433:.1e} mole fraction, exactly zero, flagged {flag_79433!r}\n"
       f"at 831.76 mbar: {at_83176:.3f} ppm, flagged {flag_83176!r}\n"
       f"{int(aloft.sum())} levels above the tabulated range, all exactly zero and all "
       f"flagged assumed: {all_zero and all_assumed}\n"
       f"provenance codes present: "
       f"{[meanings[c] for c in sorted(set(flags.tolist()))]}\n"
       f"The v0.17 rule assigns zero rather than extrapolating upward and clamping, so the "
       f"clamp crossing that put 0.213 ppm here under v0.5 no longer exists. The reason is in "
       f"the file: {root.attrs['nh3_aloft_assumption'][:90]}...")

# ---------------------------------------------------------------------------
# 3. The interior gap and the deep extrapolation
# ---------------------------------------------------------------------------
gap = float(x_nh3[int(np.argmin(np.abs(mbar - 1047.13)))]) * 1e6
deep = float(x_nh3[int(np.argmin(np.abs(mbar - 1298.48)))]) * 1e6
record(3, "the interior gap at 1047.13 mbar and the extrapolated value at 1298.48 mbar",
       abs(gap - 15.9) <= 0.1,
       f"1047.13 mbar: {gap:.3f} ppm (spec 15.9 +- 0.1), interior gap filled by interpolation\n"
       f"1298.48 mbar: {deep:.3f} ppm, the end gap filled by extrapolation "
       f"(Phase 1 found 79.3)\n"
       f"the nine tabulated values run {x_nh3[tabulated].min() * 1e6:.1f} to "
       f"{x_nh3[tabulated].max() * 1e6:.1f} ppm between {mbar[tabulated].min():.2f} and "
       f"{mbar[tabulated].max():.2f} mbar")

# ---------------------------------------------------------------------------
# 4. Mean molar mass
# ---------------------------------------------------------------------------
zero_nh3 = int(np.flatnonzero(x_nh3 == 0.0)[0])
deepest = int(np.argmax(pressure))
m_dry = mean_molar_mass(zero_nh3)
m_deep = mean_molar_mass(deepest)
record(4, "mean molar mass with ammonia zero, and at the deepest level",
       abs(m_dry - 2.1351) <= 5e-4,
       f"with NH3 = 0: {m_dry:.6f} amu (Phase 1 check 2.1351)\n"
       f"at the deepest level, {mbar[deepest]:.2f} mbar with NH3 = "
       f"{x_nh3[deepest] * 1e6:.1f} ppm: {m_deep:.6f} amu\n"
       f"the ammonia raises the mean molar mass by {m_deep - m_dry:.6f} amu, "
       f"{100 * (m_deep - m_dry) / m_dry:.4f} percent, which is the whole of its effect on the "
       f"reduction since its refractivity in this set is zero")

# ---------------------------------------------------------------------------
# 5. Mean molecular refractivity
# ---------------------------------------------------------------------------
measured = float(x_h2[zero_nh3] * script_R[0] + x_he[zero_nh3] * script_R[1]
                 + x_nh3[zero_nh3] * script_R[2])
expected = (0.94 * 136.0 + 0.06 * 35.0) * 1e-6 / LOSCHMIDT_CONSTANT
record(5, "mean molecular refractivity at a level with zero ammonia, to six figures",
       abs(measured - expected) <= 1e-6 * abs(expected),
       f"measured {measured:.6e} m3\n"
       f"expected {expected:.6e} m3, from (0.94 x 136 + 0.06 x 35) x 1e-6 / n_Loschmidt\n"
       f"relative difference {abs(measured - expected) / expected:.1e}\n"
       f"per species script_R_i: "
       + ", ".join(f"{n} {v:.6e}" for n, v in zip(names, script_R)))

# ---------------------------------------------------------------------------
# 6. The file reads as kind C and refuses a perturbed copy
# ---------------------------------------------------------------------------
tampered = HERE / "tampered_composition.nc"
shutil.copy(PRODUCT, tampered)
with netCDF4.Dataset(tampered, "a") as handle:
    handle.variables["x_H2"][0] = float(x_h2[0]) + 1e-6
refused, message = False, "not refused"
try:
    cio.read(tampered, "composition")
except CasspianSchemaError as exc:
    refused, message = True, str(exc)
record(6, "read as kind C succeeds, and a copy with one mole fraction perturbed is refused",
       refused,
       f"composition_role = {root.attrs['composition_role']!r}; "
       f"vertical_coordinate = {root.attrs['vertical_coordinate']!r}\n"
       f"a field on {dict(field.sizes)} (SPEC_04 decision O), the source column's own latitude "
       f"latitude_planetographic_deg = {root.attrs['latitude_planetographic_deg']}\n"
       f"perturbing x_H2 at one level by 1e-6:\n  {message}")

# ---------------------------------------------------------------------------
# 7. The species group
# ---------------------------------------------------------------------------
attrs = tree["species"].attrs
expected_fields = ["molar_mass_kg_mol", "refractivity_per_molecule_m3",
                   "refractivity_uncertainty_m3", "refractivity_frequency_Hz",
                   "refractivity_temperature_dependence", "is_polar", "citation"]
present = [f for f in expected_fields if f in species_group.variables]
record(7, "the species group carries the master table values and its hash",
       len(present) == len(expected_fields) and "master_table_hash" in attrs,
       f"species: {names}\n"
       f"variables present: {present}\n"
       f"master_table = {attrs['master_table']!r}\n"
       f"master_table_hash = {attrs['master_table_hash']}\n"
       f"NH3 refractivity in this set is {script_R[2]:.1e} m3, which is the Fig. 3 caption "
       f"reading: Lindal computed density from H2 and He only\n"
       f"is_polar = {species_group['is_polar'].values.tolist()}, "
       f"temperature_dependence = "
       f"{species_group['refractivity_temperature_dependence'].values.tolist()}")

# ---------------------------------------------------------------------------
# 8. SPEC_01 v0.19 amendment: the closure declaration
# ---------------------------------------------------------------------------
share = x_h2 / (x_h2 + x_he)
record(8, "v0.19: closure_rule and closure_species read back and agree with the values",
       root.attrs.get("closure_rule") == "share_of_remainder"
       and root.attrs.get("closure_species") == "H2 He"
       and float(np.ptp(share)) <= 1e-12,
       f"closure_rule = {root.attrs.get('closure_rule')!r}\n"
       f"closure_species = {root.attrs.get('closure_species')!r}\n"
       f"x_H2 / (x_H2 + x_He) = {float(share.min())!r} to {float(share.max())!r}, "
       f"spread {float(np.ptp(share)):.1e}")

# ---------------------------------------------------------------------------
# 9. SPEC_01 v0.19 amendment: unstated per molecule uncertainties are NaN
# ---------------------------------------------------------------------------
unc = species_group["refractivity_uncertainty_m3"].values
record(9, "v0.19: refractivity_uncertainty_m3 is NaN for all three species",
       bool(np.all(np.isnan(unc))) and unc.size == 3,
       f"refractivity_uncertainty_m3 = {unc.tolist()} for {names}\n"
       f"uncertainty_method = "
       f"{species_group['refractivity_uncertainty_m3'].attrs.get('uncertainty_method')!r}")

tree.close()
print()
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
