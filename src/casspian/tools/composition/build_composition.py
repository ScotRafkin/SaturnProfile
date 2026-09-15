"""Build the Lindal reduction composition, kind C.

SPEC_01 Step 8, SPEC_00 section 6.2. The composition Lindal assumed, so that the reduction
reproduces his refractivity to number density conversion: 94 percent H2 with the remainder He
and NH3, with the ammonia column of Table I filled under a declared rule.

Manuscript equations implemented: none directly. The mean molar mass and the mean molecular
refractivity this file feeds are formed in `lib.thermo` (Eq. A9) at the reduction, not here;
this tool writes the mole fractions and the per species properties they are formed from.

No physical value appears in this module. The species properties come from
`data_static/species_master.toml`, the ammonia column and the source's own stated composition
come from the raw bundle, and the fill rule and the split ratio are declared in the control
file.

**The split ratio is declared and then checked against the source.** SPEC_01 Step 8 has the
control file declare 0.94 and 0.06 as the ratio to hold. That is a choice about how to divide
the remainder, and it is not the same statement as the source's "94 percent H2"; the two happen
to coincide here. The tool refuses if the declared H2 share disagrees with the `h2_fraction`
the raw bundle transcribed from p. 1137, so the coincidence cannot drift apart unnoticed.
"""

from __future__ import annotations

import argparse
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import interp1d

from casspian.lib import io as cio
from casspian.lib.constants import CODATA_RELEASE, LOSCHMIDT_CONSTANT
from casspian.lib.control import ControlFileError, build_role, load_section

TOOL = "casspian-composition-lindal"

SECTION_KEYS = {
    # SPEC_01 v0.26: required, no default; also written as composition_role.
    "role": True,
    "raw_bundle": True,
    "species_master": True,
    "species_set": True,
    "output": True,
    "prefix": True,
    "h2_fraction_of_remainder": True,
    "he_fraction_of_remainder": True,
    "nh3_fill_rule": True,
    "vertical_coordinate": True,
    "title": False,
}

PATH_KEYS = ("raw_bundle", "species_master", "output")

#: SPEC_00 section 6.2, `refractivity_temperature_dependence`.
TEMPERATURE_FLAG_VALUES = np.array([0, 1], dtype="int8")
TEMPERATURE_FLAG_MEANINGS = "none declared_in_refractivity_temperature_note"
BOOLEAN_FLAG_VALUES = np.array([0, 1], dtype="int8")
BOOLEAN_FLAG_MEANINGS = "false true"


def fill_ammonia(pressure, nh3):
    """Fill the ammonia column under the declared rule, returning values and provenance.

    SPEC_01 v0.17 Step 8, in three parts, because the three regions are known in three
    different ways:

    * **Between the tabulated levels**, linear in pressure: `measured` where Table I gives a
      value, `interpolated` at the one interior gap.
    * **Below the deepest tabulated level**, linear extrapolation from the last two:
      `extrapolated`.
    * **Above the highest tabulated level**, exactly zero: `assumed`. Table I is blank there,
      and blank means not measured, not zero. What justifies zero is physics rather than the
      table, and the reason is written into the file as `nh3_aloft_assumption`.

    The v0.5 wording of the rule extrapolated upward and clamped at zero, which put 0.2 ppm at
    794.33 mbar and made the clamp level a puzzle (REPORT_01_step8, finding 1). There is no
    upward extrapolation and no clamp here: a value nobody measured is an assumption, and it is
    labelled as one.
    """
    pressure = np.asarray(pressure, dtype="float64")
    nh3 = np.asarray(nh3, dtype="float64")
    tabulated = np.isfinite(nh3)
    if tabulated.sum() < 2:
        raise ControlFileError("the ammonia column needs at least two tabulated values")

    lo, hi = float(pressure[tabulated].min()), float(pressure[tabulated].max())
    rule = interp1d(pressure[tabulated], nh3[tabulated], kind="linear",
                    bounds_error=False, fill_value="extrapolate")

    filled = np.asarray(rule(pressure), dtype="float64")
    aloft = pressure < lo
    filled[aloft] = 0.0

    provenance = np.where(
        aloft, "assumed",
        np.where(tabulated, "measured",
                 np.where(pressure <= hi, "interpolated", "extrapolated")),
    )
    return filled, provenance, lo, hi


def build(control_path, section: str = "composition") -> Path:
    """Build the kind C file from the control file's named section."""
    control = load_section(control_path, section, SECTION_KEYS, path_keys=PATH_KEYS)
    output = Path(control["output"])
    prefix = control["prefix"]
    role = build_role(control, control_path, section)
    set_name = control["species_set"]
    h2_share = float(control["h2_fraction_of_remainder"])
    he_share = float(control["he_fraction_of_remainder"])
    if abs(h2_share + he_share - 1.0) > 1e-12:
        raise ControlFileError(
            f"the declared remainder shares {h2_share} and {he_share} do not sum to one"
        )

    raw = cio.read(Path(control["raw_bundle"]), "raw")
    try:
        table1 = raw["table1"].dataset
        pressure = np.asarray(table1["pressure_Pa"].values, dtype="float64")
        nh3_raw = np.asarray(table1["nh3_mole_fraction"].values, dtype="float64")
        composition = raw["scalars/composition"].attrs
        latitude_group = raw["scalars/latitude"].attrs
        source_group = raw["scalars/source"].attrs
        source_h2 = float(composition["h2_fraction"])
        h2_uncertainty = float(composition["h2_uncertainty"])
        footnote = str(source_group["table1_footnote"])
        observation_date = str(source_group["observation_date"])
        latitude_planetographic = float(latitude_group["planetographic_deg"])
        latitude_value_source = str(latitude_group["value_source"])
        raw_hash = cio.input_hash_entry(Path(control["raw_bundle"]), output)
    finally:
        raw.close()

    if abs(h2_share - source_h2) > 1e-12:
        raise ControlFileError(
            f"the control file declares an H2 share of {h2_share} but the raw bundle "
            f"transcribes the source's h2_fraction as {source_h2}. The split ratio is a "
            "choice and the source statement is data; they may be made to differ "
            "deliberately, but not silently."
        )

    with open(Path(control["species_master"]), "rb") as handle:
        master = tomllib.load(handle)
    if set_name not in master.get("refractivity", {}):
        raise ControlFileError(
            f"{control['species_master']}: no refractivity set named {set_name!r}; it holds "
            f"{sorted(master.get('refractivity', {}))}"
        )
    refractivity_set = master["refractivity"][set_name]
    molar_mass = master["molar_mass"]

    nh3, nh3_provenance, p_lo, p_hi = fill_ammonia(pressure, nh3_raw)
    x_h2 = h2_share * (1.0 - nh3)
    x_he = he_share * (1.0 - nh3)

    # The species carried in this set: those the source's own statement names.
    species = ["H2", "He", "NH3"]
    fractions = {"H2": x_h2, "He": x_he, "NH3": nh3}
    provenance = {"H2": "assumed", "He": "assumed", "NH3": None}

    def attrs(units, long_name, provenance_value, **extra):
        out = {"units": units, "long_name": long_name, "provenance": provenance_value}
        for key, value in extra.items():
            if value is None or (isinstance(value, str) and not value):
                continue
            out[key] = value
        return out

    variables = {
        "pressure_Pa": (
            ("level",), pressure,
            attrs("Pa", "pressure of the source profile", "index",
                  positive="down", direction="increasing"),
        ),
    }
    for name in species:
        if name == "NH3":
            # One provenance per variable, so the dominant rule is named and the per level
            # detail goes to a companion flag.
            unique = sorted(set(nh3_provenance.tolist()))
            label = "measured" if unique == ["measured"] else "interpolated"
            variables["x_NH3"] = (
                ("level",), nh3,
                attrs("mol/mol", "ammonia mole fraction", label,
                      value_source="Table I, filled by the declared rule",
                      fill_rule=control["nh3_fill_rule"]),
            )
        else:
            variables[f"x_{name}"] = (
                ("level",), fractions[name],
                attrs("mol/mol", f"{name} mole fraction", provenance[name],
                      value_source=str(composition.get("h2_source", ""))[:400]),
            )
    variables["x_H2_uncertainty"] = (
        ("level",), np.full(pressure.shape, h2_uncertainty),
        attrs("mol/mol", "uncertainty on the H2 mole fraction", "assumed",
              uncertainty_kind="stated",
              uncertainty_method="the source's stated 94 +/- 3 percent (p. 1137)"),
    )
    variables["nh3_provenance"] = (
        ("level",),
        np.array([{"measured": 0, "interpolated": 1, "assumed": 2, "extrapolated": 3}[p]
                  for p in nh3_provenance], dtype="int8"),
        attrs("1", "how the ammonia value at this level was obtained", "index",
              flag_values=np.array([0, 1, 2, 3], dtype="int8"),
              flag_meanings="measured interpolated assumed extrapolated"),
    )

    dataset = xr.Dataset(variables)
    dataset.attrs.update({
        "title": control.get("title", f"{prefix} {role} composition"),
        "profile_or_run": prefix,
        "role": role,
        "source": str(source_group.get("citation", "")),
        "composition_role": role,
        # SPEC_01 v0.25 Step 8: the composition is an assumption of the source, valid at the
        # observation's date and declared season independent.
        "epoch": observation_date,
        "season_absent_meaning": "uniform",
        "vertical_coordinate": control["vertical_coordinate"],
        "source_statement": footnote,
        "latitude_planetocentric_absent_meaning": "point",
        "latitude_planetographic_deg": latitude_planetographic,
        "latitude_planetographic_deg_value_source": latitude_value_source,
        "species_set": set_name,
        "remainder_split": (
            f"ammonia assigned first from Table I under the declared rule, then the remainder "
            f"split {h2_share} H2 to {he_share} He"
        ),
        # SPEC_00 section 6.2 v0.13: the closure is declared, so refrac reads it rather than
        # inferring it. Ammonia is assigned first and H2 and He share the remainder, the share
        # species named first.
        "closure_rule": "share_of_remainder",
        "closure_species": "H2 He",
        "nh3_fill_rule": control["nh3_fill_rule"],
        "nh3_aloft_assumption": (
            "Above the highest tabulated level Table I is blank, which is 'not measured', not "
            "'zero'. Zero is assigned there as an assumption, justified by physics rather than "
            "by the table: Lindal states the troposphere at these levels is saturated with "
            "ammonia, and the saturation mixing ratio over NH3 ice on his own temperatures "
            "falls from a few ppm near 0.8 bar (consistent with the tabulated 2.6 ppm at "
            "831.76 mbar) to below 1e-7 by 0.5 bar and below 1e-12 at the tropopause. The "
            "assumption costs under 1e-6 of the mean molar mass and nothing in the "
            "refractivity, since the lindal1985 set carries NH3 at zero. Levels so filled "
            "carry nh3_provenance = assumed."
        ),
        "nh3_tabulated_range_Pa": np.array([p_lo, p_hi]),
        "codata_release": CODATA_RELEASE,
        "raw_bundle": raw_hash,
        "input_hashes": cio.input_hashes([
            Path(control["raw_bundle"]), Path(control["species_master"]),
            Path(control_path),
        ], output),
    })

    # ---- the species group, SPEC_00 section 6.2 --------------------------------------
    def property_of(name, key, default=0.0):
        return refractivity_set.get(name, {}).get(key, default)

    def molecular_property(name, key):
        """A property of the molecule, read from its own table beside `molar_mass`.

        `is_polar` is a property of the molecule and not of a refractivity set. It used to sit
        inside each set, where the `lindal1985` entry for NH3 omitted it, so reading it from
        that set alone would have written ammonia as nonpolar (Step 8 review, finding 2). It
        now has its own table and is read from there; a species the table is silent about is
        refused rather than defaulted.
        """
        table = master.get(key, {})
        if name not in table:
            raise ControlFileError(
                f"{control['species_master']}: no [{key}] entry for {name}. It is a property "
                "of the molecule and cannot be defaulted."
            )
        return table[name]

    per_molecule = np.array(
        [float(property_of(n, "refractivity_stp_e6")) * 1e-6 / LOSCHMIDT_CONSTANT
         for n in species], dtype="float64"
    )
    # An uncertainty the master table does not state is NaN, never 0.0 (SPEC_00 section 5,
    # v0.7; SPEC_01 v0.19): the table writes it as nan, or omits the key, and both arrive here
    # as NaN, so a stated zero and an unstated value stay distinguishable in the product.
    per_molecule_uncertainty = np.array(
        [float(property_of(n, "uncertainty_e6", float("nan"))) * 1e-6 / LOSCHMIDT_CONSTANT
         for n in species], dtype="float64"
    )
    species_group = xr.Dataset(
        {
            "molar_mass_kg_mol": (
                ("species",), np.array([float(molar_mass[n]["value"]) for n in species]),
                attrs("kg mol-1", "molar mass", "measured",
                      status="\n".join(f"{n}: {molar_mass[n]['status']}" for n in species)),
            ),
            "refractivity_per_molecule_m3": (
                ("species",), per_molecule,
                attrs("m3", "refractivity per molecule, script_R_i", "derived",
                      status="\n".join(f"{n}: {property_of(n, 'status', 'unstated')}"
                                       for n in species),
                      conversion=("refractivity_stp_e6 * 1e-6 / n_Loschmidt, with "
                                  f"n_Loschmidt = {LOSCHMIDT_CONSTANT:.9e} m-3 from "
                                  f"lib.constants, CODATA {CODATA_RELEASE}")),
            ),
            "refractivity_uncertainty_m3": (
                ("species",), per_molecule_uncertainty,
                attrs("m3", "uncertainty on the per molecule refractivity", "derived",
                      uncertainty_kind="stated",
                      uncertainty_method=("NaN where the master table states no uncertainty "
                                          "(SPEC_00 section 5); see status")),
            ),
            "refractivity_frequency_Hz": (
                ("species",),
                np.array([float(property_of(n, "frequency_Hz")) for n in species]),
                attrs("Hz", "frequency the refractivity applies to", "measured"),
            ),
            "refractivity_temperature_dependence": (
                ("species",),
                np.array([0 if property_of(n, "temperature_dependence", "none").startswith("none")
                          else 1 for n in species], dtype="int8"),
                attrs("1", "whether a temperature law is declared", "index",
                      flag_values=TEMPERATURE_FLAG_VALUES,
                      flag_meanings=TEMPERATURE_FLAG_MEANINGS),
            ),
            "is_polar": (
                ("species",),
                np.array([1 if molecular_property(n, "is_polar") else 0 for n in species],
                         dtype="int8"),
                attrs("1", "whether the molecule is polar", "index",
                      flag_values=BOOLEAN_FLAG_VALUES, flag_meanings=BOOLEAN_FLAG_MEANINGS),
            ),
        },
        coords={"species_name": (("species",), np.array(species, dtype=object))},
    )
    species_group["citation"] = (
        ("species",),
        np.array([str(property_of(n, "source", "")) for n in species], dtype=object),
    )
    species_group.attrs.update({
        # Relative to the file being written, so the record carries no machine's paths
        # (SPEC_00 v0.18 section 5).
        "master_table": cio.recorded_path(Path(control["species_master"]), output),
        "master_table_hash": cio.input_hash_entry(Path(control["species_master"]), output),
        "set": set_name,
        "refractivity_temperature_note": str(
            refractivity_set.get("NH3", {}).get("note", "")
        )[:2000],
        "nh3_note": (
            "The Lindal set carries NH3 refractivity zero. Lindal et al. (1985) Fig. 3 caption "
            "states the density was computed assuming 94 percent hydrogen with the remainder "
            "mostly helium, that is from H2 and He only, so the ammonia contributes no "
            "refractivity in their reduction. Its mole fraction still enters the mean molar "
            "mass."
        ),
    })

    cio.history_append(
        dataset,
        f"{TOOL}: ammonia column filled from Table I by the declared rule between "
        f"{p_lo:g} and {p_hi:g} Pa and clamped at zero; remainder split {h2_share} to "
        f"{he_share}; species properties from the {set_name} set of the master table, "
        f"refractivities converted to m3 per molecule by the Loschmidt constant",
    )
    return cio.write(output, dataset, "composition", groups={"species": species_group},
                     created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write the Lindal reduction composition, kind C."
    )
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="composition", help="section (default composition)")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.control, args.section)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
