"""The kind N product and the `casspian-refrac` entry point.

SPEC_02 Step 4; SPEC_00 section 6.7 (v0.13). `refrac` in one call: read the manifest, load the
six inputs, freeze the anchor, reduce, and write `<profile>_refractivity.nc` with everything it
consumed embedded, so the file can be audited alone.

This module converts the reduction's radians to the degrees a stored scalar carries (SPEC_00
section 5), once, at the write. Nothing here computes physics beyond the two extra marches of
the anchor rule spread, which are recorded as a choice and enter no number of the product.
"""

from __future__ import annotations

import argparse
import dataclasses
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib.schema import uncertainty_companion
from casspian.refrac.anchor import freeze_anchor
from casspian.refrac.reduce import reduce_profile

__all__ = ["TOOL", "build_product", "main"]

TOOL = "casspian-refrac"
METHOD = ("first order propagation of the declared input uncertainties, uncorrelated between "
          "inputs, each converted to one standard deviation by its declared kind, unstated "
          "terms left out (SPEC_02 Step 3; SPEC_00 sections 5 and 6.7)")


def _repository_root(start: Path) -> Path:
    start = Path(start).resolve()
    for candidate in [start, *start.parents]:
        if (candidate / ".git").exists():
            return candidate
    return start.parent


def _attrs(units, long_name, provenance, **extra):
    out = {"units": units, "long_name": long_name, "provenance": provenance}
    out.update({k: v for k, v in extra.items() if v is not None})
    return out


def _propagated(companion, units, long_name, scale=1.0):
    """A propagated companion as a variable payload: value and the SPEC_00 attributes."""
    return np.asarray(companion.value, dtype="float64") * scale, _attrs(
        units, long_name, "derived",
        uncertainty_kind="1sigma",
        uncertainty_method=METHOD,
        uncertainty_terms_included=" ".join(companion.included),
        uncertainty_terms_unstated=" ".join(companion.unstated),
        uncertainty_kind_conversions="\n".join(companion.conversions),
    )


def _copied(variable_or_value, attrs_from, units, long_name, default_kind="stated"):
    """A companion copied from an input, with its declared kind kept, not converted."""
    return np.float64(variable_or_value), _attrs(
        units, long_name, attrs_from.get("provenance", "derived"),
        uncertainty_kind=str(attrs_from.get("uncertainty_kind", default_kind)),
        uncertainty_method=str(attrs_from.get("uncertainty_method",
                                              "copied from the input file")),
    )


def _anchor_rule_spread(inputs, manifest):
    """The frozen pair under every other manifest anchor rule. A choice, recorded.

    SPEC_02 v0.7 Step 4 and v0.8 Step 6: the whole Step 2 fixed point is rerun under each rule
    other than the one in force, each anchored on its own kind D quantity, so the recorded
    spread is the whole effect of the choice. A rule that moves the surface also moves the
    frozen latitude, through `dphi_c/dr_anchor`.
    """
    spread = {}
    for rule, quantity in ctl.ANCHOR_QUANTITY_FOR_RULE.items():
        if rule == manifest.anchor_rule:
            continue
        frozen = freeze_anchor(inputs, dataclasses.replace(manifest, anchor_rule=rule,
                                                           anchor_quantity=quantity))
        spread[rule] = (frozen.r0_m, frozen.phi_c_rad)
    return spread


def build_product(manifest_path) -> Path:
    """Run the reduction for one manifest and write kind N. Returns the path written."""
    manifest = ctl.read_reduction_manifest(manifest_path)
    inputs = ctl.load_reduction_inputs(manifest)
    anchor = freeze_anchor(inputs, manifest)
    red = reduce_profile(inputs, manifest, anchor)
    spread = _anchor_rule_spread(inputs, manifest)

    thermo = inputs.thermo
    k = red.anchor_level
    C = red.companions
    deg = np.degrees

    level_vars = {
        "height_above_anchor_isobar_m": (
            ("level",), red.height_above_anchor_isobar_m,
            _attrs("m", "height above the anchor isobar, h - h_ref (Eq. B1)", "derived",
                   positive="up")),
        "number_density_m3": (
            ("level",), red.number_density_m3,
            _attrs("m-3", "number density, p / (k_B T) (Eq. B3.1)", "derived")),
        "refractivity": (
            ("level",), red.refractivity,
            _attrs("1", "refractivity N = n_refr - 1, unscaled (Eq. B3.3)", "derived")),
        "mean_refractivity_m3": (
            ("level",), red.mean_refractivity_m3,
            _attrs("m3", "mean refractivity per molecule, sum x_i script_R_i", "derived")),
        "mean_molar_mass_kg_mol": (
            ("level",), red.mean_molar_mass_kg_mol,
            _attrs("kg mol-1", "mean molar mass, sum x_i M_i", "derived")),
    }
    for name, units, long_name in (
        ("radius_uncertainty_m", "m", "uncertainty on the absolute radius"),
        ("height_above_anchor_isobar_uncertainty_m", "m",
         "uncertainty on the height above the anchor isobar"),
        ("number_density_uncertainty_m3", "m-3", "uncertainty on the number density"),
        ("refractivity_uncertainty", "1", "uncertainty on the refractivity"),
        ("mean_refractivity_uncertainty_m3", "m3", "uncertainty on the mean refractivity"),
        ("mean_molar_mass_uncertainty_kg_mol", "kg mol-1", "uncertainty on the mean molar mass"),
    ):
        value, attrs = _propagated(C[name], units, long_name)
        level_vars[name] = (("level",), value, attrs)

    latitude_name = "latitude_planetographic_deg"
    label_uncertainty = uncertainty_companion(latitude_name)
    scalar_vars = {
        "latitude_planetocentric_deg": (
            (), np.float64(deg(anchor.phi_c_rad)),
            _attrs("degrees_north", "planetocentric latitude of the anchor, phi_c", "derived")),
        "psi_deg": (
            (), np.float64(deg(anchor.psi_rad)),
            _attrs("degrees", "tilt of the local vertical from the radius at phi_c", "derived")),
        latitude_name: (
            (), np.float64(thermo.attrs[latitude_name]),
            _attrs("degrees_north", "planetographic latitude of the profile label, from kind T",
                   "measured",
                   value_source=str(thermo.attrs.get(f"{latitude_name}_value_source", "")),
                   swath=np.asarray(thermo.attrs.get(f"{latitude_name}_swath", []),
                                    dtype="float64"),
                   swath_source=str(thermo.attrs.get(f"{latitude_name}_swath_source", "")),
                   latitude_definition=str(thermo.attrs.get("latitude_definition", "")))),
        "anchor_isobar_pressure_Pa": (
            (), np.float64(red.pressure_Pa[k]),
            _attrs("Pa", "pressure of the anchor isobar", "index")),
        "anchor_isobar_radius_m": (
            (), np.float64(anchor.r0_m),
            _attrs("m", "absolute radius of the anchor isobar at phi_c, r0", "derived")),
        "anchor_isobar_height_m": (
            (), np.float64(red.h_ref_m),
            _attrs("m", "tabulated height of the anchor isobar, h_ref", "measured")),
    }
    for name, companion_key, units, long_name, scale in (
        ("latitude_planetocentric_uncertainty_deg", "latitude_planetocentric_uncertainty_rad",
         "degrees", "uncertainty on phi_c", float(deg(1.0))),
        ("psi_uncertainty_deg", "psi_uncertainty_rad", "degrees", "uncertainty on psi",
         float(deg(1.0))),
        ("anchor_isobar_radius_uncertainty_m", "anchor_isobar_radius_uncertainty_m", "m",
         "uncertainty on r0", 1.0),
    ):
        value, attrs = _propagated(C[companion_key], units, long_name, scale)
        scalar_vars[name] = ((), np.float64(value), attrs)
    scalar_vars[label_uncertainty] = (
        (), np.float64(thermo.attrs.get(f"{latitude_name}_uncertainty", np.nan)),
        _attrs("degrees", "uncertainty on the planetographic label, as declared in kind T",
               "measured",
               uncertainty_kind=str(thermo.attrs.get(f"{latitude_name}_uncertainty_kind",
                                                     "stated")),
               uncertainty_method="copied from kind T as declared, not converted"))
    value, attrs = _copied(thermo["pressure_uncertainty_Pa"].values[k],
                           thermo["pressure_uncertainty_Pa"].attrs, "Pa",
                           "uncertainty on the anchor isobar pressure, from kind T")
    scalar_vars["anchor_isobar_pressure_uncertainty_Pa"] = ((), value, attrs)
    value, attrs = _copied(thermo["height_uncertainty_m"].values[k],
                           thermo["height_uncertainty_m"].attrs, "m",
                           "uncertainty on the anchor isobar height, from kind T")
    scalar_vars["anchor_isobar_height_uncertainty_m"] = ((), value, attrs)

    dataset = xr.Dataset(
        {**level_vars, **scalar_vars},
        coords={"radius_m": (("level",), red.radius_m,
                             _attrs("m", "absolute radius from the center of mass (Eq. B1)",
                                    "derived", positive="up"))},
    )
    value, attrs = _propagated(C["radius_uncertainty_m"], "m",
                               "uncertainty on the absolute radius")
    dataset["radius_uncertainty_m"] = (("level",), value, attrs)

    root_dir = _repository_root(manifest.path)
    input_paths = [manifest.inputs[key] for key in ctl.INPUT_KINDS] + [manifest.path]
    dataset.attrs.update({
        "title": f"{manifest.slug} registered refractivity, kind N",
        "profile_or_run": manifest.slug,
        "role": "reduction",
        "source": (f"{TOOL} reduction of the {manifest.slug} profile under "
                   f"{manifest.path.name}; source profile: {thermo.attrs.get('source', '')}"),
        "input_hashes": cio.input_hashes(input_paths, relative_to=root_dir),
    })
    cio.history_append(
        dataset,
        f"{TOOL}: phi_c and r0 frozen on the wind geoid ({anchor.anchor_rule}), number density, "
        f"absolute radius and refractivity on {red.pressure_Pa.size} levels, companions by "
        "first order propagation",
    )

    record = {
        "fixed_point_iterates_deg": deg(anchor.iterates_rad),
        "fixed_point_iteration_count": np.int32(anchor.iteration_count),
        "fixed_point_final_step_deg": float(deg(anchor.final_step_rad)),
        "fixed_point_tolerance_deg": float(manifest.fixed_point_tolerance_deg),
        "fixed_point_closure_deg": float(deg(anchor.fixed_point_closure_rad)),
        "anchor_rule": anchor.anchor_rule,
        "anchor_surface_Pa": float(anchor.anchor_surface_Pa),
        "anchor_quantity": anchor.anchor_quantity,
        "anchor_radius_m": float(anchor.r_anchor_m),
        "north_polar_start_m": float(anchor.north_start_m),
        "polar_radius_north_m": float(anchor.polar_north_m),
        "polar_radius_south_m": float(anchor.polar_south_m),
        "polar_asymmetry_m": float(anchor.polar_asymmetry_m),
        "anchoring_residual_m": float(anchor.anchor_residual_m),
        "u_at_anchor_ms": float(anchor.u_at_phi_c_ms),
        "nowind_latitude_planetocentric_deg": float(deg(anchor.nowind_phi_c_rad)),
        "nowind_psi_deg": float(deg(anchor.nowind_psi_rad)),
        "nowind_anchor_isobar_radius_m": float(anchor.nowind_r0_m),
        "nowind_fixed_point_iterates_deg": deg(anchor.nowind_iterates_rad),
        "polar_radius_mean_m": float(0.5 * (anchor.polar_north_m + anchor.polar_south_m)),
        "equatorial_radius_marched_m": float(anchor.equator_radius_m),
        "polar_radii_note": (
            "Both polar radii, their mean and asymmetry, and the equatorial radius are those of "
            "the final march under the rule in force. Whichever quantity the rule does not "
            "anchor on is a prediction of the march, to compare with the source (SPEC_02 v0.8 "
            "Step 6)."),
        "anchor_rule_spread_note": (
            "r0 and phi_c under every other manifest anchor rule, each from the whole anchor "
            "fixed point rerun under that rule on its own kind D quantity, so the latitude "
            "shift the rule causes is included (SPEC_02 v0.7 Step 4, v0.8 Step 6). A choice of "
            "rule, not a declared uncertainty, so it enters no companion (SPEC_02 decision 5)."),
        "closure_rule": red.closure["rule"],
        "closure_species": f"{red.closure['share']} {red.closure['partner']}",
        "codata_release": red.codata_release,
        "casspian_version": cio.package_version(),
        "casspian_git_commit": cio.git_commit(),
    }
    for rule, (r0_rule, phi_rule) in spread.items():
        record[f"anchor_rule_spread_r0_{rule}_m"] = float(r0_rule)
        record[f"anchor_rule_spread_latitude_planetocentric_{rule}_deg"] = float(deg(phi_rule))
    for key, value in red.partials.items():
        record[f"partial_{key}"] = value if isinstance(value, str) else float(value)
    for key, value in red.terms.items():
        record[f"term_{key}"] = float(value)
    for key in ctl.INPUT_KINDS:
        record[f"input_sha256_{key}"] = inputs.sha256[key]
    record["input_sha256_manifest"] = manifest.sha256

    groups = {}
    for key in ctl.INPUT_KINDS:
        loaded = getattr(inputs, key)
        if isinstance(loaded, xr.DataTree):
            for node in loaded.subtree:
                suffix = node.path.strip("/")
                path = f"inputs/{key}" + (f"/{suffix}" if suffix else "")
                groups[path] = node.to_dataset(inherit=False)
        else:
            groups[f"inputs/{key}"] = loaded
    groups["manifest"] = xr.Dataset(attrs={"text": manifest.text, "sha256": manifest.sha256})
    groups["reduction_record"] = xr.Dataset(attrs=record)

    return cio.write(manifest.product, dataset, "refractivity", groups=groups, created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Reduce one profile to its kind N refractivity product.")
    parser.add_argument("manifest", help="path to <profile>_reduction.toml")
    args = parser.parse_args(argv)
    manifest = ctl.read_reduction_manifest(args.manifest)
    product = build_product(manifest.path)
    print(f"wrote {product}")
    diagnostics = manifest.diagnostics
    if diagnostics.get("figures"):
        # SPEC_02 Step 5: the standard figures, drawn from the file just written, never from
        # this run's memory. Imported here so refrac does not load matplotlib unless asked.
        from casspian.tools.plots import describe, render

        result = render(product, product.parent / "figures", diagnostics["format"],
                        int(diagnostics["dpi"]))
        for line in describe(result):
            print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
