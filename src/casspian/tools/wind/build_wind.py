"""Build a kind W wind file from a published zonal wind curve.

SPEC_01 v0.15 Step 7, SPEC_00 section 6.6. The Lindal instance: the cloud top wind Lindal
carried in his geoid, assembled from the Ingersoll and Pollard (1982) Fig. 5 solid curve,
converted to planetocentric latitude, brought to zero at both poles, and extended onto a
declared pressure grid under the source's altitude independence.

Manuscript equations implemented: none. The interpolation and the declared extensions are in
`curve.py`; this module reads, converts, assembles the file, and writes it.

**The wind and the uncertainty come from two different digitizations.** The Ingersoll and
Pollard curve is the wind: it is the authors' own smoothing of the Smith et al. (1982) cloud
tracking data, so using it removes every fitting choice from this step. The Smith points are
the uncertainty: their RMS about that curve, per 2 degree bin, is the 1 sigma. The two agree
on scale, 21.3 m/s RMS about the curve against 21.4 m/s pooled within bin scatter, which is
what makes the curve the mean curve of the data to the precision of the digitizations.

The point based path of `fit.py` is retained for a data set with no published curve, selected
by `wind_source_kind = "points"`. It is not exercised by this build.

**A tabulated wind** (SPEC_13 Step 1 item 2, section 3 ruling 1): `curve_format = "table"` reads
`curve_source` as a table with its own per-row standard deviation (`curve.read_table`), assembled
under `gap_rule = "pchip_bridge"` with `uncertainty_source = "table"` and the declared
`ring_gap_deg`. Such a source has no points file, no bins and no join, so `points_source`,
`bin_width_deg`, `min_count_per_bin` and `join_window_deg` belong to the default
`curve_format = "segments"` alone. A tabulated file's `source` is its own first line, with the
file's name.

**The data properties file carries only values that matter** (SPEC_13 v0.5 section 2a item 3): the
observation level, the latitude convention, the epoch and the season. Their sources, justification
and descriptions are its comments, and this module reads none of them.

No physical value appears in this module. The bin width, the minimum count, the gap and polar
rules, the join window and the two grids come from the control file; the observation level and
the latitude convention come from the static properties file beside the digitizations; the wind
comes from the curve.
"""

from __future__ import annotations

import argparse
import csv
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import geoid as gd
from casspian.lib import io as cio
from casspian.lib import latitude as latmod
from casspian.lib.control import ControlFileError, build_role, load_section
from casspian.tools.wind.curve import EXTENDED, AssembledCurve, TableCurve, read_curve, read_table

TOOL = "casspian-wind-from-curve"

#: The keys the data properties file's [epoch] table must carry (SPEC_01 v0.25 Step 7): the
#: epoch as an ISO date and the season (SPEC_13 v0.5 section 2a item 3: its source is a comment).
EPOCH_KEYS = ("value", "solar_longitude_deg")

SECTION_KEYS = {
    "role": True,
    "wind_source_kind": True,
    "curve_source": True,
    "points_source": False,
    "data_properties": True,
    "gravity_file": True,
    "rotation_file": True,
    "raw_bundle": True,
    "output": True,
    "prefix": True,
    "bin_width_deg": False,
    "min_count_per_bin": False,
    "gap_rule": True,
    "polar_rule": True,
    "join_window_deg": False,
    "curve_format": False,
    "uncertainty_source": False,
    "ring_gap_deg": False,
    "latitude_grid_step_deg": True,
    "vertical_structure": True,
    "pressure_grid_Pa": True,
    "title": False,
}

PATH_KEYS = (
    "curve_source", "points_source", "data_properties",
    "gravity_file", "rotation_file", "raw_bundle", "output",
)

#: The two forms of curve source, with the uncertainty source each takes and the keys each needs.
CURVE_FORMATS = {
    "segments": ("points", ("points_source", "bin_width_deg", "min_count_per_bin", "join_window_deg")),
    "table": ("table", ("ring_gap_deg",)),
}
PROVENANCE_FLAG_VALUES = np.array([0, 1, 2, 3, 4], dtype="int8")
PROVENANCE_FLAG_MEANINGS = (
    "observed interpolated parameterized extrapolated extended_by_source_assumption"
)


def read_points(path: Path):
    """Read the Smith digitized points: planetographic latitude in degrees, wind in m/s."""
    with open(Path(path), newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    latitude = np.array([float(r["latitude_planetographic_deg"]) for r in rows])
    wind = np.array([float(r["u_ms"]) for r in rows])
    return latitude, wind


def bin_points(latitude, wind, width_deg, min_count):
    """Bin in planetographic latitude. Returns centers, index, count, mean and the low flag.

    Edges fall on multiples of the width from -90, so with a width of 2 the centers are the odd
    degrees.
    """
    edges = np.arange(-90.0, 90.0 + width_deg, width_deg)
    centers = 0.5 * (edges[:-1] + edges[1:])
    index = np.clip(np.digitize(latitude, edges) - 1, 0, centers.size - 1)
    count = np.zeros(centers.shape, dtype="int64")
    mean = np.full(centers.shape, np.nan)
    for i in range(centers.size):
        members = wind[index == i]
        count[i] = members.size
        if members.size:
            mean[i] = members.mean()
    return centers, index, count, mean, (count > 0) & (count < int(min_count))


def bin_rms_about(curve, latitude, wind, index, centers, counts, min_count):
    """RMS of the points about the curve, per bin. NaN below the minimum count.

    This is the uncertainty of the wind: how far the cloud tracking points scatter from the
    curve drawn through them. A bin below the minimum count gets NaN, which SPEC_00 section 5
    makes the honest statement of no information.
    """
    residual = wind - curve(latitude)
    rms = np.full(centers.shape, np.nan)
    for i in np.flatnonzero(counts >= int(min_count)):
        members = residual[index == i]
        rms[i] = float(np.sqrt(np.mean(members**2)))
    return rms, residual


def build(control_path, section: str = "wind") -> Path:
    """Build the kind W file from the control file's named section."""
    control = load_section(control_path, section, SECTION_KEYS, path_keys=PATH_KEYS)
    kind = control["wind_source_kind"]
    if kind != "curve":
        raise ControlFileError(
            f"wind_source_kind is {kind!r}; this build is the curve path. The point based "
            "path of fit.py exists for a data set with no published curve and is not "
            "exercised here (SPEC_01 v0.15 Step 7)."
        )

    curve_format = control.get("curve_format", "segments")
    if curve_format not in CURVE_FORMATS:
        raise ControlFileError(f"curve_format = {curve_format!r}; the formats are {list(CURVE_FORMATS)}")
    uncertainty_source, needed = CURVE_FORMATS[curve_format]
    if control.get("uncertainty_source", uncertainty_source) != uncertainty_source:
        raise ControlFileError(
            f"uncertainty_source = {control['uncertainty_source']!r} with curve_format = "
            f"{curve_format!r}; that format takes {uncertainty_source!r}")
    missing = [key for key in needed if key not in control]
    if missing:
        raise ControlFileError(f"curve_format = {curve_format!r} needs {missing}")
    unused = [key for _, (_, keys) in CURVE_FORMATS.items() for key in keys
              if key not in needed and key in control]
    if unused:
        raise ControlFileError(f"curve_format = {curve_format!r} does not use {unused}")
    table_form = curve_format == "table"

    output = Path(control["output"])
    prefix = control["prefix"]
    role = build_role(control, control_path, section)
    pressure = np.asarray(control["pressure_grid_Pa"], dtype="float64")
    step = float(control["latitude_grid_step_deg"])

    with open(Path(control["data_properties"]), "rb") as handle:
        properties = tomllib.load(handle)
    epoch_table = properties.get("epoch", {})
    missing = [key for key in EPOCH_KEYS if key not in epoch_table]
    if missing:
        raise ControlFileError(
            f"{control['data_properties']}: the [epoch] table lacks {missing}; a wind file "
            "carries its epoch and season (SPEC_00 v0.17 section 5, SPEC_01 v0.25 Step 7)."
        )
    level_ref = float(properties["observation_level"]["value_Pa"])
    if not np.any(np.isclose(pressure, level_ref)):
        raise ControlFileError(
            f"the observation level {level_ref} Pa is not a node of the declared pressure grid"
        )

    gravity = cio.read(Path(control["gravity_file"]), "gravity")
    rotation = cio.read(Path(control["rotation_file"]), "rotation")
    raw = cio.read(Path(control["raw_bundle"]), "raw")
    try:
        degrees = gravity["degree"].values
        J = gravity["J"].values
        GM = float(gravity["GM_m3s2"])
        R_norm = float(gravity["normalization_radius_m"])
        rotation_name = str(rotation.attrs["system_name"])
        Omega = float(rotation["angular_rate_rad_s"])
        r_polar = float(
            raw["scalars/geodesy"].attrs["reference_geoid_polar_radius_km"]
        ) * 1e3
    finally:
        gravity.close()
        rotation.close()
        raw.close()

    if table_form:
        table = read_table(Path(control["curve_source"]))
        curve = TableCurve(table, gap_rule=control["gap_rule"], polar_rule=control["polar_rule"],
                           ring_gap_deg=control["ring_gap_deg"])
    else:
        width = float(control["bin_width_deg"])
        min_count = int(control["min_count_per_bin"])
        segments = read_curve(Path(control["curve_source"]))
        points_lat, points_u = read_points(Path(control["points_source"]))
        centers, index, count, mean, low = bin_points(points_lat, points_u, width, min_count)

        curve = AssembledCurve(
            segments,
            gap_rule=control["gap_rule"],
            polar_rule=control["polar_rule"],
            join_window_deg=control["join_window_deg"],
            bin_centers=centers,
            bin_values=mean,
        )
        rms, _ = bin_rms_about(curve, points_lat, points_u, index, centers, count, min_count)

    # The product grid is planetocentric, with both poles and the equator as nodes.
    n_steps = int(round(180.0 / step))
    phi_c_deg = np.linspace(-90.0, 90.0, n_steps + 1)
    if not np.any(np.abs(phi_c_deg) < 1e-12):
        raise ControlFileError(
            f"a latitude step of {step} deg does not put the equator on a node"
        )

    # Planetocentric to planetographic is a DIRECT evaluation; no fixed point is needed in
    # this direction (SPEC_01 v0.15 Step 7 item 4).
    def surface(phi):
        phi = np.atleast_1d(np.asarray(phi, dtype="float64"))
        r, _, _ = gd.reference_geoid(phi, r_polar, Omega, GM, J, degrees, R_norm,
                                     tol_m=1e-9, max_iter=80)
        return r

    zero_wind = lambda phi: np.zeros_like(np.asarray(phi, dtype="float64"))
    phi_g = np.degrees(
        latmod.planetographic_from_planetocentric(
            np.radians(phi_c_deg), surface, zero_wind, Omega, GM, J, degrees, R_norm
        )
    )
    phi_g[np.abs(phi_c_deg) >= 90.0] = np.sign(phi_c_deg[np.abs(phi_c_deg) >= 90.0]) * 90.0

    u_1d = curve(phi_g)
    provenance_1d = curve.provenance(phi_g)
    sigma_1d = (curve.uncertainty(phi_g) if table_form
                else curve.uncertainty(phi_g, centers, rms, count, min_count))

    n_lat, n_p = phi_c_deg.size, pressure.size
    u_total = np.repeat(u_1d[:, None], n_p, axis=1)
    u_uncertainty = np.repeat(sigma_1d[:, None], n_p, axis=1)
    provenance = np.repeat(provenance_1d[:, None], n_p, axis=1)
    reference_column = int(np.argmin(np.abs(pressure - level_ref)))
    away = np.ones(n_p, dtype=bool)
    away[reference_column] = False
    provenance[:, away] = EXTENDED

    # The three parts of SPEC_04 Step 0 deliverable 1: the reference level wind against
    # latitude, the total, and the shear along the local vertical, which is their difference.
    # This field is altitude independent, so the shear is zero and the reference level wind is
    # the curve itself; that is a property of these values, not a declaration.
    u_reference = u_total[:, reference_column].copy()
    u_shear = u_total - u_reference[:, None]

    def attrs(units, long_name, provenance_value, **extra):
        out = {"units": units, "long_name": long_name, "provenance": provenance_value}
        for key, value in extra.items():
            if value is None or (isinstance(value, str) and not value):
                continue
            out[key] = value
        return out

    dims = ("latitude_planetocentric", "pressure")
    dataset = xr.Dataset(
        {
            "u_total_ms": (dims, u_total,
                           attrs("m s-1", "total zonal wind, positive eastward", "derived",
                                 value_source=(None if table_form else
                                               "Ingersoll and Pollard 1982, Fig. 5, solid "
                                               "curve"))),
            "u_total_uncertainty_ms": (
                dims, u_uncertainty,
                attrs("m s-1", ("standard deviation of the measurements in the latitude bin"
                                if table_form else
                                "scatter of the cloud tracking points about the curve"),
                      "derived", uncertainty_kind="1sigma",
                      uncertainty_method=(None if table_form else (
                          "RMS of the Smith et al. 1982 Fig. 4 points about the Ingersoll and "
                          "Pollard 1982 Fig. 5 curve within a 2 degree planetographic bin; "
                          "NaN where the bin has fewer than the minimum count"))),
            ),
            "u_reference_ms": (
                ("latitude_planetocentric",), u_reference,
                attrs("m s-1", "zonal wind at the reference level", "derived"),
            ),
            "u_shear_ms": (dims, u_shear,
                           attrs("m s-1", "total minus the reference level wind, along the "
                                          "local vertical", "derived")),
            "value_provenance": (dims, provenance,
                                 attrs("1", "how each value was obtained", "index",
                                       flag_values=PROVENANCE_FLAG_VALUES,
                                       flag_meanings=PROVENANCE_FLAG_MEANINGS)),
            "reference_level_pressure_Pa": (
                (), level_ref,
                attrs("Pa", "level the cloud top wind is assigned to", "assumed"),
            ),
            "latitude_planetographic_deg": (
                ("latitude_planetocentric",), phi_g,
                attrs("degrees_north", "grid latitude in the source convention", "derived"),
            ),
        },
        coords={
            "latitude_planetocentric_deg": (
                ("latitude_planetocentric",), phi_c_deg,
                attrs("degrees_north", "product grid, planetocentric", "index"),
            ),
            "pressure_Pa": (
                ("pressure",), pressure,
                attrs("Pa", "declared pressure grid", "index",
                      positive="down", direction="increasing"),
            ),
        },
    )

    if table_form:
        dataset.attrs.update({
            "title": control.get("title", f"{prefix} cloud top zonal wind"),
            "profile_or_run": prefix,
            "role": role,
            "source": f"{table.header} ({Path(control['curve_source']).name})",
            "rotation_system_name": rotation_name,
            "rotation_rate_rad_s": Omega,
            "epoch": str(epoch_table["value"]),
            "solar_longitude_deg": float(epoch_table["solar_longitude_deg"]),
            "observation_level_Pa": level_ref,
            "source_latitude_convention": properties["latitude"]["convention"],
            "vertical_structure": control["vertical_structure"],
            "coverage_pressure_Pa": np.array([pressure.min(), pressure.max()]),
            "coverage_latitude_planetocentric_deg": np.array([-90.0, 90.0]),
            "curve_format": curve_format,
            "uncertainty_source": uncertainty_source,
            "gap_rule": control["gap_rule"],
            "polar_rule": control["polar_rule"],
            "ring_gap_deg": np.array([curve.ring_lo, curve.ring_hi]),
            "parameterization": (
                f"one PCHIP through all {table.latitude_deg.size} rows, which bridges the ring gap "
                f"{curve.ring_lo:g} to {curve.ring_hi:g} deg and the gaps "
                + ", ".join(f"{lo:g} to {hi:g}" for lo, hi in curve.gaps
                            if (lo, hi) != (curve.ring_lo, curve.ring_hi))
                + f" deg between their data edges, with no reflection and no join; the caps "
                f"poleward of {curve.north_max:g} and {curve.south_min:g} deg by the polar rule; "
                f"the uncertainty u_rms linear in latitude, NaN in the ring gap and poleward of "
                f"the rows; {table.dropped} rows the file marks as not data dropped"),
            "curve_segments_deg": np.array([curve.south_min, curve.north_max]),
            "latitude_conversion_inputs": cio.input_hashes(
                [Path(control["gravity_file"]), Path(control["rotation_file"])], output
            ),
            "latitude_conversion_note": (
                "grid latitudes converted from planetocentric to planetographic by the direct "
                "relation on the no wind reference geoid. The wind itself makes the geoid depart "
                "from no wind; using the no wind surface is a second order choice, recorded here "
                "rather than hidden"
            ),
            "input_hashes": cio.input_hashes([
                Path(control["curve_source"]), Path(control["data_properties"]),
                Path(control["gravity_file"]), Path(control["rotation_file"]),
                Path(control["raw_bundle"]), Path(control_path),
            ], output),
        })
        cio.history_append(
            dataset,
            f"{TOOL}: wind from {Path(control['curve_source']).name} ({table.latitude_deg.size} rows, "
            f"{table.dropped} dropped as not data) assembled under gap_rule={control['gap_rule']} "
            f"and polar_rule={control['polar_rule']}; uncertainty the table's u_rms; evaluated on "
            f"a {step} deg planetocentric grid and extended onto {n_p} pressure levels under the "
            "source's altitude independence",
        )
        return cio.write(output, dataset, "wind", created_by=TOOL)

    dataset = dataset.assign({
            "bin_rms_ms": (
                ("bin_latitude",), rms,
                attrs("m s-1", "RMS of the points about the curve in the bin", "derived"),
            ),
            "bin_count": (
                ("bin_latitude",), count,
                attrs("1", "cloud tracking points in the bin", "index"),
            ),
    })
    dataset = dataset.assign_coords({
            "bin_latitude_deg": (
                ("bin_latitude",), centers,
                attrs("degrees_north", "uncertainty bin center, planetographic", "index"),
            ),
    })
    north_min, north_max = curve.north_min, curve.north_max
    south_min, south_max = curve.south_min, curve.south_max
    dataset.attrs.update({
        "title": control.get("title", f"{prefix} cloud top zonal wind"),
        "profile_or_run": prefix,
        "role": role,
        "source": ("Ingersoll, A. P., and Pollard, D. 1982, Icarus 52, 62, Fig. 5, solid "
                   "curve; uncertainty from Smith, B. A., et al. 1982, Science 215, 504, "
                   "Fig. 4"),
        "wind_source": "ingersoll_pollard1982_fig5",
        "rotation_system_name": rotation_name,
        "rotation_rate_rad_s": Omega,
        "epoch": str(epoch_table["value"]),
        "solar_longitude_deg": float(epoch_table["solar_longitude_deg"]),
        "method": "cloud tracking, smoothed by the source",
        "observation_level_Pa": level_ref,
        "source_latitude_convention": properties["latitude"]["convention"],
        "vertical_structure": control["vertical_structure"],
        "coverage_pressure_Pa": np.array([pressure.min(), pressure.max()]),
        "coverage_latitude_planetocentric_deg": np.array([-90.0, 90.0]),
        "gap_rule": control["gap_rule"],
        "polar_rule": control["polar_rule"],
        "join_window_deg": np.asarray(control["join_window_deg"], dtype="float64"),
        "parameterization": (
            f"ring gap {curve.join_hi:g} to {north_min:g} deg filled by reflecting the "
            f"northern segment about the equator (the source's own dashed curve); the strip "
            f"within {north_min:g} deg of the equator bridged by an even cubic Hermite "
            f"matched in value and slope; blended linearly into the southern segment across "
            f"{control['join_window_deg']} deg"
        ),
        "curve_segments_deg": np.array([south_min, south_max, north_min, north_max]),
        "latitude_conversion_inputs": cio.input_hashes(
            [Path(control["gravity_file"]), Path(control["rotation_file"])], output
        ),
        "latitude_conversion_note": (
            "grid latitudes converted from planetocentric to planetographic by the direct "
            "relation on the no wind reference geoid. The wind itself makes the geoid depart "
            "from no wind; using the no wind surface is a second order choice, recorded here "
            "rather than hidden"
        ),
        "binning": (f"width {width} deg in planetographic latitude, edges on multiples of the "
                    f"width from -90; bins below {min_count} points carry NaN uncertainty"),
        "source_caption_curve": (
            "Zonal velocity profile (left) and its curvature or second derivative (right) for "
            "Saturn from Voyager 2 in late August, 1981. The data are from Smith et al. (1982) "
            "and are based on far fewer points than the similar curves for Jupiter. Data are "
            "missing at southern equatorial latitudes because that region was obscured by "
            "rings and ring shadow at the time of the Voyager encounters."
        ),
        "source_caption_points": (
            "Zonal (eastward) winds in the reference frame of Saturn's magnetic field (5). "
            "Each point refers to a single cloud feature that has been followed for one "
            "Saturnian rotation. Points from 7 deg to 20 deg S are from Voyager 1 (8)."
        ),
        "input_hashes": cio.input_hashes([
            Path(control["curve_source"]), Path(control["points_source"]),
            Path(control["data_properties"]), Path(control["gravity_file"]),
            Path(control["rotation_file"]), Path(control["raw_bundle"]),
            Path(control_path),
        ], output),
    })
    cio.history_append(
        dataset,
        f"{TOOL}: wind from the Ingersoll and Pollard 1982 Fig. 5 solid curve "
        f"({segments.north[0].size + segments.south[0].size} samples) assembled under "
        f"gap_rule={control['gap_rule']} and polar_rule={control['polar_rule']}; uncertainty "
        f"from {points_lat.size} Smith et al. 1982 Fig. 4 points binned at {width} deg; "
        f"evaluated on a {step} deg planetocentric grid and extended onto {n_p} pressure "
        f"levels under the source's altitude independence",
    )
    return cio.write(output, dataset, "wind", created_by=TOOL)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Write a kind W wind file from a published wind curve."
    )
    parser.add_argument("control", help="path to the TOML control file")
    parser.add_argument("--section", default="wind", help="section to read (default wind)")
    args = parser.parse_args(argv)
    print(f"wrote {build(args.control, args.section)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
