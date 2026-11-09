"""`casspian-section`: an anchor transferred to every latitude, as a latitude-pressure section.

SPEC_20 v0.4 section 2 deliverable 2, REPORT_20_preexecution section 2, REVIEW_20_step0. No new physics.
The section is a byproduct of two transfers, one from the anchor toward each pole (Step 0 ruling):

1. A base run supplies the namelist and the inputs (composition, gravity, rotation and the wind, so the
   base run's wind is the section's). Its namelist is copied into the section's work directory with
   only the anchor and the target changed, its inputs read where they are.
2. Each march is the model's own `forward.transfer.chain`, on a mesh whose exact nodes include every
   section latitude on its side of the anchor, so the curves pass through each one.
3. At each section latitude the production is redone with `forward.transfer.produce_at_target`, the
   function that makes a run's delivered profile, so a section latitude delivers what a run ending
   there would. Reading the isobar label alone would miss that by up to 0.3 K (REPORT_20_preexecution
   section 2.1).

The lapse rate is `-dT/dz` along the local vertical, the centered difference on the column's height.

Writes one netCDF file, the product, which nothing reads; the control file is committed and the
product ignored (REVIEW_20_step0 ruling 5).

Manuscript equations implemented: none. No physical value appears in this module.
"""

from __future__ import annotations

import argparse
import math
import os
import re
import shutil
import tomllib
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import estimate as fe
from casspian.forward import propagate as fprop
from casspian.forward import transfer as ft
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf

TOOL = "casspian-section"


def read_control(path) -> dict:
    path = Path(path).resolve()
    with open(path, "rb") as handle:
        c = tomllib.load(handle)["section"]
    here = path.parent
    root = (here / c["repository_root"]).resolve()
    return {"path": path, "here": here, "root": root, "name": c["name"], "base_run": c["base_run"],
            "anchor_slug": c["anchor_slug"], "anchor_path": (root / c["anchor"]).resolve(),
            "step_deg": float(c["latitude_step_deg"]), "limit_deg": float(c["latitude_limit_deg"]),
            "composition": (root / c["composition"]).resolve() if "composition" in c else None,
            "output": here / c["output"]}


def grid(control) -> np.ndarray:
    """Planetocentric latitudes, every step from the south limit to the north."""
    count = int(round(2 * control["limit_deg"] / control["step_deg"]))
    return np.linspace(-control["limit_deg"], control["limit_deg"], count + 1)


def namelist(control, target_deg, work: Path) -> Path:
    """The base run's namelist with only the anchor and the target changed, its inputs where they are."""
    base, root = control["base_run"], control["root"]
    source = root / "forward" / base
    name = f"{control['name']}_{'north' if target_deg > 0 else 'south'}"
    directory = work / name
    if directory.exists():
        shutil.rmtree(directory)
    directory.mkdir(parents=True)
    text = (source / f"{base}.toml").read_text(encoding="utf-8")
    inputs = os.path.relpath(source / "inputs", directory).replace(os.sep, "/")
    anchor = os.path.relpath(control["anchor_path"], directory).replace(os.sep, "/")
    text = text.replace(f'"inputs/{base}_', f'"{inputs}/{base}_').replace(f'"{base}', f'"{name}')
    text = re.sub(r'(?m)^slug\s*=.*$', f'slug   = "{control["anchor_slug"]}"', text)
    text = re.sub(r'(?m)^path\s*=.*$', f'path   = "{anchor}"', text)
    text = re.sub(r'(?m)^latitude_planetocentric_deg\s*=.*$', f"latitude_planetocentric_deg = {target_deg!r}", text)
    if control["composition"] is not None:
        # The run's composition must be on the anchor's levels (forward.transfer.place produces the anchor on
        # them), so an anchor other than the base run's brings its own, the same mixture on its levels.
        composition = os.path.relpath(control["composition"], directory).replace(os.sep, "/")
        text = re.sub(r'(?m)^composition\s*=.*$', f'composition = "{composition}"', text)
    text = text.replace("figures = true", "figures = false")
    # The datum only places a run's altitude, which the section does not carry (it carries the height above
    # the anchor isobar); 100 mbar is a level of every anchor, where 1 bar is not of the Voyager 1 ingress.
    text = re.sub(r'(?m)^datum_isobar_Pa\s*=.*$', "datum_isobar_Pa = 1.0e4", text)
    path = directory / f"{name}.toml"
    path.write_bytes(text.encode("utf-8"))
    return path


def march(path: Path, latitudes_deg) -> dict:
    """One transfer from the anchor to the namelist's target, produced at each of `latitudes_deg`."""
    nl = ctl.read_run_namelist(path)
    inputs = ctl.load_run_inputs(nl)
    anchors = tuple(fprop.propagate(a, nl.solar_longitude_deg) for a in inputs.anchors)
    gauge = fe.gauge_latitude(anchors)
    target = math.radians(float(nl.target_latitude_deg))
    p_b = ft.boundary_pressure(anchors) if nl.p_b_Pa is None else float(nl.p_b_Pa)
    field = wf.WindField(inputs.wind)
    spacing = math.radians(float(nl.grid["latitude_spacing_deg"]))
    Phi_spacing = float(nl.grid["geopotential_spacing_m2s2"])
    at = np.radians(np.asarray(latitudes_deg, dtype="float64"))
    first = [ft.place(a, inputs, field, p_b) for a in anchors]
    required = np.concatenate([[math.radians(a.latitude_planetocentric_deg) for a in anchors], [gauge, target], at])
    mesh = lm.build_mesh(required, [p.geopotential_m2s2 for p in first], spacing, Phi_spacing)
    loop = nl.numerics["outer_loop"]
    carried = ft.chain(inputs, anchors, gauge_latitude_rad=gauge, target_latitude_rad=target,
                       gauge_isobar_Pa=nl.gauge_isobar_Pa, p_b=p_b, latitude_spacing_rad=spacing,
                       geopotential_spacing_m2s2=Phi_spacing, relative_tolerance_ln_p=loop["relative_tolerance_ln_p"],
                       max_iterations=loop["max_iterations"],
                       kernel_uncertainty_per_rad=float(nl.estimation["kernel_uncertainty_per_rad"]),
                       datum_isobar_Pa=float(nl.datum_isobar_Pa), field=field, mesh=mesh)
    state, estimate = carried.state, carried.estimate
    curves, _, along = ft.carry_to(state, inputs, gauge, target, estimate.geopotential_m2s2, estimate.C,
                                   estimate.label_pressure_Pa)
    out = {}
    for phi_deg, phi in zip(latitudes_deg, at):
        j = int(np.flatnonzero(curves.latitude_rad == float(phi))[0])
        produced = ft.produce_at_target(state, inputs, float(phi), curves.geopotential_m2s2[:, j], along[:, j],
                                        estimate.label_pressure_Pa, p_b, field, float(nl.datum_isobar_Pa))
        out[float(phi_deg)] = produced
    return {"profiles": out, "passes": len(state.passes) if hasattr(state, "passes") else None,
            "anchor_latitude_deg": float(anchors[0].latitude_planetocentric_deg), "mesh_nodes": int(mesh.latitude_rad.size)}


def section(control, latitudes_deg=None, work: Path = None) -> xr.Dataset:
    """The section on `latitudes_deg` (the control file's grid by default)."""
    latitudes = np.asarray(grid(control) if latitudes_deg is None else latitudes_deg, dtype="float64")
    work = work or control["here"] / "work"
    with xr.open_datatree(control["anchor_path"]) as anchor:
        anchor_lat = float(anchor.to_dataset(inherit=False)["latitude_planetocentric_deg"])
    limit = control["limit_deg"]
    profiles, record = {}, {}
    for target, side in ((-limit, latitudes < anchor_lat), (limit, latitudes >= anchor_lat)):
        if not side.any():
            continue
        result = march(namelist(control, target, work), latitudes[side])
        profiles.update(result["profiles"])
        record["north" if target > 0 else "south"] = result["mesh_nodes"]
    shutil.rmtree(work, ignore_errors=True)
    return assemble(control, latitudes, profiles, anchor_lat, record)


def assemble(control, latitudes, profiles, anchor_lat, record) -> xr.Dataset:
    first = profiles[float(latitudes[0])]
    labels = np.asarray(first["label_pressure_Pa"], dtype="float64")

    def stack(key):
        rows = []
        for lat in latitudes:
            p = profiles[float(lat)]
            if not np.array_equal(np.asarray(p["label_pressure_Pa"]), labels):
                raise ValueError(f"the levels at {lat} deg differ from those at {latitudes[0]} deg")
            rows.append(np.asarray(p[key], dtype="float64"))
        return np.stack(rows)

    T, z = stack("temperature_K"), stack("height_above_anchor_isobar_m")
    lapse = -np.stack([np.gradient(t, h) for t, h in zip(T, z)]) * 1e3
    dims = ("latitude", "level")
    ds = xr.Dataset(
        {
            "temperature_K": (dims, T, {"units": "K", "long_name": "delivered temperature"}),
            "pressure_Pa": (dims, stack("pressure_Pa"), {"units": "Pa", "long_name": "delivered pressure"}),
            "pressure_identity_residual": (dims, stack("pressure_identity_residual"),
                                           {"units": "1", "long_name": "delivered pressure over its label, minus one"}),
            "density_kg_m3": (dims, stack("density_kg_m3"), {"units": "kg m-3", "long_name": "mass density"}),
            "number_density_m3": (dims, stack("number_density_m3"), {"units": "m-3", "long_name": "number density"}),
            "mean_molar_mass_kg_mol": (dims, stack("mean_molar_mass_kg_mol"), {"units": "kg mol-1", "long_name": "mean molar mass"}),
            "u_ms": (dims, stack("u_column_ms"), {"units": "m s-1", "long_name": "zonal wind of the run on the labels"}),
            "height_above_anchor_isobar_m": (dims, z, {"units": "m", "long_name": "height above the anchor isobar along the local vertical"}),
            "lapse_rate_K_km": (dims, lapse, {"units": "K km-1", "long_name": "-dT/dz along the local vertical, centered"}),
        },
        coords={"latitude_planetocentric_deg": (("latitude",), latitudes, {"units": "degree"}),
                "pressure_label_Pa": (("level",), labels, {"units": "Pa", "long_name": "the anchor's isobar labels"})},
    )
    ds.attrs.update({
        "title": f"{control['name']}: the {control['anchor_slug']} anchor transferred to every latitude under {control['base_run']}'s inputs (SPEC_20)",
        "anchor": control["anchor_slug"], "anchor_latitude_planetocentric_deg": anchor_lat,
        "base_run": control["base_run"], "created_by": TOOL, "casspian_git_commit": cio.git_commit(),
        "input_hashes": cio.input_hashes([control["path"], control["anchor_path"]], control["output"]),
        "mesh_nodes": str(record), "method": "two transfers, one toward each pole, produced at each latitude by forward.transfer.produce_at_target",
    })
    return ds


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog=TOOL, description="An anchor transferred to every latitude, as a section.")
    parser.add_argument("control", help="the section's control file")
    args = parser.parse_args(argv)
    control = read_control(args.control)
    ds = section(control)
    control["output"].parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(control["output"])
    print(f"wrote {control['output']}: {ds.sizes['latitude']} latitudes, {ds.sizes['level']} levels")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
