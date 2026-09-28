"""The wind fields the SPEC_04 Step 5 acceptance needs beside the closure file's.

`cylinder_extended` rebuilds decision P's cylinder-extended wind by the Step 2 construction, and
`sheared` writes the sheared synthetic file of the Step 4 expected values and reads it back through
the kind W reader. Both are taken from `tests/step04_4/accept_step04_4.py`, where they were
measured against decisions P and Q and reported in REPORT_04_step4 checks 4 and 11, so that this
suite builds its own states and stands alone. Neither touches a registered input.
"""

import math

import numpy as np

from casspian.forward import transfer as tr
from casspian.lib import geoid
from casspian.lib import io as cio
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf


def cylinder_extended(inputs, anchor, closure_field, p_b, constants, gauge_isobar_Pa,
                      latitude_spacing_rad, extra_latitudes=()):
    """Decision P's cylinder-extended wind, per hemisphere, as a fixed point.

    The wind at the reference level carried along cylinders of constant `s = r cos(phi)`, the
    cylinder radius taken from the model's own geometry, which depends on the wind, which is why it
    is a fixed point. Returns the field and what the construction gives on the anchor's column, the
    quantity the Step 2 acceptance compares with decision P's stated values.
    """
    phi_a = math.radians(float(anchor.latitude_planetocentric_deg))
    gauge_level = int(anchor.gauge_level_index)
    p_tab = np.asarray(anchor.label_pressure_Pa, dtype="float64")
    lat_file = np.radians(np.asarray(inputs.wind["latitude_planetocentric_deg"].values,
                                    dtype="float64"))
    p_file = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
    u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
    p_ref = float(inputs.wind["reference_level_pressure_Pa"])
    ref_col = int(np.flatnonzero(p_file == p_ref)[0])
    Phi_closure = tr.place(anchor, inputs, closure_field, p_b).geopotential_m2s2

    rising = np.argsort(np.log(p_tab))
    _x = np.log(p_tab)[rising]
    _y = Phi_closure[rising]
    _flat = tr.IsobarMap(geopotential_m2s2=np.sort(Phi_closure)[None, :],
                         ln_pressure=np.log(p_tab)[np.argsort(Phi_closure)][None, :])

    def Phi_of_p(p):
        """The anchor's flat-isobar map, linear in `ln p`, continued at the slope of its last
        interval, which is the rule `forward.transfer.IsobarMap` applies in the other direction."""
        lp = np.log(np.asarray(p, dtype="float64"))
        low = (_y[1] - _y[0]) / (_x[1] - _x[0])
        high = (_y[-1] - _y[-2]) / (_x[-1] - _x[-2])
        return np.where(lp < _x[0], _y[0] + low * (lp - _x[0]),
                        np.where(lp > _x[-1], _y[-1] + high * (lp - _x[-1]), np.interp(lp, _x, _y)))

    def on_file_grid(Phi):
        """`p(Phi)` of the flat map, held inside the wind file's own pressure range.

        The round trip through the map is not exact to the bit at its ends, so the pressure can land
        one ulp outside the grid and `lib.windfield` refuses it, correctly. The excursion is
        asserted to be at that scale, so this is a clamp of one ulp and nothing structural.
        """
        p = np.exp(_flat.ln_p_at(0, Phi))
        low, high = float(p_file[0]), float(p_file[-1])
        excursion = max(float(np.max(low - np.minimum(p, low))) / low,
                        float(np.max(np.maximum(p, high) - high)) / high)
        assert excursion < 1e-12, f"the flat map leaves the file's grid by {excursion:.3e} relative"
        return np.clip(p, low, high)

    fine = np.arange(float(lat_file[0]), float(lat_file[-1]) + 0.5 * latitude_spacing_rad,
                     latitude_spacing_rad)
    lat_inv = np.unique(np.clip(
        np.concatenate([fine, lat_file, [phi_a], np.asarray(extra_latitudes, dtype="float64")]),
        float(lat_file[0]), float(lat_file[-1])))
    on_file = np.searchsorted(lat_inv, lat_file)

    Phi_k, field_now, passes = Phi_closure.copy(), closure_field, 0
    while True:
        Phi_p = Phi_of_p(p_file)
        nodes = np.unique(np.concatenate([Phi_p, Phi_k, [0.0]]))
        at_pressure = np.searchsorted(nodes, Phi_p)
        r0_inv = geoid.through_anchor(
            lat_inv, phi_a, float(np.asarray(anchor.radius_m)[gauge_level]),
            lambda x: field_now.wind_at(np.asarray(x, dtype="float64"), gauge_isobar_Pa),
            *constants).radius
        at_level = np.searchsorted(nodes, Phi_k)
        radius_inv = np.empty((lat_inv.size, p_file.size))
        anchor_column = None
        for i, latitude in enumerate(lat_inv):
            built_column = lm.column(
                latitude, nodes, float(r0_inv[i]),
                (lambda Phi, la=latitude: field_now.wind_at(la, on_file_grid(Phi))),
                *constants)
            radius_inv[i] = built_column.radius_m[at_pressure]
            if latitude == phi_a:
                anchor_column = built_column
        s_reference = radius_inv[:, ref_col] * np.cos(lat_inv)

        def branch(mask):
            s, q = s_reference[mask], lat_inv[mask]
            rising_s = np.argsort(s)
            return s[rising_s], q[rising_s]

        s_north, phi_north = branch(lat_inv >= 0.0)
        s_south, phi_south = branch(lat_inv <= 0.0)
        s_grid = radius_inv[on_file] * np.cos(lat_file)[:, None]
        in_north = lat_file >= 0.0
        phi_star = np.empty_like(s_grid)
        phi_star[in_north] = np.interp(s_grid[in_north], s_north, phi_north)
        phi_star[~in_north] = np.interp(s_grid[~in_north], s_south, phi_south)
        u_cylinder = np.interp(phi_star, lat_file, u_reference_file)
        u_cylinder[0], u_cylinder[-1] = 0.0, 0.0
        # The wind on the construction itself, at the anchor's own levels: the anchor column's
        # radius there turned into a cylinder radius and inverted, with no file grid in between.
        on_construction = np.interp(
            np.interp(anchor_column.radius_m[at_level] * math.cos(phi_a), s_north, phi_north),
            lat_file, u_reference_file)
        change = float(np.max(np.abs(u_cylinder - field_now.u_total_ms)))
        passes += 1
        built = inputs.wind.copy(deep=True)
        built["u_total_ms"] = (("latitude_planetocentric", "pressure"), u_cylinder,
                              dict(inputs.wind["u_total_ms"].attrs))
        field_now = wf.WindField(built)
        Phi_k = tr.place(anchor, inputs, field_now, p_b).geopotential_m2s2
        if change < 1.0e-6 or passes >= 10:
            break

    return {
        "field": field_now,
        "passes": passes,
        "latitudes": int(lat_inv.size),
        "on_construction": np.asarray(on_construction, dtype="float64"),
        "u_read_back": np.asarray(field_now.wind_at(np.full(p_tab.shape, phi_a), p_tab),
                                  dtype="float64"),
        "final_change_ms": change,
    }


def sheared(inputs, path, beta=0.1):
    """The sheared synthetic wind of the Step 4 expected values, written and read back.

    `u_total(phi, p) = u_reference(phi) [1 + beta ln(p_ref / p)]` above the reference level and
    `u_reference` below it, written as a kind W file in its three parts so that the reader checks
    the sum identity and the poles.
    """
    lat_file = np.asarray(inputs.wind["latitude_planetocentric_deg"].values, dtype="float64")
    p_file = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
    u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
    p_ref = float(inputs.wind["reference_level_pressure_Pa"])
    grid = np.broadcast_to(u_reference_file[:, None], (lat_file.size, p_file.size))
    factor = np.where(p_file < p_ref, 1.0 + beta * np.log(p_ref / p_file), 1.0)
    u_sheared = grid * factor[None, :]
    written = inputs.wind.copy(deep=True)
    written["u_total_ms"] = (("latitude_planetocentric", "pressure"), u_sheared,
                             dict(inputs.wind["u_total_ms"].attrs))
    written["u_shear_ms"] = (("latitude_planetocentric", "pressure"), u_sheared - grid,
                             dict(inputs.wind["u_shear_ms"].attrs))
    written.attrs["title"] = "Sheared synthetic wind for the SPEC_04 Step 5 acceptance"
    out = cio.write(path, written, "wind", created_by="accept_step04_5")
    handle = cio.read(out, "wind")
    try:
        back = handle.load()
    finally:
        handle.close()
    sum_identity = float(np.max(np.abs(
        np.asarray(back["u_total_ms"].values)
        - np.asarray(back["u_reference_ms"].values)[:, None]
        - np.asarray(back["u_shear_ms"].values))))
    return {"field": wf.WindField(back), "path": out, "sum_identity_ms": sum_identity,
            "beta": float(beta), "reference_pressure_Pa": p_ref}
