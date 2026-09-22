"""The run's wind on the model's geometry. SPEC_04 Step 1 deliverable 2.

Kind W gives `u_total` on a grid of planetocentric latitude and pressure along the local
vertical. The model needs `u` at a point of its own geometry, which is a latitude and the
pressure the isobar map gives that point, so this module is the one place where the file's grid
is read off. Decision L fixes the rule: linear between the file's nodes, and every derivative
the transfer uses is the interpolant's own, not a separate fit. The same rule registered the
anchor, in `refrac.anchor.wind_of_latitude`, which is the reference level of this field.

Nothing here is extended or declared. The wind tool's fit, its parameterization and its prose
provenance are the tool's business (SPEC_04 section 1, decision A); what the file's grid does
not cover, this refuses.

Every latitude here is in radians, as SPEC_00 section 3.1 requires inside `lib`; pressures are
in Pa and the interpolation in pressure is linear in `ln p`.
"""

import numpy as np

LATITUDE_NAME = "latitude_planetocentric_deg"
PRESSURE_NAME = "pressure_Pa"
TOTAL_NAME = "u_total_ms"


class WindField:
    """`u_total` on the file's grid, with the reading rule of decision L.

    Built once from a kind W dataset and then asked for values, because the outer loop of
    SPEC_00 section 7.2 asks for the whole mesh once per pass and the file's grid does not
    change between passes.

    * `wind_at(phi, p)`: `u_total`, linear in latitude and linear in `ln p`.
    * `wind_derivatives(phi, p)`: `(du/dphi)_p` and `(du/dln p)_phi` of that same interpolant.
    * `wind_on_mesh(latitude, pressure)`: `u` at every node of a mesh column set.

    All three refuse a point outside the grid rather than extrapolating.
    """

    def __init__(self, wind):
        latitude = np.radians(np.asarray(wind[LATITUDE_NAME].values, dtype="float64"))
        pressure = np.asarray(wind[PRESSURE_NAME].values, dtype="float64")
        u = np.asarray(wind[TOTAL_NAME].values, dtype="float64")
        if u.shape != (latitude.size, pressure.size):
            raise ValueError(
                f"{TOTAL_NAME} has shape {u.shape}; kind W carries it on "
                f"(latitude, pressure) = {(latitude.size, pressure.size)}"
            )
        # The file's own order is not guaranteed ascending in either axis, and every index
        # below assumes it is. The file is sorted here once rather than at every call.
        in_latitude, in_pressure = np.argsort(latitude), np.argsort(pressure)
        self.latitude_rad = latitude[in_latitude]
        self.pressure_Pa = pressure[in_pressure]
        self.ln_pressure = np.log(self.pressure_Pa)
        self.u_total_ms = u[np.ix_(in_latitude, in_pressure)]
        self.reference_pressure_Pa = float(wind["reference_level_pressure_Pa"])

    # -- the grid ----------------------------------------------------------------------

    @property
    def latitude_bounds_rad(self):
        """The grid's own extent in latitude. No coverage attribute is read."""
        return float(self.latitude_rad[0]), float(self.latitude_rad[-1])

    @property
    def pressure_bounds_Pa(self):
        """The grid's own extent in pressure. No coverage attribute is read."""
        return float(self.pressure_Pa[0]), float(self.pressure_Pa[-1])

    def _cell(self, phi, p):
        """The bracketing cell and the two interpolation weights, refusing outside the grid.

        The index rule is `searchsorted(..., side="right") - 1`, so a point that falls exactly on
        a node takes the interval to its north in latitude and to higher pressure in pressure,
        which is what SPEC_04 Step 1 deliverable 2 fixes for the derivatives. At the last node of
        either axis there is no such interval and the previous one is used, which is the only one
        there; the value is the node's own either way.
        """
        phi = np.asarray(phi, dtype="float64")
        p = np.asarray(p, dtype="float64")
        phi, p = np.broadcast_arrays(phi, p)
        self._refuse_outside(phi, p)
        ln_p = np.log(p)

        def bracket(nodes, value):
            index = np.searchsorted(nodes, value, side="right") - 1
            index = np.clip(index, 0, nodes.size - 2)
            width = nodes[index + 1] - nodes[index]
            return index, (value - nodes[index]) / width

        j, weight_phi = bracket(self.latitude_rad, phi)
        k, weight_p = bracket(self.ln_pressure, ln_p)
        return j, k, weight_phi, weight_p

    def _refuse_outside(self, phi, p):
        """Refuse a point the file's grid does not reach. Decision N: extrapolation is refused."""
        low_phi, high_phi = self.latitude_bounds_rad
        low_p, high_p = self.pressure_bounds_Pa
        outside = (phi < low_phi) | (phi > high_phi)
        if outside.any():
            worst = float(np.degrees(phi[outside].flat[0]))
            raise ValueError(
                f"the wind is asked for at planetocentric latitude {worst} deg, outside the "
                f"file's grid [{np.degrees(low_phi)}, {np.degrees(high_phi)}] deg; the field is "
                "not extended (SPEC_04 Step 1 deliverable 2)"
            )
        outside = (p < low_p) | (p > high_p)
        if outside.any():
            worst = float(p[outside].flat[0])
            raise ValueError(
                f"the wind is asked for at {worst} Pa, outside the file's grid "
                f"[{low_p}, {high_p}] Pa; the field is not extended "
                "(SPEC_04 Step 1 deliverable 2)"
            )

    # -- the field ---------------------------------------------------------------------

    def wind_at(self, phi, p):
        """`u_total` at planetocentric latitude `phi` (rad) and pressure `p` (Pa), in m/s.

        Linear between the file's latitude nodes and linear in `ln p` between its pressure
        nodes (decision L). `phi` and `p` broadcast against each other and the result has the
        broadcast shape.
        """
        j, k, weight_phi, weight_p = self._cell(phi, p)
        u = self.u_total_ms
        south = u[j, k] * (1.0 - weight_p) + u[j, k + 1] * weight_p
        north = u[j + 1, k] * (1.0 - weight_p) + u[j + 1, k + 1] * weight_p
        return south * (1.0 - weight_phi) + north * weight_phi

    def wind_derivatives(self, phi, p):
        """`((du/dphi)_p, (du/dln p)_phi)` of the interpolant `wind_at` uses.

        Each is the interpolant's own derivative and so is constant across the cell in the
        direction it differences: `(du/dphi)_p` is the difference of the two `ln p` interpolated
        columns over the latitude interval, and `(du/dln p)_phi` the difference of the two
        latitude interpolated rows over the `ln p` interval. Units are m/s per radian and m/s
        per unit `ln p`. At a node the cell taken is the one to the north, or to higher pressure,
        as `_cell` states.
        """
        j, k, weight_phi, weight_p = self._cell(phi, p)
        u = self.u_total_ms
        south = u[j, k] * (1.0 - weight_p) + u[j, k + 1] * weight_p
        north = u[j + 1, k] * (1.0 - weight_p) + u[j + 1, k + 1] * weight_p
        shallow = u[j, k] * (1.0 - weight_phi) + u[j + 1, k] * weight_phi
        deep = u[j, k + 1] * (1.0 - weight_phi) + u[j + 1, k + 1] * weight_phi
        d_phi = self.latitude_rad[j + 1] - self.latitude_rad[j]
        d_ln_p = self.ln_pressure[k + 1] - self.ln_pressure[k]
        return (north - south) / d_phi, (deep - shallow) / d_ln_p

    def wind_on_mesh(self, latitude_rad, pressure_Pa):
        """`u` at every node of a set of columns, from each node's latitude and pressure.

        `latitude_rad` is one latitude per column and `pressure_Pa` is the pressure the isobar
        map gives every node, with the columns along its first axis. The mesh object itself
        arrives at Step 2; what this needs of it is the column latitudes and the map, which is
        what it takes.
        """
        latitude = np.asarray(latitude_rad, dtype="float64")
        pressure = np.asarray(pressure_Pa, dtype="float64")
        if latitude.ndim != 1 or pressure.ndim < 1 or pressure.shape[0] != latitude.size:
            raise ValueError(
                f"the pressure map has shape {pressure.shape} and there are {latitude.size} "
                "column latitudes; the map's first axis is the columns"
            )
        shape = (latitude.size,) + (1,) * (pressure.ndim - 1)
        return self.wind_at(latitude.reshape(shape), pressure)

    def reference_wind(self, phi):
        """`u_total` at the file's reference level, the callable the reduction anchored on.

        Provided so that the transfer and the reduction read one field: this is `wind_at` at
        `reference_level_pressure_Pa` and reproduces `refrac.anchor.wind_of_latitude` on a file
        whose reference level is a pressure node, which kind W requires.
        """
        return self.wind_at(phi, np.full(np.shape(phi), self.reference_pressure_Pa))
