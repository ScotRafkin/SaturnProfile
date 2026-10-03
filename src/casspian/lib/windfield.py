"""The run's wind on the model's geometry. SPEC_04 Step 1 deliverable 2, SPEC_07 Step 1.

Kind W gives `u_total` on a grid of planetocentric latitude and pressure along the local
vertical. The model needs `u` at a point of its own geometry, which is a latitude and the
pressure the isobar map gives that point, so this module is the one place where the file's grid
is read off. Decision L2 (SPEC_07) fixes the rule: a shape preserving cubic (PCHIP, Fritsch and
Carlson, in the form `scipy.interpolate.PchipInterpolator` uses) between the file's nodes, first
in latitude on every pressure row, then in `ln p` at the point. Its first derivatives are
continuous, it passes through every node, and it never leaves a cell's corner values. Every
derivative the transfer uses is the interpolant's own, not a separate fit. The same rule
registers the anchor: `refrac.anchor.wind_of_latitude` is this field's `reference_wind`.

Nothing here is extended or declared. The wind tool's fit, its parameterization and its prose
provenance are the tool's business (SPEC_04 section 1, decision A); what the file's grid does
not cover, this refuses.

Every latitude here is in radians, as SPEC_00 section 3.1 requires inside `lib`; pressures are
in Pa and the interpolation in pressure is in `ln p`.
"""

import numpy as np
from scipy.interpolate import PchipInterpolator

LATITUDE_NAME = "latitude_planetocentric_deg"
PRESSURE_NAME = "pressure_Pa"
TOTAL_NAME = "u_total_ms"


class WindField:
    """`u_total` on the file's grid, with the reading rule of decision L2.

    Built once from a kind W dataset and then asked for values, because the outer loop of
    SPEC_00 section 7.2 asks for the whole mesh once per pass and the file's grid does not
    change between passes. The latitude step's coefficients are formed here, once.

    * `wind_at(phi, p)`: `u_total`, PCHIP in latitude, then PCHIP in `ln p`.
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
        if self.pressure_Pa.size < 3 or self.latitude_rad.size < 3:
            raise ValueError(
                f"the wind's grid is {self.latitude_rad.size} latitudes by {self.pressure_Pa.size} "
                "pressures; decision L2 needs at least three nodes on each axis"
            )
        # Step 1 of decision L2: every pressure row interpolated in latitude, once.
        self._rows = PchipInterpolator(self.latitude_rad, self.u_total_ms, axis=0)
        self._rows_dphi = self._rows.derivative()

    # -- the grid ----------------------------------------------------------------------

    @property
    def latitude_bounds_rad(self):
        """The grid's own extent in latitude. No coverage attribute is read."""
        return float(self.latitude_rad[0]), float(self.latitude_rad[-1])

    @property
    def pressure_bounds_Pa(self):
        """The grid's own extent in pressure. No coverage attribute is read."""
        return float(self.pressure_Pa[0]), float(self.pressure_Pa[-1])

    def _points(self, phi, p):
        """The points broadcast, refused outside the grid, and the `ln p` interval of each.

        The index rule is `searchsorted(..., side="right") - 1`, so a point that falls exactly on
        a node takes the interval to higher pressure; at the last node there is no such interval
        and the previous one is used, the node's own value being returned either way.
        """
        phi = np.asarray(phi, dtype="float64")
        p = np.asarray(p, dtype="float64")
        phi, p = np.broadcast_arrays(phi, p)
        self._refuse_outside(phi, p)
        x = np.log(p)
        k = np.clip(np.searchsorted(self.ln_pressure, x, side="right") - 1,
                    0, self.ln_pressure.size - 2)
        return phi, x, k

    def _rows_at(self, phi):
        """`v_k(phi)` and `v_k'(phi)` on every pressure row, at each point: `(..., rows)`.

        Evaluated once per distinct latitude, since the model asks for whole columns at one
        latitude; each value is the same arithmetic whatever else is asked with it. At the last
        latitude node the row is the file's own, which the polynomial of the last interval
        reproduces only to rounding.
        """
        unique, inverse = np.unique(phi, return_inverse=True)
        values, slopes = self._rows(unique), self._rows_dphi(unique)
        last = unique == self.latitude_rad[-1]
        values[last] = self.u_total_ms[-1]
        inverse = inverse.reshape(phi.shape)
        return values[inverse], slopes[inverse]

    def _read(self, phi, p, derivatives):
        """Step 2 of decision L2: PCHIP in `ln p` on the four rows about each point's interval.

        In local form, the piece as the node value plus powers of `s = x - x_k`, with scipy's
        slopes and coefficients (`PchipInterpolator._find_derivatives`, `_edge_case` and
        `CubicHermiteSpline`), so a node returns the file's value bit for bit and a column with
        no vertical shear the same value at every pressure. With `derivatives`, also
        `(du/dln p)` and `(du/dphi)`, the second the same expression differentiated in latitude
        through the rows' derivatives and the slopes' own (the tangent of the harmonic mean,
        zero on its zero branch).
        """
        phi, x, k = self._points(phi, p)
        values, slopes = self._rows_at(phi)
        nodes = self.ln_pressure
        last = nodes.size - 1
        index = np.clip(k[..., None] + np.arange(-1, 3), 0, last)
        v = np.take_along_axis(values, index, axis=-1)
        dv = np.take_along_axis(slopes, index, axis=-1)
        h = nodes[np.minimum(index[..., 1:], last)] - nodes[index[..., :-1]]
        safe = np.where(h == 0.0, 1.0, h)
        d = (v[..., 1:] - v[..., :-1]) / safe
        dd = (dv[..., 1:] - dv[..., :-1]) / safe
        # Columns 0, 1, 2 of h, d, dd are the intervals k - 1, k, k + 1; at the grid's ends the
        # clipped ones are unused, since the edge rule takes their place.
        m_k, dm_k = _node_slope(k == 0, h[..., 1], h[..., 2], d[..., 1], d[..., 2], dd[..., 1],
                                dd[..., 2], h[..., 0], d[..., 0], dd[..., 0])
        m_n, dm_n = _node_slope(k + 1 == last, h[..., 1], h[..., 0], d[..., 1], d[..., 0],
                                dd[..., 1], dd[..., 0], h[..., 2], d[..., 2], dd[..., 2],
                                reverse=True)
        width = h[..., 1]
        s = x - nodes[k]
        u = _local(v[..., 1], m_k, m_n, d[..., 1], width, s)
        # At the last node the piece of the last interval reaches it only to rounding.
        u = np.where(x == nodes[last], v[..., 2], u)
        if not derivatives:
            return u
        du_dx = _local_slope(m_k, m_n, d[..., 1], width, s)
        du_dphi = _local(dv[..., 1], dm_k, dm_n, dd[..., 1], width, s)
        return u, du_dphi, du_dx

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

        Decision L2: PCHIP in latitude, then PCHIP in `ln p`. `phi` and `p` broadcast against
        each other and the result has the broadcast shape.
        """
        return self._read(phi, p, derivatives=False)

    def wind_derivatives(self, phi, p):
        """`((du/dphi)_p, (du/dln p)_phi)` of the interpolant `wind_at` uses.

        Both continuous across every node (decision L2): `(du/dln p)` is the derivative of the
        cubic in `ln p`, and `(du/dphi)` that of the same expression in latitude, carried
        through the rows' PCHIP derivatives and the derivatives of the slopes. Units are m/s per
        radian and m/s per unit `ln p`.
        """
        _, du_dphi, du_dx = self._read(phi, p, derivatives=True)
        return du_dphi, du_dx

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


def _harmonic(h_prev, h_next, d_prev, d_next, dd_prev, dd_next):
    """The interior PCHIP slope at a node and its derivative in latitude.

    scipy's weighted harmonic mean of the two secants, `w1 = 2 h_next + h_prev`,
    `w2 = h_next + 2 h_prev`, and zero where the secants differ in sign or either is zero. The
    derivative is the tangent of the same formula, `m^2 (w1 d_prev' / d_prev^2 + w2 d_next' /
    d_next^2) / (w1 + w2)`, and zero on the zero branch.
    """
    zero = (np.sign(d_next) != np.sign(d_prev)) | (d_next == 0) | (d_prev == 0)
    w1 = 2 * h_next + h_prev
    w2 = h_next + 2 * h_prev
    with np.errstate(divide="ignore", invalid="ignore"):
        mean = (w1 / d_prev + w2 / d_next) / (w1 + w2)
        m = 1.0 / mean
        dm = m * m * (w1 * dd_prev / (d_prev * d_prev) + w2 * dd_next / (d_next * d_next)) / (w1 + w2)
    return np.where(zero, 0.0, m), np.where(zero, 0.0, dm)


def _edge(h0, h1, m0, m1, dm0, dm1):
    """scipy's three point slope at the grid's first or last node, shape preserving, and its
    derivative in latitude on the branch taken."""
    d = ((2 * h0 + h1) * m0 - h0 * m1) / (h0 + h1)
    dd = ((2 * h0 + h1) * dm0 - h0 * dm1) / (h0 + h1)
    flat = np.sign(d) != np.sign(m0)
    clamp = ~flat & (np.sign(m0) != np.sign(m1)) & (np.abs(d) > 3.0 * np.abs(m0))
    return (np.where(flat, 0.0, np.where(clamp, 3.0 * m0, d)),
            np.where(flat, 0.0, np.where(clamp, 3.0 * dm0, dd)))


def _node_slope(at_edge, h_here, h_beyond, d_here, d_beyond, dd_here, dd_beyond,
                h_other, d_other, dd_other, reverse=False):
    """The slope at one end of a point's interval: the edge rule where that end is the grid's
    first (or, `reverse`, last) node, the harmonic mean elsewhere. `here` is the point's own
    interval, `beyond` the next one away from the other end, `other` the one on the far side."""
    edge = _edge(h_here, h_beyond, d_here, d_beyond, dd_here, dd_beyond)
    if reverse:
        inner = _harmonic(h_here, h_other, d_here, d_other, dd_here, dd_other)
    else:
        inner = _harmonic(h_other, h_here, d_other, d_here, dd_other, dd_here)
    return np.where(at_edge, edge[0], inner[0]), np.where(at_edge, edge[1], inner[1])


def _local(y, m0, m1, slope, width, s):
    """A cubic Hermite piece in local form, scipy's coefficients, at `s = x - x_k`."""
    t = (m0 + m1 - 2 * slope) / width
    c0, c1 = t / width, (slope - m0) / width - t
    return y + m0 * s + c1 * (s * s) + c0 * (s * s * s)


def _local_slope(m0, m1, slope, width, s):
    """The derivative of `_local` in `x`."""
    t = (m0 + m1 - 2 * slope) / width
    c0, c1 = t / width, (slope - m0) / width - t
    return m0 + 2 * c1 * s + 3 * c0 * (s * s)
