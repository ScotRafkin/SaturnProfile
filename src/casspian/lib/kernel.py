"""The kernels on the working mesh. SPEC_04 Step 3.

Three quantities live here, all on the `(phi_i, Phi_j)` nodes of `lib.mesh`:

* the shear kernel `S = 2 Omega_abs r (du/dZ)_R` of Eq. A15, with
  `(du/dZ)_R = sin(phi) (du/dr)_phi + (cos(phi) / r) (du/dphi)_r`;
* its vertical integral `I(phi, Phi) = int_0^Phi (S / g) dPhi'` of Eq. A16, from the reference
  surface at each latitude;
* the transfer kernel `K = S / g + (d ln(R_bar / m_bar) / dphi)_p` of Eq. A27, the second term
  from the run's composition on the isobar.

**The wind is never resampled by the mesh (SPEC_06).** `IsobarKernel` holds the smooth geometry
of the mesh, formed once per state; `shear_at` evaluates Eq. A15 at any point with the wind read
there from the wind file's interpolant and the geometry bilinear on the mesh. The transfer reads it
on the isobar, at the curve's label (`shear_on_isobar`); the tracing integrates it in each column
with the wind file's pressure nodes as breakpoints of a two-point Gauss rule (`ColumnIntegral`).
Every change of slope the wind file can hold sits at those nodes, so the shear is used at the wind
file's resolution whatever the mesh spacing. The node kernel `shear_kernel` and its trapezoid
`shear_integral` remain, for the record and for the synthetic tests.

**Where the derivatives come from.** `u` is never differentiated on the mesh. The wind file is
the only thing that knows how `u` varies, and `lib.windfield` returns its interpolant's own
partials in the file's coordinates, `(phi, ln p)`, linear between nodes (decision L). Those are
in the file's coordinates, at fixed pressure, and A15 needs them at fixed radius, so they are
converted with the isobar map's slopes:

    (du/dphi)_r = (du/dphi)_p + (du/dln p)_phi (dln p/dphi)_r
    (du/dr)_phi = (du/dln p)_phi (dln p/dr)_phi

and only those slopes are taken by centered differences on the mesh, where the pressure map is
smooth. A wind file with no shear has `(du/dln p)_phi = 0` at every node, so both mesh differences
multiply zero and `S` reduces to `2 Omega_abs cos(phi) u'(phi)` exactly, with no difference taken.

Pure functions, NumPy in and out, no file access. Every angle is in radians.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from casspian.lib.gravity import omega_abs
from casspian.lib.reduction import mean_over_species

__all__ = ["isobar_slopes", "shear_kernel", "shear_integral", "composition_term",
           "transfer_kernel", "IsobarKernel", "isobar_kernel", "shear_at", "shear_on_isobar",
           "ColumnIntegral", "column_integral"]


def isobar_slopes(mesh, radius_m, pressure_Pa, g_radial_ms2):
    """`((dln p/dphi)_r, (dln p/dr)_phi)` on the nodes, by centered differences on the mesh.

    The mesh's own coordinates are `(phi, Phi)`, and both slopes are wanted in `(phi, r)`, so the
    mesh differences are changed over. Along a column `dr/dPhi = 1 / g` (Step 2), so

        (dln p/dr)_phi  = g (dln p/dPhi)_phi
        (dln p/dphi)_r  = (dln p/dphi)_Phi - g (dr/dphi)_Phi (dln p/dPhi)_phi

    the second because moving along a surface of constant `r` costs
    `(dPhi/dphi)_r = -(dr/dphi)_Phi / (dr/dPhi)_phi`. Each array is on the `(latitude,
    geopotential)` nodes.
    """
    ln_p = np.log(np.asarray(pressure_Pa, dtype="float64"))
    radius = np.asarray(radius_m, dtype="float64")
    g = np.asarray(g_radial_ms2, dtype="float64")
    if not (ln_p.shape == radius.shape == g.shape == mesh.shape):
        raise ValueError(
            f"the mesh is {mesh.shape} and the node arrays are {radius.shape}, {ln_p.shape} "
            f"and {g.shape}"
        )
    d_ln_p_d_phi = mesh.d_dlatitude(ln_p)
    d_r_d_phi = mesh.d_dlatitude(radius)
    d_ln_p_d_Phi = np.gradient(ln_p, mesh.geopotential_m2s2, axis=1, edge_order=2)
    return (d_ln_p_d_phi - g * d_r_d_phi * d_ln_p_d_Phi, g * d_ln_p_d_Phi)


def shear_kernel(mesh, radius_m, pressure_Pa, g_radial_ms2, u_ms, field, Omega):
    """`S` on the nodes, Eq. A15. SPEC_04 Step 3 deliverable 1.

    `field` is the run's `lib.windfield.WindField`, which supplies the only derivatives of `u`
    taken anywhere: its interpolant's partials at the node's own `(phi, p)`. `Omega` is the
    rotation rate of the frame the wind is declared in; `Omega_abs` is formed from `u` at the
    node by Eq. A2.

    Returns `(S, du_dZ)` so that a caller can check the axial derivative on its own, which the
    synthetic states of the acceptance do.
    """
    latitude = np.asarray(mesh.latitude_rad, dtype="float64")[:, None]
    radius = np.asarray(radius_m, dtype="float64")
    pressure = np.asarray(pressure_Pa, dtype="float64")
    u = np.asarray(u_ms, dtype="float64")

    d_u_d_phi_p, d_u_d_ln_p = field.wind_derivatives(np.broadcast_to(latitude, radius.shape),
                                                     pressure)
    d_ln_p_d_phi_r, d_ln_p_d_r = isobar_slopes(mesh, radius, pressure, g_radial_ms2)
    d_u_d_phi_r = d_u_d_phi_p + d_u_d_ln_p * d_ln_p_d_phi_r
    d_u_d_r_phi = d_u_d_ln_p * d_ln_p_d_r

    d_u_d_Z = np.sin(latitude) * d_u_d_r_phi + (np.cos(latitude) / radius) * d_u_d_phi_r
    return 2.0 * omega_abs(u, radius, latitude, Omega) * radius * d_u_d_Z, d_u_d_Z


def shear_integral(mesh, s_over_g):
    """`I` on the nodes, Eq. A16: the trapezoid of `S / g` in `Phi` from the `Phi = 0` node.

    Taken outward from the gauge in both directions, as `lib.geopotential.geopotential` takes the
    geopotential, so `I` is exactly zero on the `Phi = 0` node of every column and no node near it
    carries the round-off of a long cumulative sum.
    """
    values = np.asarray(s_over_g, dtype="float64")
    if values.shape != mesh.shape:
        raise ValueError(f"the mesh is {mesh.shape} and S / g is {values.shape}")
    Phi = mesh.geopotential_m2s2
    origin = mesh.gauge_geopotential_index
    increments = 0.5 * (values[:, :-1] + values[:, 1:]) * np.diff(Phi)[None, :]
    out = np.empty_like(values)
    out[:, origin] = 0.0
    out[:, origin + 1:] = np.cumsum(increments[:, origin:], axis=1)
    out[:, :origin] = -np.cumsum(increments[:, :origin][:, ::-1], axis=1)[:, ::-1]
    return out


def composition_term(composition, label_pressure_Pa):
    """`(d ln(R_bar / m_bar) / dphi)_p` at each isobar label. SPEC_04 Step 3 deliverable 3.

    Kind C is a field on `(level, latitude)` with pressure the vertical coordinate (decision O),
    so `R_bar` and `m_bar` are formed on the whole field, interpolated in `ln p` to each label,
    and differenced in latitude on the file's own grid by centered differences. Returns
    `(latitude_rad, term)` with `term` on `(label, latitude)`.

    A file whose columns are identical gives exactly zero, which is what the transfer composition
    is: the difference of equal numbers, not a small number.
    """
    root = composition.to_dataset(inherit=False)
    species = composition["species"].to_dataset(inherit=False)
    names = [str(name) for name in species["species_name"].values]
    x = np.stack([np.asarray(root[f"x_{name}"].values, dtype="float64") for name in names])
    R_i = np.asarray(species["refractivity_per_molecule_m3"].values, dtype="float64")
    M_i = np.asarray(species["molar_mass_kg_mol"].values, dtype="float64")
    R_bar = mean_over_species(x, R_i)
    m_bar = mean_over_species(x, M_i)

    latitude = np.radians(np.asarray(root["latitude_planetocentric_deg"].values, dtype="float64"))
    pressure = np.asarray(root["pressure_Pa"].values, dtype="float64")
    order = np.argsort(latitude)
    latitude = latitude[order]
    ratio = np.log(R_bar / m_bar)[:, order]

    rising = np.argsort(np.log(pressure))
    labels = np.atleast_1d(np.asarray(label_pressure_Pa, dtype="float64"))
    on_labels = np.stack([
        np.interp(np.log(labels), np.log(pressure)[rising], ratio[rising, j])
        for j in range(latitude.size)], axis=1)

    if latitude.size < 3:
        return latitude, np.zeros_like(on_labels)
    return latitude, _three_point_derivative(on_labels, latitude)


def _three_point_derivative(values, nodes):
    """The three-point derivative for unequal spacing, written as differences.

    The same formula `np.gradient` applies, rearranged so that every term carries a difference of
    neighbouring values rather than a weighted sum of the values themselves. The two forms are
    algebraically identical and differ only in floating point, and the difference matters here:
    `np.gradient` given a coordinate array takes its non-uniform path, whose coefficients
    `a + b + c` do not cancel to zero, so a field that is constant in latitude comes back with a
    derivative of about 1e-12 rather than zero. SPEC_04 Step 3 asks the composition term to be
    zero where the columns are identical, and a difference of equal numbers is exactly zero.

    Second order in the interior and at both ends. `values` is differenced along its last axis.
    """
    values = np.asarray(values, dtype="float64")
    nodes = np.asarray(nodes, dtype="float64")
    h = np.diff(nodes)
    h1, h2 = h[:-1], h[1:]
    forward = np.diff(values, axis=-1)
    out = np.empty_like(values)
    out[..., 1:-1] = ((h1 ** 2 * forward[..., 1:] + h2 ** 2 * forward[..., :-1])
                      / (h1 * h2 * (h1 + h2)))
    # The ends take the one-sided second-order form, the first divided difference corrected by
    # the second, which is again a difference of differences.
    second_low = (forward[..., 1] / h[1] - forward[..., 0] / h[0]) / (h[0] + h[1])
    out[..., 0] = forward[..., 0] / h[0] - h[0] * second_low
    second_high = (forward[..., -1] / h[-1] - forward[..., -2] / h[-2]) / (h[-1] + h[-2])
    out[..., -1] = forward[..., -1] / h[-1] + h[-1] * second_high
    return out


def transfer_kernel(mesh, s_over_g, phi_rad, geopotential_m2s2,
                    composition_latitude_rad=None, composition_slope=None):
    """`K` at a point of a curve, Eq. A27. SPEC_04 Step 3 deliverable 4.

    `S / g` is taken bilinearly on the mesh, since the curve does not pass through its nodes, and
    the composition term linearly in latitude on the composition file's grid, since it lives
    there and on the isobar label rather than on the mesh. Omitting the composition arguments
    gives the shear term alone, which is what a run whose composition is uniform in latitude
    needs and what it would compute.
    """
    phi = np.asarray(phi_rad, dtype="float64")
    Phi = np.asarray(geopotential_m2s2, dtype="float64")
    phi, Phi = np.broadcast_arrays(phi, Phi)
    values = np.asarray(s_over_g, dtype="float64")

    def cell(nodes, point):
        index = np.clip(np.searchsorted(nodes, point, side="right") - 1, 0, nodes.size - 2)
        return index, (point - nodes[index]) / (nodes[index + 1] - nodes[index])

    i, a = cell(mesh.latitude_rad, phi)
    j, b = cell(mesh.geopotential_m2s2, Phi)
    shear = (values[i, j] * (1.0 - a) * (1.0 - b) + values[i + 1, j] * a * (1.0 - b)
             + values[i, j + 1] * (1.0 - a) * b + values[i + 1, j + 1] * a * b)
    if composition_slope is None:
        return shear
    return shear + np.interp(phi, np.asarray(composition_latitude_rad, dtype="float64"),
                             np.asarray(composition_slope, dtype="float64"))


@dataclass(frozen=True)
class IsobarKernel:
    """What the shear term on an isobar needs, formed once per state. SPEC_06 Step 1.

    The mesh and, on its nodes, the radius, the radial gravity and the two isobar slopes of
    `isobar_slopes`, which are smooth fields of the pressure map; the run's `WindField`, which is
    read at the curve's point and label; and `Omega`, the rotation rate of the wind's frame.
    """

    mesh: object
    radius_m: np.ndarray
    g_radial_ms2: np.ndarray
    d_ln_p_d_phi_r: np.ndarray
    d_ln_p_d_r: np.ndarray
    field: object
    Omega: float


def isobar_kernel(mesh, radius_m, pressure_Pa, g_radial_ms2, field, Omega) -> IsobarKernel:
    """The `IsobarKernel` of a state: the slopes formed once on the nodes, as `shear_kernel` forms
    them."""
    radius = np.asarray(radius_m, dtype="float64")
    g = np.asarray(g_radial_ms2, dtype="float64")
    d_ln_p_d_phi_r, d_ln_p_d_r = isobar_slopes(mesh, radius, pressure_Pa, g)
    return IsobarKernel(mesh=mesh, radius_m=radius, g_radial_ms2=g,
                        d_ln_p_d_phi_r=d_ln_p_d_phi_r, d_ln_p_d_r=d_ln_p_d_r, field=field,
                        Omega=float(Omega))


def shear_on_isobar(kernel: IsobarKernel, phi_rad, geopotential_m2s2, label_pressure_Pa):
    """`S / g` at curve points `(phi, Phi)` on the isobar `p_label`, Eq. A15. SPEC_06 Step 1.

    `shear_at` with the wind read at the curve's label rather than at the pressure the map gives
    there (deliverable 5). A `(level,)` label against `(level, latitude)` points is given as
    `label[:, None]`.
    """
    return shear_at(kernel, phi_rad, geopotential_m2s2, label_pressure_Pa)


def shear_at(kernel: IsobarKernel, phi_rad, geopotential_m2s2, pressure_Pa):
    """`S / g` at points `(phi, Phi)` with the wind read at pressure `p`, Eq. A15. SPEC_06.

    The wind's part is read exactly at `(phi, p)` from the wind file's interpolant: `u`,
    `(du/dphi)_p` and `(du/dln p)_phi` (decision L; at a wind file node the cell toward higher
    pressure and toward the north, `WindField.wind_derivatives`' convention). The geometry, `r`,
    `g` and the two isobar slopes, is read bilinearly on the mesh at `(phi, Phi)`. The rest is
    `shear_kernel`'s arithmetic. The three arguments broadcast.
    """
    phi, Phi, p_label = np.broadcast_arrays(np.asarray(phi_rad, dtype="float64"),
                                            np.asarray(geopotential_m2s2, dtype="float64"),
                                            np.asarray(pressure_Pa, dtype="float64"))
    mesh = kernel.mesh
    radius = transfer_kernel(mesh, kernel.radius_m, phi, Phi)
    g = transfer_kernel(mesh, kernel.g_radial_ms2, phi, Phi)
    d_ln_p_d_phi_r = transfer_kernel(mesh, kernel.d_ln_p_d_phi_r, phi, Phi)
    d_ln_p_d_r = transfer_kernel(mesh, kernel.d_ln_p_d_r, phi, Phi)

    u = kernel.field.wind_at(phi, p_label)
    d_u_d_phi_p, d_u_d_ln_p = kernel.field.wind_derivatives(phi, p_label)
    d_u_d_phi_r = d_u_d_phi_p + d_u_d_ln_p * d_ln_p_d_phi_r
    d_u_d_r_phi = d_u_d_ln_p * d_ln_p_d_r

    d_u_d_Z = np.sin(phi) * d_u_d_r_phi + (np.cos(phi) / radius) * d_u_d_phi_r
    return 2.0 * omega_abs(u, radius, phi, kernel.Omega) * radius * d_u_d_Z / g


# The two-point Gauss rule on a piece `[a, b]`: the points `mid -/+ half / sqrt(3)`, each weighted
# `half`. Exact for cubics, so it stays adequate if the wind's interpolant becomes cubic.
_GAUSS_OFFSET = 1.0 / np.sqrt(3.0)


def _piecewise(kernel, isobar_map, index, start, stop):
    """`int S/g dPhi` in column `index` from `start` to `stop` (arrays), piece by piece.

    The breakpoints are the geopotentials where the column crosses a wind file pressure node
    strictly between the two ends, found by inverting the column's pressure map; on each piece a
    two-point Gauss rule, with `S/g` from `shear_at` at the pressure the map gives at the point.
    """
    start, stop = np.broadcast_arrays(np.asarray(start, dtype="float64"),
                                      np.asarray(stop, dtype="float64"))
    phi = float(kernel.mesh.latitude_rad[index])
    nodes = kernel.field.ln_pressure
    ln_a, ln_b = isobar_map.ln_p_at(index, start), isobar_map.ln_p_at(index, stop)
    low, high = np.minimum(ln_a, ln_b), np.maximum(ln_a, ln_b)
    first = np.searchsorted(nodes, low, side="right")
    count = np.searchsorted(nodes, high, side="left") - first
    most = int(count.max()) if count.size else 0
    inner = np.empty((most,) + start.shape)
    for m in range(most):
        at = np.clip(first + m, 0, nodes.size - 1)
        inner[m] = np.where(m < count, isobar_map.geopotential_at(index, nodes[at]), stop)
    # Padded with `stop`, sorted, and reversed where the integral runs downward, so the edges
    # run from `start` to `stop` in order and the padding gives pieces of zero length.
    inner = np.sort(inner, axis=0)
    inner = np.where(start > stop, inner[::-1], inner)
    edges = np.concatenate([start[None], inner, stop[None]])
    mid, half = 0.5 * (edges[1:] + edges[:-1]), 0.5 * (edges[1:] - edges[:-1])
    points = np.concatenate([mid - half * _GAUSS_OFFSET, mid + half * _GAUSS_OFFSET])
    values = shear_at(kernel, phi, points, np.exp(isobar_map.ln_p_at(index, points)))
    pieces = half * (values[:mid.shape[0]] + values[mid.shape[0]:])
    return pieces.sum(axis=0)


@dataclass(frozen=True)
class ColumnIntegral:
    """`I(phi, Phi)`, Eq. A16, from the shear read at the wind file's resolution. SPEC_06 v0.3.

    `at_nodes` is `I` on the mesh's nodes, zero on the `Phi = 0` node of every column as
    `shear_integral` makes it, each cell integrated by `_piecewise`. Between nodes in a column,
    `I` is the node below plus the partial cell up to the point, integrated the same way; between
    mesh latitudes, linear in latitude between the two columns. `isobar_map` is the state's map:
    `ln_p_at(index, Phi)` and its inverse `geopotential_at(index, ln_p)`.
    """

    kernel: IsobarKernel
    isobar_map: object
    at_nodes: np.ndarray

    def in_column(self, index, geopotential_m2s2):
        """`I` in column `index` at any geopotentials."""
        Phi = np.asarray(geopotential_m2s2, dtype="float64")
        nodes = self.kernel.mesh.geopotential_m2s2
        j = np.clip(np.searchsorted(nodes, Phi, side="right") - 1, 0, nodes.size - 2)
        return self.at_nodes[index, j] + _piecewise(self.kernel, self.isobar_map, index,
                                                    nodes[j], Phi)

    def at(self, phi_rad, geopotential_m2s2):
        """`I` at one latitude and any geopotentials; linear in latitude between columns."""
        latitude = self.kernel.mesh.latitude_rad
        phi = float(phi_rad)
        i = int(np.clip(np.searchsorted(latitude, phi, side="right") - 1, 0, latitude.size - 2))
        a = (phi - latitude[i]) / (latitude[i + 1] - latitude[i])
        if a == 0.0:
            return self.in_column(i, geopotential_m2s2)
        if a == 1.0:
            return self.in_column(i + 1, geopotential_m2s2)
        return ((1.0 - a) * self.in_column(i, geopotential_m2s2)
                + a * self.in_column(i + 1, geopotential_m2s2))


def column_integral(kernel: IsobarKernel, isobar_map) -> ColumnIntegral:
    """The state's `ColumnIntegral`: every column's cells integrated piece by piece, then summed
    outward from the `Phi = 0` node in both directions, as `shear_integral` sums them."""
    mesh = kernel.mesh
    Phi = mesh.geopotential_m2s2
    cells = np.stack([_piecewise(kernel, isobar_map, i, Phi[:-1], Phi[1:])
                      for i in range(mesh.latitude_rad.size)])
    origin = mesh.gauge_geopotential_index
    out = np.empty(mesh.shape)
    out[:, origin] = 0.0
    out[:, origin + 1:] = np.cumsum(cells[:, origin:], axis=1)
    out[:, :origin] = -np.cumsum(cells[:, :origin][:, ::-1], axis=1)[:, ::-1]
    return ColumnIntegral(kernel=kernel, isobar_map=isobar_map, at_nodes=out)
