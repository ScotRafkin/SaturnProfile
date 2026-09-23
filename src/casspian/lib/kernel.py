"""The kernels on the working mesh. SPEC_04 Step 3.

Three quantities live here, all on the `(phi_i, Phi_j)` nodes of `lib.mesh`:

* the shear kernel `S = 2 Omega_abs r (du/dZ)_R` of Eq. A15, with
  `(du/dZ)_R = sin(phi) (du/dr)_phi + (cos(phi) / r) (du/dphi)_r`;
* its vertical integral `I(phi, Phi) = int_0^Phi (S / g) dPhi'` of Eq. A16, from the reference
  surface at each latitude;
* the transfer kernel `K = S / g + (d ln(R_bar / m_bar) / dphi)_p` of Eq. A27, the second term
  from the run's composition on the isobar.

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

import numpy as np

from casspian.lib.gravity import omega_abs
from casspian.lib.reduction import mean_over_species

__all__ = ["isobar_slopes", "shear_kernel", "shear_integral", "composition_term",
           "transfer_kernel"]


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
