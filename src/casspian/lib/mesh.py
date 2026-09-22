"""The working mesh and the columns on it. SPEC_04 Step 2.

The transfer is carried out on a mesh of planetocentric latitude and geopotential, `(phi_i,
Phi_j)`, and on the column at each latitude node: the radius, the height along the local vertical,
and the gravity the two are integrated with. Pure functions and one small dataclass, NumPy in and
out, no file access and no Saturn quantity of its own; every angle is in radians.

**The mesh is built from the data, not declared** (SPEC_04 v0.9). The namelist gives two spacings
and nothing else: the latitude nodes span the anchors, the gauge latitude and the target with one
spacing of margin each side, and the geopotential nodes span the anchors' own `Phi_k` with one
spacing beyond each. Nothing about the extent is checked, because there is nothing a user declared
to check it against. When Step 4's tracing finds a curve reaching an edge, `extend` adds whole
spacings on that side and the columns are rebuilt on the new nodes.

**The columns.** At each latitude node the radius follows `dr/dPhi = 1 / g` from the reference
surface at `Phi = 0`, outward in both directions, and the height along the local vertical follows
`dz/dPhi = 1 / |g_eff|` on the same integration (decision B). `g` is positive inward and `Phi`
increases upward (`lib.gravity`, `lib.geopotential`), so both slopes are positive and no hemisphere
needs a sign of its own. The two differ by `1 / cos psi`, which is the whole of the difference
between a radial altitude and a field-line altitude.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from casspian.lib.gravity import g_eff_vector

__all__ = ["Mesh", "Column", "build_mesh", "column", "build_columns", "SIDES"]

#: The sides `extend` names. The specification says "extended on that side by what the curve
#: needs"; a curve `Phi_k(phi)` can run off either end of either axis, so all four are named.
SIDES = ("south", "north", "below", "above")


@dataclass(frozen=True)
class Mesh:
    """The `(phi_i, Phi_j)` nodes, with the latitudes that had to be exact nodes recorded."""

    latitude_rad: np.ndarray
    geopotential_m2s2: np.ndarray
    latitude_spacing_rad: float
    geopotential_spacing_m2s2: float
    #: The anchor latitudes, the gauge latitude and the target, in the order given. Each is a
    #: node of `latitude_rad` exactly, and the mesh carries them so that a caller can find them
    #: again without matching floats of its own.
    required_latitude_rad: tuple

    @property
    def shape(self):
        return self.latitude_rad.size, self.geopotential_m2s2.size

    @property
    def gauge_geopotential_index(self) -> int:
        """The node where `Phi = 0`, which the build guarantees exists."""
        return int(np.flatnonzero(self.geopotential_m2s2 == 0.0)[0])

    def latitude_index(self, phi_rad) -> int:
        """The node index of a latitude that is one of the mesh's own nodes."""
        found = np.flatnonzero(self.latitude_rad == float(phi_rad))
        if found.size == 0:
            raise ValueError(
                f"{np.degrees(float(phi_rad))} degrees is not a node of this mesh; the anchors, "
                "the gauge latitude and the target are its exact nodes"
            )
        return int(found[0])

    def d_dlatitude(self, values, axis=0):
        """The centered difference in latitude, three-point for unequal spacing.

        Inserting the anchors, the gauge latitude and the target leaves cells shorter than the
        declared spacing on either side of them, so a difference that assumed uniform spacing
        would be first order there. `np.gradient` against the node coordinates is the three-point
        formula this deliverable asks for, second order on unequal spacing in the interior and,
        with `edge_order=2`, second order one-sided at the two edges as well. The edge order is
        given explicitly because `np.gradient` defaults to 1, which is first order one-sided and
        is not exact even on a quadratic: measured on a quadratic over nodes like these, the
        default leaves 1.5e-1 at each end against 1.1e-14 in the interior. The margin keeps the
        edges away from any traced curve, but a first-order edge would still enter the mesh-wide
        finite-difference fields Step 3 forms.
        """
        return np.gradient(np.asarray(values, dtype="float64"), self.latitude_rad, axis=axis,
                           edge_order=2)

    def extend(self, side: str, amount: float) -> "Mesh":
        """A mesh with whole spacings added on one side, enough to cover `amount`.

        `amount` is what a traced curve needs beyond the present edge, in radians on the latitude
        axis and in m2/s2 on the geopotential axis; it is rounded up to whole spacings, so the
        existing nodes are kept exactly and the new ones continue the same uniform run.
        """
        if side not in SIDES:
            raise ValueError(f"unknown side {side!r}; the sides are {SIDES}")
        if not np.isfinite(amount) or amount <= 0.0:
            raise ValueError(f"extend needs a positive amount; got {amount!r}")
        if side in ("south", "north"):
            nodes, spacing = self.latitude_rad, self.latitude_spacing_rad
        else:
            nodes, spacing = self.geopotential_m2s2, self.geopotential_spacing_m2s2
        steps = int(np.ceil(float(amount) / spacing))
        if side in ("south", "below"):
            added = nodes[0] - spacing * np.arange(steps, 0, -1)
            nodes = np.concatenate([added, nodes])
        else:
            added = nodes[-1] + spacing * np.arange(1, steps + 1)
            nodes = np.concatenate([nodes, added])
        if side in ("south", "north"):
            return replace(self, latitude_rad=nodes)
        return replace(self, geopotential_m2s2=nodes)


@dataclass(frozen=True)
class Column:
    """One latitude's column: the two altitudes and the gravity they were integrated with."""

    latitude_rad: float
    geopotential_m2s2: np.ndarray
    reference_radius_m: float
    radius_m: np.ndarray
    z_local_vertical_m: np.ndarray
    g_radial_ms2: np.ndarray
    G_phi_ms2: np.ndarray
    g_magnitude_ms2: np.ndarray
    psi_rad: np.ndarray
    u_ms: np.ndarray


def _spacing(name, value):
    value = float(value)
    if not np.isfinite(value) or value <= 0.0:
        raise ValueError(f"{name} must be a positive number; got {value!r}")
    return value


def build_mesh(required_latitude_rad, anchor_geopotential_m2s2, latitude_spacing_rad,
               geopotential_spacing_m2s2) -> Mesh:
    """The working mesh. SPEC_04 Step 2 deliverable 1.

    `required_latitude_rad` is every latitude that must be an exact node: the anchors, the gauge
    latitude (fixed before the mesh, decision K) and the target. `anchor_geopotential_m2s2` is
    every anchor's `Phi_k` from Step 1 deliverable 3, in one array or several.

    Latitude: uniform at the spacing from one spacing south of the southernmost required latitude
    to one spacing north of the northernmost, with every required latitude inserted exactly, so
    the run is uniform except in the two cells each insertion splits. Geopotential: whole multiples
    of the spacing, so `Phi = 0` is a node by construction, from one spacing below the lowest
    anchor `Phi_k` to one spacing above the highest.
    """
    latitude_spacing = _spacing("latitude_spacing_rad", latitude_spacing_rad)
    geopotential_spacing = _spacing("geopotential_spacing_m2s2", geopotential_spacing_m2s2)

    required = np.unique(np.asarray(required_latitude_rad, dtype="float64").ravel())
    if required.size == 0:
        raise ValueError("the mesh needs at least one required latitude")
    low, high = float(required[0]) - latitude_spacing, float(required[-1]) + latitude_spacing
    steps = int(np.ceil((high - low) / latitude_spacing))
    uniform = low + latitude_spacing * np.arange(steps + 1)
    latitude = np.unique(np.concatenate([uniform, required]))

    # One anchor's Phi_k, or a list of them: M anchors each bring their own levels.
    per_anchor = (anchor_geopotential_m2s2 if isinstance(anchor_geopotential_m2s2, (list, tuple))
                  else [anchor_geopotential_m2s2])
    Phi = np.concatenate([np.atleast_1d(np.asarray(a, dtype="float64")).ravel()
                          for a in per_anchor])
    if Phi.size == 0:
        raise ValueError("the mesh needs at least one anchor's geopotential")
    below = int(np.floor(float(Phi.min()) / geopotential_spacing)) - 1
    above = int(np.ceil(float(Phi.max()) / geopotential_spacing)) + 1
    if not below <= 0 <= above:
        raise ValueError(
            f"the anchors' geopotential runs {Phi.min()} to {Phi.max()} m2/s2 and does not "
            "bracket the gauge, where Phi = 0; the gauge isobar must be inside the anchors"
        )
    geopotential = geopotential_spacing * np.arange(below, above + 1, dtype="float64")

    return Mesh(
        latitude_rad=latitude,
        geopotential_m2s2=geopotential,
        latitude_spacing_rad=latitude_spacing,
        geopotential_spacing_m2s2=geopotential_spacing,
        required_latitude_rad=tuple(float(v) for v in
                                    np.asarray(required_latitude_rad, dtype="float64").ravel()),
    )


def column(phi_c, geopotential_m2s2, reference_radius_m, u_of_geopotential,
           Omega, GM, J, degrees, R_norm) -> Column:
    """One column by RK4 on the geopotential nodes. SPEC_04 Step 2 deliverable 1.

    `dr/dPhi = 1 / g(r, phi, u)` and `dz/dPhi = 1 / |g_eff|(r, phi, u)`, integrated together from
    `(reference_radius_m, 0)` outward in both directions, so `Phi = 0` carries the reference
    radius exactly and `z = 0` there. The nodes need not be uniform and need not be the mesh's:
    the acceptance integrates an anchor's own `Phi_k`, which is how the column is compared with
    what the reduction tabulated.

    `u_of_geopotential(Phi)` returns the zonal wind in m/s at this latitude, which the caller
    forms from the wind field and the isobar map; it is called at the RK4 stages as well as the
    nodes, since the stages are where the scheme needs it. The `u` recorded on the nodes is that
    same callable at the nodes, which is what `wind_on_mesh` gives them.
    """
    Phi = np.asarray(geopotential_m2s2, dtype="float64")
    if Phi.ndim != 1 or Phi.size < 2:
        raise ValueError(f"the column needs at least two geopotential nodes; got {Phi.shape}")
    if not np.all(np.diff(Phi) > 0.0):
        raise ValueError("the geopotential nodes must be strictly increasing")
    start = np.flatnonzero(Phi == 0.0)
    if start.size == 0:
        raise ValueError("the geopotential nodes must contain Phi = 0, where the column starts")
    origin = int(start[0])
    phi = float(phi_c)
    constants = (Omega, GM, J, degrees, R_norm)

    def slopes(radius, value):
        """`(1/g, 1/|g_eff|)` at a radius and a geopotential, the two right-hand sides."""
        u = float(np.asarray(u_of_geopotential(value), dtype="float64").reshape(()))
        g, _, magnitude, _ = g_eff_vector(u, radius, phi, *constants)
        return 1.0 / float(g), 1.0 / float(magnitude)

    radius = np.empty(Phi.shape, dtype="float64")
    height = np.empty(Phi.shape, dtype="float64")
    radius[origin], height[origin] = float(reference_radius_m), 0.0

    # An initial value problem in Phi: each node needs the one before it, so this is a loop over
    # nodes and not whole-array arithmetic (SPEC_00 section 3.1, as `lib.geoid`'s march is).
    for direction in (1, -1):
        r, z = float(reference_radius_m), 0.0
        index = origin
        while 0 <= index + direction < Phi.size:
            here, there = float(Phi[index]), float(Phi[index + direction])
            h = there - here
            k1r, k1z = slopes(r, here)
            k2r, k2z = slopes(r + 0.5 * h * k1r, here + 0.5 * h)
            k3r, k3z = slopes(r + 0.5 * h * k2r, here + 0.5 * h)
            k4r, k4z = slopes(r + h * k3r, there)
            r = r + (h / 6.0) * (k1r + 2.0 * k2r + 2.0 * k3r + k4r)
            z = z + (h / 6.0) * (k1z + 2.0 * k2z + 2.0 * k3z + k4z)
            index += direction
            radius[index], height[index] = r, z

    u = np.asarray([float(np.asarray(u_of_geopotential(v), dtype="float64").reshape(()))
                    for v in Phi], dtype="float64")
    g, G_phi, magnitude, psi = g_eff_vector(u, radius, np.full(Phi.shape, phi), *constants)
    return Column(
        latitude_rad=phi,
        geopotential_m2s2=Phi,
        reference_radius_m=float(reference_radius_m),
        radius_m=radius,
        z_local_vertical_m=height,
        g_radial_ms2=g,
        G_phi_ms2=G_phi,
        g_magnitude_ms2=magnitude,
        psi_rad=psi,
        u_ms=u,
    )


def build_columns(mesh: Mesh, reference_radius_m, wind_of, Omega, GM, J, degrees, R_norm):
    """The column at every latitude node of a mesh, in node order.

    `reference_radius_m` is `r0` on the mesh's latitude nodes, the Step 1 surface through the
    first anchor. `wind_of(index)` returns that column's `u_of_geopotential` callable, which the
    caller builds from the wind field and the isobar map; the columns are rebuilt on every pass of
    the outer loop because that map changes.
    """
    r0 = np.asarray(reference_radius_m, dtype="float64")
    if r0.shape != mesh.latitude_rad.shape:
        raise ValueError(
            f"the reference radius has shape {r0.shape} and the mesh has "
            f"{mesh.latitude_rad.shape} latitude nodes"
        )
    return [column(float(mesh.latitude_rad[i]), mesh.geopotential_m2s2, float(r0[i]),
                   wind_of(i), Omega, GM, J, degrees, R_norm)
            for i in range(mesh.latitude_rad.size)]
