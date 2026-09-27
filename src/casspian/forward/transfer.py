"""The isobar tracing, the transfer of refractivity along it, and the outer loop. SPEC_04 Step 4.

Three things live here, in the order a run uses them:

* `trace`: the characteristic `dPhi/dphi = -I(phi, Phi)` of Eq. A24 and B6.2, by RK4 on the
  mesh's latitude nodes with `I` bilinear between them. Each level of an anchor is an isobar,
  and its curve is where that isobar sits at every latitude.
* `transfer`: Eq. A28, `ln N = ln N_k + int K dphi'`, the trapezoid of the transfer kernel
  along each curve on the same nodes.
* `outer_loop`: the fixed point of SPEC_00 section 7.2. The wind is data along the local
  vertical against pressure, so placing it on the model's geometry needs the isobar map, and
  the map is what the tracing produces. The loop starts from flat isobars, builds the columns
  and the kernels under the map it has, traces, forms the map the curves give, and repeats
  until the map stops moving.

**What the map is.** At each latitude node the map is `ln p` against `Phi` through that
latitude's knots: the isobar labels and the geopotential the curves put them at. On the first
pass the knots are the anchors' own productions at their own latitudes, the same at every
latitude, which is the barotropic guess. Between anchors the knots are those of the nearest
anchor in latitude, the rule the barotropic guess states, kept for the traced map so that the
two are one rule. Outside the knots the map is continued at the slope of its last interval, as
Step 2 decision 8 continues the flat map: a clamp there puts several mesh nodes at one
pressure, which is the fault SPEC_04 section 14 ruling 2 named.

**What is refused (decision N).** Two numerical failures this step names: a pair of curves that
cross, and the outer loop reaching `max_iterations` without converging. A curve that reaches the
mesh's `Phi` edge is not a failure: the trace stops, returns the side and the excess, and the
loop extends the mesh by what the curve needs and repeats the pass, recorded.

Pure functions and NumPy arrays; the only file access is through the loaded `inputs`. Every
angle is in radians.
"""

from __future__ import annotations

from dataclasses import dataclass
from types import MappingProxyType

import numpy as np

from casspian.forward import production as fp
from casspian.lib import geoid
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf

__all__ = ["Curves", "IsobarMap", "Placement", "TransferState", "boundary_pressure", "place",
           "barotropic_map", "trace", "transfer", "outer_loop"]


def _constants(inputs):
    """`(Omega, GM, J, degree, R_norm)`, the argument run `lib.gravity` and `lib.mesh` take."""
    gravity, rotation = inputs.gravity, inputs.rotation
    return (
        float(rotation["angular_rate_rad_s"]),
        float(gravity["GM_m3s2"]),
        np.asarray(gravity["J"].values, dtype="float64"),
        np.asarray(gravity["degree"].values),
        float(gravity["normalization_radius_m"]),
    )


def _interpolate_continued(x, nodes, values):
    """Linear interpolation continued at the slope of the last interval outside the nodes.

    `np.interp` clamps instead, which would put every point beyond the outermost knot at one
    value. `nodes` is strictly increasing.
    """
    x = np.asarray(x, dtype="float64")
    low = (values[1] - values[0]) / (nodes[1] - nodes[0])
    high = (values[-1] - values[-2]) / (nodes[-1] - nodes[-2])
    return np.where(
        x < nodes[0], values[0] + low * (x - nodes[0]),
        np.where(x > nodes[-1], values[-1] + high * (x - nodes[-1]),
                 np.interp(x, nodes, values)))


# ---------------------------------------------------------------------------------------------
# The anchors on the model's geometry
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Placement:
    """One anchor produced on its own levels under the run's wind. SPEC_04 Step 1 deliverable 3.

    `label_pressure_Pa` is the produced pressure of each level, the label of its isobar; the
    tabulated pressures are what the first pass supplies to the wind, and the produced ones what
    every pass after it does.
    """

    slug: str
    latitude_rad: float
    reference_radius_m: float
    geopotential_m2s2: np.ndarray
    label_pressure_Pa: np.ndarray
    ln_N: np.ndarray
    production: object


def boundary_pressure(anchors) -> float:
    """`p_b`, the label of the topmost isobar of the highest anchor (decision F).

    The pressure differential is exact, so the top isobar is one surface planet wide and `p_b` is
    one scalar for the run. `anchors` are `lib.control.LoadedAnchor` objects, whose labels are
    the tabulated pressures as they arrive.
    """
    return float(min(float(np.asarray(a.label_pressure_Pa, dtype="float64").min())
                     for a in anchors))


def place(anchor, inputs, field, p_b, label_pressure_Pa=None) -> Placement:
    """Produce one anchor at its own latitude under the run's wind along its column.

    `field` is the run's `lib.windfield.WindField`; the wind on the column is
    `wind_at(phi_i, p_k)` with `p_k` the labels supplied, the anchor's tabulated pressures on
    the first pass of the outer loop and its produced pressures after (Step 1 deliverable 3).
    """
    labels = (np.asarray(anchor.label_pressure_Pa, dtype="float64")
              if label_pressure_Pa is None else
              np.asarray(label_pressure_Pa, dtype="float64"))
    phi = float(np.radians(anchor.latitude_planetocentric_deg))
    u_column = np.asarray(field.wind_at(np.full(labels.shape, phi), labels), dtype="float64")
    production = fp.produce(fp.profile_from_anchor(anchor.tree), inputs,
                            anchor.gauge_level_index, float(p_b), u_column=u_column)
    return Placement(
        slug=anchor.slug,
        latitude_rad=phi,
        reference_radius_m=float(np.asarray(anchor.radius_m,
                                            dtype="float64")[anchor.gauge_level_index]),
        geopotential_m2s2=np.asarray(production.geopotential.geopotential_m2s2, dtype="float64"),
        label_pressure_Pa=np.asarray(production.pressure_Pa, dtype="float64"),
        ln_N=np.asarray(anchor.ln_N, dtype="float64"),
        production=production,
    )


# ---------------------------------------------------------------------------------------------
# The isobar map
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class IsobarMap:
    """`p(phi, Phi)` on a mesh: `ln p` against `Phi` through each latitude's knots.

    Both are one array per latitude node, ascending in `Phi`; a latitude's knots are the isobar
    labels and where the curves put them there, and between anchors they are the nearest
    anchor's. They are held per latitude rather than as one rectangle because two anchors need
    not carry the same number of levels.
    """

    geopotential_m2s2: tuple
    ln_pressure: tuple

    def ln_p_at(self, index: int, geopotential_m2s2):
        """`ln p` at one latitude node, for any `Phi`, continued outside the knots."""
        return _interpolate_continued(geopotential_m2s2, self.geopotential_m2s2[index],
                                      self.ln_pressure[index])

    def on_nodes(self, mesh) -> np.ndarray:
        """`ln p` at every `(phi_i, Phi_j)` node of `mesh`."""
        return np.stack([self.ln_p_at(i, mesh.geopotential_m2s2)
                         for i in range(mesh.latitude_rad.size)])


def _nearest_anchor(mesh, latitudes) -> np.ndarray:
    """The index of the anchor nearest in latitude to each mesh node; ties take the first."""
    distance = np.abs(np.asarray(latitudes, dtype="float64")[None, :]
                      - mesh.latitude_rad[:, None])
    return np.argmin(distance, axis=1)


def _knots(geopotential_m2s2, label_pressure_Pa):
    """One anchor's knots, sorted ascending in `Phi`."""
    Phi = np.asarray(geopotential_m2s2, dtype="float64")
    ln_p = np.log(np.asarray(label_pressure_Pa, dtype="float64"))
    rising = np.argsort(Phi)
    return Phi[rising], ln_p[rising]


def barotropic_map(mesh, placements) -> IsobarMap:
    """The first pass's map: flat isobars, each anchor's own production carried in latitude.

    Every latitude node takes the knots of the anchor nearest to it, so at one anchor the map is
    that anchor's production exactly.
    """
    which = _nearest_anchor(mesh, [p.latitude_rad for p in placements])
    knots = [_knots(p.geopotential_m2s2, p.label_pressure_Pa) for p in placements]
    return IsobarMap(
        geopotential_m2s2=tuple(knots[k][0] for k in which),
        ln_pressure=tuple(knots[k][1] for k in which),
    )


def _map_from_curves(mesh, curves, label_pressure_Pa) -> IsobarMap:
    """The map the traced curves give: at each latitude, the labels where the curves put them.

    `curves` is one `Curves` per anchor, each covering every latitude node of the mesh, and
    `label_pressure_Pa` the matching labels. Between anchors the nearest anchor's curves are
    taken, as the barotropic guess takes its knots.
    """
    which = _nearest_anchor(mesh, [c.start_latitude_rad for c in curves])
    Phi, ln_p = [None] * mesh.latitude_rad.size, [None] * mesh.latitude_rad.size
    for k, curve in enumerate(curves):
        if curve.latitude_rad.size != mesh.latitude_rad.size:
            raise ValueError(
                f"the curves of anchor {k} cover {curve.latitude_rad.size} latitudes and the "
                f"mesh has {mesh.latitude_rad.size}; the map is formed on the whole mesh"
            )
        # The curves cannot cross, so the order of the levels in `Phi` is the same at every
        # latitude and the anchor's own is taken once.
        rising = np.argsort(curve.geopotential_m2s2[:, 0])
        labels = np.log(np.asarray(label_pressure_Pa[k], dtype="float64"))[rising]
        for row in np.flatnonzero(which == k):
            Phi[row] = curve.geopotential_m2s2[rising, row]
            ln_p[row] = labels
    return IsobarMap(geopotential_m2s2=tuple(Phi), ln_pressure=tuple(ln_p))


# ---------------------------------------------------------------------------------------------
# The tracing
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Curves:
    """One anchor's isobars across a run of latitude nodes. SPEC_04 Step 4 deliverable 1.

    `latitude_rad` is ascending and `geopotential_m2s2` is `(level, latitude)` in the order the
    levels were given. `reached` is `None`, or `(side, excess)` naming the mesh edge a curve
    reached and how far past it the curve went, in m2/s2; the nodes beyond that point are not
    in the arrays.
    """

    start_latitude_rad: float
    latitude_rad: np.ndarray
    geopotential_m2s2: np.ndarray
    reached: object = None


def _refuse_if_crossing(Phi, order, latitude_rad):
    """Refuse where the curves are no longer in the order they started in."""
    bad = np.flatnonzero(np.diff(Phi[order]) <= 0.0)
    if bad.size == 0:
        return
    first = int(bad[0])
    raise ValueError(
        f"isobars {int(order[first])} and {int(order[first + 1])} cross at "
        f"{np.degrees(float(latitude_rad)):.6f} degrees; the wind field is not balanced there"
    )


def trace(mesh, shear_integral, geopotential_m2s2, phi_from, phi_to) -> Curves:
    """The isobars through `geopotential_m2s2` at `phi_from`, carried to `phi_to`.

    RK4 of `dPhi/dphi = -I(phi, Phi)` (Eq. A24, B6.2) on the mesh's own latitude nodes, with `I`
    bilinear between them. Both latitudes are nodes of the mesh. The levels may be given in any
    order and come back in it.

    A curve that leaves the mesh's geopotential range stops the trace, which returns what it has
    with the side and the excess in `reached`; the caller extends the mesh and repeats. A pair of
    curves that cross is refused: the isobars of one atmosphere cannot exchange places, and where
    the traced ones do, the wind field is not in balance with the state that produced them.
    """
    start = mesh.latitude_index(phi_from)
    stop = mesh.latitude_index(phi_to)
    Phi = np.asarray(geopotential_m2s2, dtype="float64").copy()
    if Phi.ndim != 1:
        raise ValueError(f"trace takes one geopotential per level; got shape {Phi.shape}")
    nodes = mesh.latitude_rad
    order = np.argsort(Phi)
    _refuse_if_crossing(Phi, order, nodes[start])

    step = 1 if stop >= start else -1
    low, high = float(mesh.geopotential_m2s2[0]), float(mesh.geopotential_m2s2[-1])
    kept = [Phi.copy()]
    latitudes = [float(nodes[start])]
    reached = None
    index = start
    while index != stop:
        here, there = float(nodes[index]), float(nodes[index + step])
        h = there - here
        # `transfer_kernel` with no composition is the bilinear value of a mesh field, which is
        # what the stages of the scheme need between the nodes.
        k1 = -lk.transfer_kernel(mesh, shear_integral, here, Phi)
        k2 = -lk.transfer_kernel(mesh, shear_integral, here + 0.5 * h, Phi + 0.5 * h * k1)
        k3 = -lk.transfer_kernel(mesh, shear_integral, here + 0.5 * h, Phi + 0.5 * h * k2)
        k4 = -lk.transfer_kernel(mesh, shear_integral, there, Phi + h * k3)
        Phi = Phi + (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
        index += step
        # The excess is the largest overshoot over the levels, on each side.
        below, above = float(np.max(low - Phi)), float(np.max(Phi - high))
        if below > 0.0 or above > 0.0:
            reached = ("below", below) if below >= above else ("above", above)
            break
        _refuse_if_crossing(Phi, order, nodes[index])
        kept.append(Phi.copy())
        latitudes.append(float(nodes[index]))

    latitude = np.asarray(latitudes, dtype="float64")
    curves = np.stack(kept, axis=1)
    if step < 0:
        latitude, curves = latitude[::-1], curves[:, ::-1]
    return Curves(start_latitude_rad=float(nodes[start]), latitude_rad=latitude,
                  geopotential_m2s2=curves, reached=reached)


def _joined(south: Curves, north: Curves) -> Curves:
    """One anchor's two traces, south and north of it, as one set of curves over both."""
    if south.reached is not None or north.reached is not None:
        return south if south.reached is not None else north
    latitude = np.concatenate([south.latitude_rad[:-1], north.latitude_rad])
    curves = np.concatenate([south.geopotential_m2s2[:, :-1], north.geopotential_m2s2], axis=1)
    return Curves(start_latitude_rad=south.start_latitude_rad, latitude_rad=latitude,
                  geopotential_m2s2=curves, reached=None)


# ---------------------------------------------------------------------------------------------
# The transfer
# ---------------------------------------------------------------------------------------------
def transfer(mesh, s_over_g, curves: Curves, ln_N, composition_latitude_rad=None,
             composition_slope=None) -> np.ndarray:
    """`ln N` along each curve, Eq. A28. SPEC_04 Step 4 deliverable 1.

    `ln_N` is the value at the curve's own latitude, one per level, and the result is
    `(level, latitude)` on the curve's nodes. `K = S/g + (d ln(R_bar/m_bar)/dphi)_p` is formed
    at each node of each curve by `lib.kernel.transfer_kernel` and integrated by the trapezoid in
    latitude from the start node outward, so the value there is `ln N` unchanged and the integral
    carries its own sign whichever way the curve runs.

    `composition_slope` is `(level, latitude)` on `composition_latitude_rad`, one row per level,
    as `lib.kernel.composition_term` returns it for that level's label. Omitted, `K` is the shear
    term alone, which is what a composition uniform in latitude gives.
    """
    ln_N = np.asarray(ln_N, dtype="float64")
    phi = curves.latitude_rad
    Phi = curves.geopotential_m2s2
    if ln_N.shape != (Phi.shape[0],):
        raise ValueError(f"{Phi.shape[0]} curves and {ln_N.shape} starting values")
    K = lk.transfer_kernel(mesh, s_over_g, phi[None, :], Phi)
    if composition_slope is not None:
        slope = np.atleast_2d(np.asarray(composition_slope, dtype="float64"))
        grid = np.asarray(composition_latitude_rad, dtype="float64")
        K = K + np.stack([np.interp(phi, grid, slope[k]) for k in range(slope.shape[0])])

    start = int(np.flatnonzero(phi == curves.start_latitude_rad)[0])
    increments = 0.5 * (K[:, :-1] + K[:, 1:]) * np.diff(phi)[None, :]
    out = np.empty_like(K)
    out[:, start] = 0.0
    if start + 1 < phi.size:
        out[:, start + 1:] = np.cumsum(increments[:, start:], axis=1)
    if start > 0:
        out[:, :start] = -np.cumsum(increments[:, :start][:, ::-1], axis=1)[:, ::-1]
    return ln_N[:, None] + out


# ---------------------------------------------------------------------------------------------
# The outer loop
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class TransferState:
    """What one converged pass of the outer loop leaves on the mesh."""

    mesh: object
    reference_radius_m: np.ndarray
    columns: tuple
    radius_m: np.ndarray
    pressure_Pa: np.ndarray
    g_radial_ms2: np.ndarray
    u_ms: np.ndarray
    s_over_g: np.ndarray
    shear_integral: np.ndarray
    isobar_map: IsobarMap
    placements: tuple
    curves: tuple
    record: MappingProxyType


def _stack(columns, name):
    return np.stack([getattr(c, name) for c in columns])


def outer_loop(inputs, anchors, *, gauge_latitude_rad, target_latitude_rad, gauge_isobar_Pa,
               p_b, latitude_spacing_rad, geopotential_spacing_m2s2, relative_tolerance_ln_p,
               max_iterations, field=None, mesh=None) -> TransferState:
    """The fixed point of the wind on the model's geometry. SPEC_04 Step 4 deliverable 1.

    One pass: place every anchor under the wind the map gives along its column, build the
    reference surface and the columns, form the kernels, trace every anchor's levels across the
    whole mesh, and read the map off the curves. The pass repeats until the largest `|d ln p|`
    on the mesh between the map that built it and the map it produced is below
    `relative_tolerance_ln_p`, and refuses at `max_iterations`, naming the residual.

    A file whose columns do not vary with pressure has a wind on the mesh that does not depend on
    the map, so the second pass reproduces the first exactly and the residual is zero; the code
    does not special-case it, and the record says how many passes were taken.

    `mesh` is the mesh to work on. Left out, the loop builds it from the anchors' levels and the
    latitudes that must be nodes, which is what a run does. Given, it is used as it stands and
    still grown when a curve or an anchor's levels need more: two runs that must be compared level
    for level have to be on one mesh, which is what the M = 2 identity test needs (SPEC_04 v0.15
    Step 4).
    """
    field = wf.WindField(inputs.wind) if field is None else field
    constants = _constants(inputs)
    required = np.unique(np.asarray(
        [float(np.radians(a.latitude_planetocentric_deg)) for a in anchors]
        + [float(gauge_latitude_rad), float(target_latitude_rad)], dtype="float64"))
    if mesh is not None:
        missing = [float(x) for x in required if not np.any(mesh.latitude_rad == x)]
        if missing:
            raise ValueError(
                f"the mesh given has no node at {np.degrees(missing)} degrees; the anchors, the "
                "gauge latitude and the target are its exact nodes (SPEC_04 Step 2)"
            )

    labels = [np.asarray(a.label_pressure_Pa, dtype="float64") for a in anchors]
    isobar_map, extensions = None, []
    passes = []
    for iteration in range(1, int(max_iterations) + 1):
        placements = tuple(place(a, inputs, field, p_b, labels[i])
                           for i, a in enumerate(anchors))
        if mesh is None:
            mesh = lm.build_mesh(required, [p.geopotential_m2s2 for p in placements],
                                 latitude_spacing_rad, geopotential_spacing_m2s2)
        else:
            mesh = _cover_levels(mesh, placements, extensions)
        if isobar_map is None:
            isobar_map = barotropic_map(mesh, placements)

        # Extending in geopotential leaves the latitude nodes alone, so the map's knots, which
        # are one row per latitude, stay as they are.
        #
        # The excess a stopped trace reports is the overshoot at the node it stopped on, and it
        # says nothing about how much further the curve climbs over the latitudes it has not
        # reached; extending by it alone makes the pass repeat once per node, and each repeat
        # rebuilds every column. So a second stop on the same side within a pass adds twice what
        # the one before it added. Both the excess and what was added are recorded.
        added = {"below": 0.0, "above": 0.0}
        while True:
            state = _geometry(inputs, mesh, isobar_map, placements, field, constants,
                              gauge_isobar_Pa)
            traced, reached = _trace_all(mesh, state["shear_integral"], placements)
            if reached is None:
                break
            side, excess = reached
            amount = max(float(excess) + geopotential_spacing_m2s2, 2.0 * added[side])
            added[side] = amount
            extensions.append({"pass": iteration, "side": side, "excess_m2s2": float(excess),
                               "added_m2s2": amount})
            mesh = mesh.extend(side, amount)

        new_map = _map_from_curves(mesh, traced, [p.label_pressure_Pa for p in placements])
        residual = float(np.max(np.abs(new_map.on_nodes(mesh) - isobar_map.on_nodes(mesh))))
        passes.append({"pass": iteration, "residual_ln_p": residual,
                       "latitudes": int(mesh.latitude_rad.size),
                       "geopotential_nodes": int(mesh.geopotential_m2s2.size)})
        if residual < float(relative_tolerance_ln_p):
            # The state returned is the one the curves were traced on, and `isobar_map` is the
            # map that built it; the map the curves give differs from it by less than the
            # tolerance, which is what convergence means here.
            return TransferState(
                mesh=mesh, isobar_map=isobar_map, placements=placements, curves=tuple(traced),
                record=MappingProxyType({
                    "passes": int(iteration),
                    "residual_ln_p": residual,
                    "history": tuple(MappingProxyType(p) for p in passes),
                    "extensions": tuple(MappingProxyType(e) for e in extensions),
                    "relative_tolerance_ln_p": float(relative_tolerance_ln_p),
                }),
                **state)
        isobar_map = new_map
        labels = [p.label_pressure_Pa for p in placements]

    raise ValueError(
        f"the outer loop did not converge in {max_iterations} passes: the largest |d ln p| on "
        f"the mesh is {passes[-1]['residual_ln_p']:.3e} against a tolerance of "
        f"{float(relative_tolerance_ln_p):.3e} (SPEC_04 Step 4 deliverable 1)"
    )


def _cover_levels(mesh, placements, extensions):
    """Extend the mesh where a pass's production moved an anchor's levels outside it."""
    Phi = np.concatenate([p.geopotential_m2s2 for p in placements])
    for side, excess in (("below", float(mesh.geopotential_m2s2[0] - Phi.min())),
                         ("above", float(Phi.max() - mesh.geopotential_m2s2[-1]))):
        if excess > 0.0:
            extensions.append({"pass": "levels", "side": side, "excess_m2s2": excess})
            mesh = mesh.extend(side, excess + mesh.geopotential_spacing_m2s2)
    return mesh


def _geometry(inputs, mesh, isobar_map, placements, field, constants, gauge_isobar_Pa):
    """The reference surface, the columns, the pressure on the nodes and the kernels."""
    first = placements[0]
    reference = geoid.through_anchor(
        mesh.latitude_rad, first.latitude_rad, first.reference_radius_m,
        lambda x: field.wind_at(np.asarray(x, dtype="float64"), gauge_isobar_Pa),
        *constants).radius
    columns = lm.build_columns(
        mesh, reference,
        (lambda i: (lambda Phi, index=i, latitude=float(mesh.latitude_rad[i]):
                    field.wind_at(latitude, np.exp(isobar_map.ln_p_at(index, Phi))))),
        *constants)
    radius = _stack(columns, "radius_m")
    g = _stack(columns, "g_radial_ms2")
    u = _stack(columns, "u_ms")
    pressure = np.exp(isobar_map.on_nodes(mesh))
    S, _ = lk.shear_kernel(mesh, radius, pressure, g, u, field, constants[0])
    s_over_g = S / g
    return {
        "reference_radius_m": reference,
        "columns": tuple(columns),
        "radius_m": radius,
        "pressure_Pa": pressure,
        "g_radial_ms2": g,
        "u_ms": u,
        "s_over_g": s_over_g,
        "shear_integral": lk.shear_integral(mesh, s_over_g),
    }


def _trace_all(mesh, shear_integral, placements):
    """Every anchor's levels traced to both edges of the mesh, or the first edge reached."""
    traced = []
    for placement in placements:
        south = trace(mesh, shear_integral, placement.geopotential_m2s2,
                      placement.latitude_rad, mesh.latitude_rad[0])
        north = trace(mesh, shear_integral, placement.geopotential_m2s2,
                      placement.latitude_rad, mesh.latitude_rad[-1])
        joined = _joined(south, north)
        if joined.reached is not None:
            return traced, joined.reached
        traced.append(joined)
    return traced, None
