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

import math

import numpy as np
import xarray as xr

from casspian.forward import estimate as fe
from casspian.forward import production as fp
from casspian.forward import propagate as fprop
from casspian.lib import geoid
from casspian.lib import control as ctl
from casspian.lib import io as cio
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf
from casspian.lib.schema import PROFILE_EPOCH_WITHOUT_DATE

#: SPEC_04 section 16 ruling 2: what `D_ij` cannot be read below. An anchor read back from a kind N
#: file is placed by the field-line integral over its own tabulated levels while the trace that
#: carried it is the characteristic, and the two differ by that much for no physical reason. The
#: product states it once so that a reader of `D_ij` has it beside them.
IDENTITY_FLOOR = 6.5e-8
IDENTITY_FLOOR_NOTE = (
    "two anchors whose isobars coincide disagree in D_ij at this level for no physical reason: an "
    "anchor read from a kind N file is placed by the field-line integral over its own tabulated "
    "levels (decision G) while the trace that carried it is the characteristic. Measured at Step 4 "
    "as 6.484e-08, mesh-independent, against the tracing's own identity of 3.553e-15. It is "
    "invisible against the anchors' declared uncertainties (2.3e-2 for Lindal); D_ij below it says "
    "nothing about the atmosphere (SPEC_04 section 16 ruling 2)."
)

__all__ = ["Arrived", "Curves", "IDENTITY_FLOOR", "IsobarMap", "Placement", "TransferResult",
           "TransferState",
           "barotropic_map", "boundary_pressure", "carry_to", "outer_loop", "place",
           "produce_at_target", "run", "trace", "transfer"]


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


# ---------------------------------------------------------------------------------------------
# The run: the production at the target, the altitude and the datum, the product. SPEC_04 Step 5
# ---------------------------------------------------------------------------------------------
@dataclass(frozen=True)
class Arrived:
    """One anchor carried to the gauge latitude, and its curves across the mesh."""

    anchor: object
    curves: Curves
    ln_N_along: np.ndarray
    arrival: object


@dataclass(frozen=True)
class TransferResult:
    """What a transfer run leaves: the product path, the state, the estimate and the figures."""

    product: object
    state: TransferState
    estimate: object
    arrived: tuple
    target: MappingProxyType
    figures: object


def carry_to(state: TransferState, inputs, phi_from, phi_to, geopotential_m2s2, ln_N,
             label_pressure_Pa):
    """Trace from one latitude to another and transfer `ln N` along the curves.

    The two halves of Eq. A24 and A28 in the order a run uses them, with the composition term on
    the labels the curves carry. Returns `(curves, column index of phi_to, ln N along the curves)`.
    """
    curves = trace(state.mesh, state.shear_integral, geopotential_m2s2, phi_from, phi_to)
    if curves.reached is not None:
        side, excess = curves.reached
        raise ValueError(
            f"a curve reached the mesh's {side} edge by {excess} m2/s2 while carrying from "
            f"{np.degrees(float(phi_from))} to {np.degrees(float(phi_to))} degrees; the mesh is "
            "grown inside the outer loop and this trace is after it (SPEC_04 Step 5)"
        )
    latitude_rad, slope = lk.composition_term(inputs.composition, label_pressure_Pa)
    along = transfer(state.mesh, state.s_over_g, curves, ln_N, latitude_rad, slope)
    return curves, int(np.flatnonzero(curves.latitude_rad == float(phi_to))[0]), along


def produce_at_target(state: TransferState, inputs, target_latitude_rad,
                      arrival_geopotential_m2s2, ln_N, label_pressure_Pa, p_b, field,
                      datum_isobar_Pa):
    """The production, the altitude and the datum at the target. Step 5 deliverables 1 and 2.

    The levels come back top down, which is the order `pressure_from_top` integrates in and the
    order kind `profile` stores. `R_bar` and `m_bar` are the run's composition at the target
    latitude on the isobar labels (section 15 ruling 4); `radius_m` and the height along the local
    vertical are the target column's, integrated on the arrival levels themselves rather than
    interpolated between the mesh's nodes, which is what SPEC_04 section 16 ruling 1 settled for
    the anchor written at Step 4 and holds for the same reason here.

    `arrival_geopotential_m2s2` is where the isobars arrive **at the target**, which is what the
    curves carry there and not what they left the gauge at. The distinction is the isobar shift,
    and it is the whole of the pressure identity: the refractivity is larger at 10 N by 1.1e-2 in
    `ln N`, the stack compresses by as much in the mass between isobars, and the hydrostatic
    integral returns each label. Entered with the gauge's geopotential instead, the produced
    pressure misses its label by exactly that 1.1e-2, which is how this was caught.
    """
    Phi_in = np.asarray(arrival_geopotential_m2s2, dtype="float64")
    order = np.argsort(Phi_in)[::-1]
    Phi = Phi_in[order]
    labels = np.asarray(label_pressure_Pa, dtype="float64")[order]
    N = np.exp(np.asarray(ln_N, dtype="float64")[order])

    index = state.mesh.latitude_index(target_latitude_rad)
    nodes = np.unique(np.concatenate([Phi, [0.0]]))
    column = lm.column(
        float(target_latitude_rad), nodes, float(state.reference_radius_m[index]),
        (lambda Phi_at, i=index, la=float(target_latitude_rad):
         field.wind_at(la, np.exp(state.isobar_map.ln_p_at(i, Phi_at)))),
        *_constants(inputs))
    at_level = np.searchsorted(nodes, Phi)
    radius = column.radius_m[at_level]
    z_lv = column.z_local_vertical_m[at_level]

    R_bar, m_bar = fp.mean_properties_on_labels(inputs.composition,
                                               math.degrees(float(target_latitude_rad)), labels)
    n, density, layer_mass, pressure, temperature = fp.produce_on_geopotential(
        N, Phi, R_bar, m_bar, p_b)

    datum_Phi = fp.datum_geopotential(pressure, Phi, datum_isobar_Pa)
    datum_z = float(np.interp(datum_Phi, column.geopotential_m2s2, column.z_local_vertical_m))
    datum_radius = float(np.interp(datum_Phi, column.geopotential_m2s2, column.radius_m))
    gauge = int(np.flatnonzero(Phi == 0.0)[0])
    u_column = np.asarray(field.wind_at(np.full(Phi.shape, float(target_latitude_rad)), labels),
                          dtype="float64")
    return MappingProxyType({
        "order": order,
        "geopotential_m2s2": Phi,
        "label_pressure_Pa": labels,
        "refractivity": N,
        "radius_m": radius,
        "height_above_anchor_isobar_m": z_lv,
        "altitude_m": z_lv - datum_z,
        "mean_refractivity_m3": R_bar,
        "mean_molar_mass_kg_mol": m_bar,
        "number_density_m3": n,
        "density_kg_m3": density,
        "layer_mass_Pa": layer_mass,
        "pressure_Pa": pressure,
        "temperature_K": temperature,
        "pressure_identity_residual": pressure / labels - 1.0,
        "u_column_ms": u_column,
        "psi_rad": column.psi_rad[at_level],
        "datum_geopotential_m2s2": datum_Phi,
        "datum_radius_m": datum_radius,
        "datum_altitude_reference_m": datum_z,
        "gauge_level_index": gauge,
        "boundary_pressure_Pa": float(p_b),
        "column": column,
    })


def _transfer_product(namelist, inputs, anchors, state, estimate, target_run, curves_to_target,
                      product):
    """Kind `profile` in transfer mode. SPEC_04 Step 5 deliverable 3.

    The root is the delivered profile at the target on the arrival levels, top down. What closure
    mode carries and this does not is the tabulated pair, which a transferred profile has nothing
    to be compared against; what this carries and closure does not is the isobar label of every
    level and the pressure identity against it, the gauge's own `C` and geopotential beside the
    target's, the altitude above the datum, and the four groups of the transfer.
    """
    level = ("level",)
    order = target_run["order"]
    gauge_index = target_run["gauge_level_index"]
    target_deg = float(namelist.target_latitude_deg)
    unstated = ", ".join(fp._declared_terms(anchors[0].tree.to_dataset(inherit=False))
                         + ["kernel_error", "propagation"])
    reason = ("not propagated: a first order companion carrying only the terms that pass through "
              "the transfer would be read as a total; the Monte Carlo wrapper supplies the "
              "product's uncertainty (SPEC_03 decision 5, SPEC_04 Step 5 deliverable 3)")

    def companion(units, of):
        return (level, np.full(order.size, np.nan),
                fp._attrs(units, f"uncertainty on {of}", "modeled", uncertainty_kind="1sigma",
                          uncertainty_method="not propagated",
                          uncertainty_terms_unstated=unstated, uncertainty_note=reason))

    def scalar_companion(units, of):
        """The same, for a scalar: kind profile carries a companion for every modeled variable."""
        return ((), np.nan,
                fp._attrs(units, f"uncertainty on {of}", "modeled", uncertainty_kind="1sigma",
                          uncertainty_method="not propagated",
                          uncertainty_terms_unstated=unstated, uncertainty_note=reason))

    variables = {
        "radius_m": (level, target_run["radius_m"],
                     fp._attrs("m", "absolute planetocentric radius of the level on the target's "
                               "column (B8 carries the radial distance)", "modeled",
                               positive="up")),
        "height_above_anchor_isobar_m": (
            level, target_run["height_above_anchor_isobar_m"],
            fp._attrs("m", "height above the gauge isobar along the local vertical at the target",
                      "modeled", positive="up")),
        "altitude_m": (level, target_run["altitude_m"],
                       fp._attrs("m", "delivered altitude above the datum isobar along the local "
                                 "vertical (decision B, B7.3 and B7.4)", "modeled",
                                 positive="up")),
        "altitude_uncertainty_m": companion("m", "altitude_m"),
        # Modeled here, not copied: at the target the radius and the height along the local vertical
        # come from the target's own column, so each carries the companion closure mode has no need
        # of (SPEC_03 decision 5, and the schema is what caught it).
        "radius_uncertainty_m": companion("m", "radius_m"),
        "height_above_anchor_isobar_uncertainty_m": companion(
            "m", "height_above_anchor_isobar_m"),
        "refractivity": (level, target_run["refractivity"],
                         fp._attrs("1", "refractivity transferred to the target along the isobars "
                                   "(Eq. A28)", "modeled")),
        "refractivity_uncertainty": companion("1", "refractivity"),
        "refractivity_gauge": (level, np.exp(estimate.C)[order],
                               fp._attrs("1", "exp C, the anchor constant at the gauge latitude "
                                         "(Eq. A30)", "modeled")),
        "refractivity_gauge_uncertainty": companion("1", "refractivity_gauge"),
        "geopotential_gauge_uncertainty_m2s2": companion("m2 s-2", "geopotential_gauge_m2s2"),
        "geopotential_gauge_m2s2": (level, estimate.geopotential_m2s2[order],
                                    fp._attrs("m2 s-2", "the arrival geopotential at the gauge "
                                              "latitude, where the estimate was made", "modeled",
                                              positive="up")),
        "pressure_label_Pa": (level, target_run["label_pressure_Pa"],
                              fp._attrs("Pa", "the label of each isobar, carried from the gauge "
                                        "(decision H)", "index")),
        "pressure_identity_residual": (
            level, target_run["pressure_identity_residual"],
            fp._attrs("1", "p_produced / p_label - 1, the internal check of the whole chain "
                      "(decision H); derived from two variables of this file and so carrying no "
                      "companion of its own", "derived")),
        "u_column_ms": (level, target_run["u_column_ms"],
                        fp._attrs("m s-1", "the run's zonal wind along the target's column",
                                  "modeled")),
        "u_column_uncertainty_ms": companion("m s-1", "u_column_ms"),
        "mean_refractivity_m3": (level, target_run["mean_refractivity_m3"],
                                 fp._attrs("m3", "mean refractivity per molecule from the run's "
                                           "kind C at the target, on the isobar labels",
                                           "modeled")),
        "mean_refractivity_uncertainty_m3": companion("m3", "mean_refractivity_m3"),
        "mean_molar_mass_kg_mol": (level, target_run["mean_molar_mass_kg_mol"],
                                   fp._attrs("kg mol-1", "mean molar mass from the run's kind C at "
                                             "the target, on the isobar labels", "modeled")),
        "mean_molar_mass_uncertainty_kg_mol": companion("kg mol-1", "mean_molar_mass_kg_mol"),
        "number_density_m3": (level, target_run["number_density_m3"],
                              fp._attrs("m-3", "number density, N / R_bar (Eq. B4)", "modeled")),
        "number_density_uncertainty_m3": companion("m-3", "number_density_m3"),
        "pressure_Pa": (level, target_run["pressure_Pa"],
                        fp._attrs("Pa", "hydrostatic pressure, summed from the boundary at the top "
                                  "(Eq. B5)", "modeled",
                                  boundary_pressure_Pa=target_run["boundary_pressure_Pa"],
                                  positive="down", direction="increasing")),
        "pressure_uncertainty_Pa": companion("Pa", "pressure_Pa"),
        "temperature_K": (level, target_run["temperature_K"],
                          fp._attrs("K", "temperature, p R_bar / (k_B N) (Eq. B6)", "modeled")),
        "temperature_uncertainty_K": companion("K", "temperature_K"),
        "geopotential_uncertainty_m2s2": companion("m2 s-2", "geopotential_m2s2"),
        "latitude_planetocentric_deg": ((), target_deg,
                                        fp._attrs("degrees_north", "the target latitude", "index")),
        "gauge_latitude_planetocentric_deg": (
            (), math.degrees(estimate.gauge_latitude_rad),
            fp._attrs("degrees_north", "the gauge latitude, the weighted centroid of the "
                      "construction anchors (decision K)", "modeled")),
        "psi_deg": ((), float(np.degrees(target_run["psi_rad"][gauge_index])),
                    fp._attrs("degrees", "tilt of the local vertical at the target, on the gauge "
                              "isobar", "modeled")),
        "gauge_isobar_Pa": ((), float(namelist.gauge_isobar_Pa),
                            fp._attrs("Pa", "the declared gauge isobar, where Phi = 0", "index")),
        "gauge_level_index": ((), np.int32(gauge_index),
                              fp._attrs("1", "the level of the gauge isobar", "index")),
        "boundary_pressure_Pa": ((), target_run["boundary_pressure_Pa"],
                                 fp._attrs("Pa", "the boundary pressure p_b, the label of the "
                                           "topmost isobar of the highest anchor (decision F)",
                                           "index")),
        "boundary_level_index": ((), np.int32(0),
                                 fp._attrs("1", "the level of the boundary", "index")),
        "datum_isobar_Pa": ((), float(namelist.datum_isobar_Pa),
                            fp._attrs("Pa", "the declared datum isobar, from which the delivered "
                                      "altitude is measured", "index")),
        "datum_geopotential_m2s2": ((), target_run["datum_geopotential_m2s2"],
                                    fp._attrs("m2 s-2", "the geopotential where the produced "
                                              "pressure equals the datum isobar", "modeled")),
        "datum_radius_m": ((), target_run["datum_radius_m"],
                           fp._attrs("m", "the target column's radius at the datum", "modeled")),
        "reference_surface_radius_m": (
            (), float(state.reference_radius_m[state.mesh.latitude_index(
                math.radians(target_deg))]),
            fp._attrs("m", "r0 at the target, the gauge isobar of the reference surface (Eq. B3)",
                      "modeled")),
        "gauge_latitude_planetocentric_uncertainty_deg": scalar_companion(
            "degrees_north", "gauge_latitude_planetocentric_deg"),
        "psi_uncertainty_deg": scalar_companion("degrees", "psi_deg"),
        "datum_geopotential_uncertainty_m2s2": scalar_companion("m2 s-2",
                                                                "datum_geopotential_m2s2"),
        "datum_radius_uncertainty_m": scalar_companion("m", "datum_radius_m"),
        "reference_surface_radius_uncertainty_m": scalar_companion(
            "m", "reference_surface_radius_m"),
    }
    dataset = xr.Dataset(variables, coords={
        "geopotential_m2s2": (level, target_run["geopotential_m2s2"],
                              fp._attrs("m2 s-2", "geopotential along the local vertical, zero on "
                                        "the gauge isobar", "modeled", positive="up",
                                        gauge_isobar_Pa=float(namelist.gauge_isobar_Pa)))})
    paths = ([entry.path for entry in namelist.anchors]
             + [namelist.inputs[k] for k in ctl.RUN_INPUT_KINDS] + [namelist.path])
    dataset.attrs = {
        "title": f"{namelist.name} transferred profile, kind profile",
        "profile_or_run": namelist.name,
        "role": "forward",
        "source": (f"{fp.TOOL} transfer of "
                   f"{', '.join(entry.slug for entry in namelist.anchors)} to {target_deg} "
                   f"degrees under the run {namelist.name}; {namelist.description}"),
        "epoch": namelist.date if namelist.date else PROFILE_EPOCH_WITHOUT_DATE,
        "solar_longitude_deg": float(namelist.solar_longitude_deg),
        "solar_longitude_source": "declared in the namelist",
        "anchor_solar_longitudes_deg": np.array(
            [float(a.tree.attrs["solar_longitude_deg"]) for a in anchors]),
        "boundary_pressure_Pa": target_run["boundary_pressure_Pa"],
        "latitude_planetocentric_absent_meaning": "point",
        "mode": namelist.mode,
        "input_hashes": cio.input_hashes(paths, product),
    }
    cio.history_append(
        dataset,
        f"{fp.TOOL}: {len(anchors)} anchor(s) transferred to {target_deg} degrees on "
        f"{order.size} arrival levels; isobars traced by RK4 on the working mesh, ln N by the "
        f"trapezoid along them (Eq. A24, A28); the estimate at "
        f"{math.degrees(estimate.gauge_latitude_rad):.6f} degrees (Eq. A30); pressure by the "
        f"hydrostatic integral from p_b = {target_run['boundary_pressure_Pa']!r} Pa; altitude "
        f"above the {float(namelist.datum_isobar_Pa)!r} Pa datum along the local vertical")
    return dataset


def _copied(nodes, vertical: str) -> dict:
    """A copied subtree's nodes with its own vertical dimension renamed away from `level`.

    A transfer product's `level` is the union of the anchors' arrival levels, which at M >= 2 is
    neither the anchor's level count nor an input's, and an `xarray` DataTree child cannot carry a
    dimension of the parent's name with another size: the M = 2 product could be written and not
    read back, which is how this was found. So a copied file's own vertical dimension is named for
    whose levels it is, at M = 1 as well, and the product's structure does not change with the
    number of anchors. The variables, their values and their attributes are the copied file's
    unchanged, which is what "verbatim" means here.
    """
    renamed = {}
    for name, node in nodes.items():
        if "level" in node.dims and "level" not in node.variables:
            node = node.rename_dims({"level": vertical})
        renamed[name] = node
    return renamed


def _transfer_groups(namelist, inputs, anchors, state, estimate, arrived, target_run,
                     curves_to_target):
    """The groups a transfer product carries beyond the inputs and the namelist.

    `anchors/<slug>` is the kind N file verbatim, as closure mode writes it, with this run's own
    per-anchor arrays in a subgroup of it rather than mixed into the copy: what the anchor carries
    is the anchor's, and what the transfer made of it is the transfer's. `reference_surface` and
    `isobars` are on the mesh's latitude nodes, `estimate` on the union levels.
    """
    groups = {}
    for index, anchor in enumerate(anchors):
        slug = anchor.slug
        groups.update(_copied(fp._nodes(anchor.tree, f"anchors/{slug}"),
                              "anchor_level"))
        present = estimate.present[index]
        groups[f"anchors/{slug}/transfer"] = xr.Dataset(
            {
                "C_i": (("union_level",), np.where(present, estimate.C_i[index], np.nan),
                        fp._attrs("1", "this anchor's estimate of the anchor constant at the gauge "
                                  "latitude, ln N transferred there (Eq. A29); absent where the "
                                  "anchor does not span the level", "modeled")),
                "sigma_ln_N_measurement": (
                    ("anchor_level",), np.asarray(anchor.sigma_ln_N_measurement, dtype="float64"),
                    fp._attrs("1", "the anchor's own uncertainty of ln N as reduced, scaled by its "
                              "measurement_uncertainty_scale (decision J)", "derived")),
                "sigma_ln_N_season": (
                    ("anchor_level",), np.asarray(anchor.sigma_ln_N_season, dtype="float64"),
                    fp._attrs("1", "the seasonal term of the anchor's uncertainty (decision J)",
                              "derived", season_term=anchor.season_term)),
                "sigma_ln_N_on_union": (("union_level",),
                                        np.where(present, estimate.sigma_ln_N[index], np.nan),
                                        fp._attrs("1", "the anchor's uncertainty on the union "
                                                  "levels, as the estimate weighted it", "derived")),
                "weight_in_estimate": (("union_level",),
                                       np.where(present, estimate.estimate_weight[index], np.nan),
                                       fp._attrs("1", "1 / (sigma_i^2 + P_i), the weight of Eq. A30",
                                                 "modeled")),
                "present_on_union": (("union_level",), present.astype("int8"),
                                     fp._attrs("1", "whether the anchor spans each union level",
                                               "index")),
                "geopotential_arrival_m2s2": (
                    ("anchor_level",), arrived[index].curves.geopotential_m2s2[
                        :, state.mesh.latitude_index(estimate.gauge_latitude_rad)],
                    fp._attrs("m2 s-2", "where this anchor's isobars arrive at the gauge latitude",
                              "modeled", positive="up")),
            },
            attrs={
                "role_weight": float(anchor.weight),
                "role": "construction" if anchor.is_construction else "validation",
                "measurement_uncertainty_scale": float(anchor.measurement_uncertainty_scale),
                "latitude_planetocentric_deg": float(anchor.latitude_planetocentric_deg),
                "transfer_variance_P_i": float(estimate.transfer_variance[index]),
                "reference_surface_residual_m": float(
                    state.reference_radius_m[state.mesh.latitude_index(
                        math.radians(anchor.latitude_planetocentric_deg))]
                    - np.asarray(anchor.radius_m, dtype="float64")[anchor.gauge_level_index]),
                "anchor_solar_longitude_deg": ("uniform" if anchor.anchor_season_deg is None
                                               else float(anchor.anchor_season_deg)),
                "run_solar_longitude_deg": float(anchor.run_season_deg),
                "season_matches_run": str(bool(anchor.season_matches_run)),
                **{f"propagation_{k}": str(v) for k, v in dict(anchor.propagation).items()},
            })

    latitude = np.degrees(state.mesh.latitude_rad)
    groups["reference_surface"] = xr.Dataset(
        {"reference_surface_radius_m": (("latitude_planetocentric",), state.reference_radius_m,
                                        fp._attrs("m", "r0 on the mesh's latitude nodes, the gauge "
                                                  "isobar marched by Eq. B3 from the first anchor",
                                                  "modeled")),
         "latitude_planetocentric_deg": (("latitude_planetocentric",), latitude,
                                         fp._attrs("degrees_north", "the mesh's latitude nodes",
                                                   "index"))})

    isobars = {"latitude_planetocentric_deg": (("latitude_planetocentric",), latitude,
                                               fp._attrs("degrees_north", "the mesh's latitude "
                                                         "nodes", "index"))}
    for index, anchor in enumerate(anchors):
        # Each anchor's curves are on its own levels, which need not be one another's and are not
        # the union the root carries: one dimension per anchor, named for it.
        vertical = f"level_{anchor.slug}"
        isobars[f"geopotential_{anchor.slug}_m2s2"] = (
            (vertical, "latitude_planetocentric"), arrived[index].curves.geopotential_m2s2,
            fp._attrs("m2 s-2", f"Phi_k(phi) for {anchor.slug}, the traced isobars (Eq. A24)",
                      "modeled", positive="up"))
        isobars[f"ln_refractivity_{anchor.slug}"] = (
            (vertical, "latitude_planetocentric"), arrived[index].ln_N_along,
            fp._attrs("1", f"ln N_k(phi) for {anchor.slug} along its isobars (Eq. A28)", "modeled"))
    isobars["geopotential_gauge_to_target_m2s2"] = (
        ("union_level", "latitude_gauge_to_target"), curves_to_target.geopotential_m2s2,
        fp._attrs("m2 s-2", "the curves carrying C from the gauge latitude to the target",
                  "modeled", positive="up"))
    isobars["latitude_gauge_to_target_deg"] = (
        ("latitude_gauge_to_target",), np.degrees(curves_to_target.latitude_rad),
        fp._attrs("degrees_north", "the nodes those curves cover", "index"))
    groups["isobars"] = xr.Dataset(isobars)

    differences = {}
    for (left, right), value in estimate.differences.items():
        differences[f"D_{left}_{right}"] = (
            ("union_level",), value,
            fp._attrs("1", f"C_{left} - C_{right} on the levels both cover (Eq. A33)", "modeled"))
    groups["estimate"] = xr.Dataset(
        {
            "geopotential_m2s2": (("union_level",), estimate.geopotential_m2s2,
                                  fp._attrs("m2 s-2", "the union of the construction anchors' "
                                            "arrival levels at the gauge latitude", "modeled",
                                            positive="up")),
            "C": (("union_level",), estimate.C,
                  fp._attrs("1", "the anchor constant, the inverse variance mean of Eq. A30",
                            "modeled")),
            "C_variance": (("union_level",), estimate.variance,
                           fp._attrs("1", "1 / sum w_i, a diagnostic and not the product's "
                                     "uncertainty", "modeled")),
            "C_variance_inflated": (("union_level",), estimate.variance_inflated,
                                    fp._attrs("1", "the variance inflated by the reduced "
                                              "chi-square where it exceeds one", "modeled")),
            "chi_square_reduced": (("union_level",), estimate.chi_square_reduced,
                                   fp._attrs("1", "the reduced chi-square of the C_i about C, "
                                             "absent where fewer than two anchors are present",
                                             "modeled")),
            "label_pressure_Pa": (("union_level",), estimate.label_pressure_Pa,
                                  fp._attrs("Pa", "the isobar label of each union level, the "
                                            "weighted mean of the anchors' labels", "modeled")),
            "label_spread_ln_p": (("union_level",), estimate.label_spread_ln_p,
                                  fp._attrs("1", "the spread of the anchors' labels on each level, "
                                            "in ln p (decision H at M >= 2)", "modeled")),
            **differences,
        },
        attrs={
            "gauge_latitude_planetocentric_deg": math.degrees(estimate.gauge_latitude_rad),
            "identity_floor": IDENTITY_FLOOR,
            "identity_floor_note": IDENTITY_FLOOR_NOTE,
            **{k: str(v) for k, v in dict(estimate.record).items()},
        })
    return groups


def _transfer_record(namelist, inputs, anchors, state, estimate, target_run) -> dict:
    """What the run did, for the product's `transfer_record`. SPEC_04 Step 5 deliverable 3."""
    residual = target_run["pressure_identity_residual"]
    shift = {}
    for index, anchor in enumerate(anchors):
        curves = state.curves[index]
        at_anchor = state.mesh.latitude_index(curves.start_latitude_rad)
        at_target = state.mesh.latitude_index(math.radians(float(namelist.target_latitude_deg)))
        moved = (curves.geopotential_m2s2[:, at_target]
                 - curves.geopotential_m2s2[:, at_anchor])
        worst = int(np.argmax(np.abs(moved)))
        shift[f"isobar_shift_largest_{anchor.slug}_m2s2"] = float(moved[worst])
        shift[f"isobar_shift_largest_{anchor.slug}_level"] = worst
    return {
        "mode": namelist.mode,
        "run": namelist.name,
        "anchors": ", ".join(a.slug for a in anchors),
        "target_latitude_planetocentric_deg": float(namelist.target_latitude_deg),
        "gauge_latitude_planetocentric_deg": math.degrees(estimate.gauge_latitude_rad),
        "gauge_isobar_Pa": float(namelist.gauge_isobar_Pa),
        "datum_isobar_Pa": float(namelist.datum_isobar_Pa),
        "boundary_pressure_Pa": target_run["boundary_pressure_Pa"],
        "mesh_latitude_spacing_deg": float(namelist.grid["latitude_spacing_deg"]),
        "mesh_geopotential_spacing_m2s2": float(namelist.grid["geopotential_spacing_m2s2"]),
        "mesh_latitude_nodes": int(state.mesh.shape[0]),
        "mesh_geopotential_nodes": int(state.mesh.shape[1]),
        "mesh_latitude_range_deg": (f"{np.degrees(state.mesh.latitude_rad[0]):.6f} to "
                                    f"{np.degrees(state.mesh.latitude_rad[-1]):.6f}"),
        "mesh_geopotential_range_m2s2": (f"{state.mesh.geopotential_m2s2[0]:.1f} to "
                                         f"{state.mesh.geopotential_m2s2[-1]:.1f}"),
        "outer_loop_passes": int(state.record["passes"]),
        "outer_loop_residuals_ln_p": ", ".join(
            f"{h['residual_ln_p']:.6e}" for h in state.record["history"]),
        "outer_loop_relative_tolerance_ln_p": float(state.record["relative_tolerance_ln_p"]),
        "mesh_extensions": (", ".join(
            f"{e['side']} by {e['added_m2s2']:.1f} on pass {e['pass']}"
            for e in state.record["extensions"]) or "none"),
        "shear_kernel_largest_abs_s_over_g_per_rad": float(np.max(np.abs(state.s_over_g))),
        "pressure_identity_largest_abs": float(np.max(np.abs(residual))),
        "pressure_identity_largest_level": int(np.argmax(np.abs(residual))),
        "pressure_identity_rms": float(np.sqrt(np.mean(residual ** 2))),
        "gauge_latitude_rule": str(namelist.estimation["gauge_latitude_rule"]),
        "kernel_uncertainty_per_rad": float(namelist.estimation["kernel_uncertainty_per_rad"]),
        "model_error_correlation_length_deg": float(
            namelist.estimation["model_error_correlation_length_deg"]),
        "identity_floor": IDENTITY_FLOOR,
        "wind_solar_longitude_deg": str(inputs.seasons.get("wind", "not declared")),
        "composition_solar_longitude_deg": str(inputs.seasons.get("composition", "not declared")),
        "run_solar_longitude_deg": float(namelist.solar_longitude_deg),
        "latitude_drift_neglected": (
            "the drift of each anchor along its own vertical is neglected in latitude and "
            "recorded: 0.027 degrees at the top of the Lindal profile (SPEC_04 section 1, "
            "decision G)"),
        **shift,
        **{k: str(v) for k, v in dict(state.record).items() if k not in ("history", "extensions")},
    }


@dataclass(frozen=True)
class Carried:
    """What a run leaves between the outer loop and the product."""

    state: TransferState
    estimate: object
    arrived: tuple
    curves_to_target: Curves
    target: MappingProxyType


def chain(inputs, anchors, *, gauge_latitude_rad, target_latitude_rad, gauge_isobar_Pa, p_b,
          latitude_spacing_rad, geopotential_spacing_m2s2, relative_tolerance_ln_p,
          max_iterations, kernel_uncertainty_per_rad, datum_isobar_Pa, field=None,
          mesh=None) -> Carried:
    """The outer loop, every anchor to the gauge, the estimate there, `C` to the target, and the
    production at the target. SPEC_04 Step 5 deliverable 4 between the inputs and the product.

    `run` is this with the namelist read, the product written and the figures rendered. It is a
    function of its own rather than part of `run` so that a run built on a mesh given to the outer
    loop can be measured, which the M = 2 identity needs (SPEC_04 v0.15 Step 4 ruling 1) and which
    `run`, taking everything from the namelist, has no argument for.
    """
    anchors = tuple(anchors)
    field = wf.WindField(inputs.wind) if field is None else field
    state = outer_loop(
        inputs, list(anchors), gauge_latitude_rad=gauge_latitude_rad,
        target_latitude_rad=target_latitude_rad, gauge_isobar_Pa=gauge_isobar_Pa, p_b=p_b,
        latitude_spacing_rad=latitude_spacing_rad,
        geopotential_spacing_m2s2=geopotential_spacing_m2s2,
        relative_tolerance_ln_p=relative_tolerance_ln_p, max_iterations=max_iterations,
        field=field, mesh=mesh)

    arrived, arrivals = [], []
    at_gauge_node = state.mesh.latitude_index(gauge_latitude_rad)
    for index, anchor in enumerate(anchors):
        placement = state.placements[index]
        curves, at_gauge, along = carry_to(
            state, inputs, placement.latitude_rad, gauge_latitude_rad,
            placement.geopotential_m2s2, anchor.ln_N, placement.label_pressure_Pa)
        arrivals.append(fe.Arrival(
            slug=anchor.slug, latitude_rad=placement.latitude_rad, weight=anchor.weight,
            geopotential_m2s2=state.curves[index].geopotential_m2s2[:, at_gauge_node],
            C_i=along[:, at_gauge], sigma_ln_N=fe.sigma_ln_N(anchor),
            label_pressure_Pa=placement.label_pressure_Pa))
        arrived.append(Arrived(anchor=anchor, curves=state.curves[index],
                               ln_N_along=transfer(
                                   state.mesh, state.s_over_g, state.curves[index], anchor.ln_N,
                                   *lk.composition_term(inputs.composition,
                                                        placement.label_pressure_Pa)),
                               arrival=arrivals[-1]))

    estimate = fe.estimate(arrivals, float(kernel_uncertainty_per_rad), gauge_latitude_rad)
    curves_to_target, at_target, carried = carry_to(
        state, inputs, gauge_latitude_rad, target_latitude_rad, estimate.geopotential_m2s2,
        estimate.C, estimate.label_pressure_Pa)
    target_run = produce_at_target(
        state, inputs, target_latitude_rad, curves_to_target.geopotential_m2s2[:, at_target],
        carried[:, at_target], estimate.label_pressure_Pa, p_b, field, float(datum_isobar_Pa))
    return Carried(state=state, estimate=estimate, arrived=tuple(arrived),
                   curves_to_target=curves_to_target, target=target_run)


def write_product(namelist, inputs, anchors, state, estimate, arrived, target_run,
                  curves_to_target):
    """Write kind `profile` in transfer mode and return the path. Step 5 deliverable 3.

    The root, the groups of deliverable 3, the four inputs, the namelist and the record, in the
    order kind `profile` carries them. `run` writes its product with this, and so does the M = 2
    acceptance run, which is built on a mesh given to the outer loop rather than one the loop
    builds (SPEC_04 v0.15 Step 4 ruling 1) and therefore cannot go through `run`: one path, so
    that what a suite writes is what a run writes.
    """
    product = namelist.product
    dataset = _transfer_product(namelist, inputs, anchors, state, estimate, target_run,
                                curves_to_target, product)
    groups = _transfer_groups(namelist, inputs, anchors, state, estimate, arrived, target_run,
                              curves_to_target)
    for key in ctl.RUN_INPUT_KINDS:
        loaded = getattr(inputs, key)
        if isinstance(loaded, xr.DataTree):
            groups.update(_copied(fp._nodes(loaded, f"inputs/{key}"), "input_level"))
        else:
            groups.update(_copied({f"inputs/{key}": loaded}, "input_level"))
    groups["namelist"] = xr.Dataset(attrs={
        "text": namelist.text, "sha256": namelist.sha256,
        "resolved": fp._resolved_namelist(namelist, inputs, product)})
    groups["transfer_record"] = xr.Dataset(attrs=_transfer_record(
        namelist, inputs, anchors, state, estimate, target_run))
    return cio.write(product, dataset, "profile", groups=groups, created_by=fp.TOOL)


def run(namelist_path) -> TransferResult:
    """The transfer run of SPEC_04 Step 5 deliverable 4.

    Load, propagate every anchor through the hook, fix the gauge latitude before the mesh
    (decision K), run the outer loop, carry every anchor to the gauge, estimate the anchor
    constant there, carry `C` to the target, produce, place the datum, write kind `profile`, and
    render the figures when the namelist asks. Every file is read and closed before the next is
    opened, which `lib.control` does for the inputs and this does for nothing else: the product is
    written once at the end (section 15 ruling 3).
    """
    namelist = ctl.read_run_namelist(namelist_path)
    if namelist.mode != "transfer":
        raise ctl.ControlFileError(
            f"{namelist.path}: mode {namelist.mode!r} is not this driver's; "
            "forward.production.run is the closure.")
    inputs = ctl.load_run_inputs(namelist)
    anchors = tuple(fprop.propagate(a, namelist.solar_longitude_deg) for a in inputs.anchors)
    gauge_latitude_rad = fe.gauge_latitude(anchors)
    target_latitude_rad = math.radians(float(namelist.target_latitude_deg))
    p_b = (boundary_pressure(anchors) if namelist.p_b_Pa is None else float(namelist.p_b_Pa))
    field = wf.WindField(inputs.wind)
    loop = namelist.numerics["outer_loop"]

    carried = chain(
        inputs, anchors, gauge_latitude_rad=gauge_latitude_rad,
        target_latitude_rad=target_latitude_rad, gauge_isobar_Pa=namelist.gauge_isobar_Pa,
        p_b=p_b, latitude_spacing_rad=math.radians(float(namelist.grid["latitude_spacing_deg"])),
        geopotential_spacing_m2s2=float(namelist.grid["geopotential_spacing_m2s2"]),
        relative_tolerance_ln_p=loop["relative_tolerance_ln_p"],
        max_iterations=loop["max_iterations"],
        kernel_uncertainty_per_rad=float(namelist.estimation["kernel_uncertainty_per_rad"]),
        datum_isobar_Pa=float(namelist.datum_isobar_Pa), field=field)
    state, estimate, target_run = carried.state, carried.estimate, carried.target

    written = write_product(namelist, inputs, anchors, state, estimate, carried.arrived,
                            target_run, carried.curves_to_target)

    figures = None
    diagnostics = dict(namelist.diagnostics)
    if diagnostics.get("figures"):
        # Drawn from the file just written, never from this run's memory (SPEC_02 Step 5).
        from casspian.tools.plots import render

        figures = render(written, namelist.output_directory / "figures", diagnostics["format"],
                         int(diagnostics["dpi"]))
    return TransferResult(product=written, state=state, estimate=estimate,
                          arrived=carried.arrived, target=target_run, figures=figures)
