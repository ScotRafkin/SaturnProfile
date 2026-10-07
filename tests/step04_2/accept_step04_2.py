"""Acceptance checks for SPEC_04 v0.10 Step 2: `lib.mesh`, the columns, `z_lv`, and the
cylinder-extended wind file of decision P.

Every check prints its measured value; checks beyond the specification are labeled so. The anchor
and the run's inputs are the clean swept products of the Step 1 acceptance commit `780a2de`, so
nothing is relaxed and no `-dirty` refusal is touched.

This step registers nothing. The cylinder-extended wind is written under `reports/step04_2/`,
which git ignores, and read back through `wind_at` like any other file; the closure run's own
inputs are not touched.

Checks 1a and 1b are the specification's uniform-gravity check exactly as it is written, and both
fail. They are reported as failures with the measured value and the reason, and the statements
that are exactly true are checked beside them (1c, 1d). Nothing was loosened.
"""

import math
import sys
from pathlib import Path

import numpy as np

from casspian.forward import production as fp
from casspian.lib import control as ctl
from casspian.lib import geoid
from casspian.lib import io as cio
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf
from casspian.lib.gravity import g_eff_vector
from casspian.lib.schema import CasspianSchemaError

HERE = Path("reports/step04_2")
CLOSURE = Path("forward/lindal_closure")
TOOL = "tests/step04_2/accept_step04_2.py"

HERE.mkdir(parents=True, exist_ok=True)
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def refusal(callable_, *args, **kwargs):
    try:
        callable_(*args, **kwargs)
    except Exception as exc:
        return f"{type(exc).__name__}: {exc}"
    return None


def worst_relative(got, want):
    got, want = np.asarray(got, dtype="float64"), np.asarray(want, dtype="float64")
    scale = np.where(want == 0.0, 1.0, np.abs(want))
    return float(np.max(np.abs(got - want) / scale))


# ---------------------------------------------------------------------------
# The run, the anchor, the closure geometry
# ---------------------------------------------------------------------------
namelist = ctl.read_run_namelist(CLOSURE / "lindal_closure.toml")
inputs = ctl.load_run_inputs(namelist)
anchor = inputs.anchor
anchor_root = anchor.to_dataset(inherit=False)

phi_a = math.radians(float(anchor_root["latitude_planetocentric_deg"].values))
r_a = float(anchor_root["anchor_isobar_radius_m"].values)
h_levels = np.asarray(anchor_root["height_above_anchor_isobar_m"].values, dtype="float64")
p_tab = np.asarray(anchor["inputs"]["thermo"]["pressure_Pa"].values, dtype="float64")
gauge_level = int(np.flatnonzero(p_tab == namelist.gauge_isobar_Pa)[0])

Omega = float(inputs.rotation["angular_rate_rad_s"])
GM = float(inputs.gravity["GM_m3s2"])
J = np.asarray(inputs.gravity["J"].values)
DEGREES = np.asarray(inputs.gravity["degree"].values)
R_NORM = float(inputs.gravity["normalization_radius_m"])
CONSTANTS = (Omega, GM, J, DEGREES, R_NORM)

closure_field = wf.WindField(inputs.wind)
lat_file = np.radians(np.asarray(inputs.wind["latitude_planetocentric_deg"].values,
                                 dtype="float64"))
p_file = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
p_ref = float(inputs.wind["reference_level_pressure_Pa"])
reference_column = int(np.flatnonzero(p_file == p_ref)[0])

profile = fp.profile_from_anchor(anchor)
p_b = float(p_tab[0])
closure_column = np.asarray(
    closure_field.wind_at(np.full(p_tab.shape, phi_a), p_tab), dtype="float64")
closure_production = fp.produce(profile, inputs, gauge_level, p_b, u_column=closure_column)
Phi_closure = closure_production.geopotential.geopotential_m2s2

LATITUDE_SPACING = math.radians(0.05)
GEOPOTENTIAL_SPACING = 5.0e3
target = math.radians(10.0)
mesh = lm.build_mesh([phi_a, target], Phi_closure, LATITUDE_SPACING, GEOPOTENTIAL_SPACING)


def surface_on(latitudes, field):
    """The Step 1 reference surface through the anchor, under a field's gauge-isobar wind."""
    return geoid.through_anchor(
        latitudes, phi_a, r_a,
        lambda x: field.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
        *CONSTANTS).radius


def isobar_maps(Phi_k):
    """`Phi(p)` linear in `ln p` with linear extrapolation, and its inverse `p(Phi)`.

    The flat-isobar map of the anchor, which decision P names. The file's pressure grid reaches
    1 Pa and 1e6 Pa while the anchor's levels span 19.95 to 129848 Pa, so 23 of the 61 columns
    fall outside it; `np.interp` would clamp them all onto the two end values, which is why the
    ends are continued at the slope of the last interval instead.
    """
    order = np.argsort(np.log(p_tab))
    x, y = np.log(p_tab)[order], Phi_k[order]
    low = (y[1] - y[0]) / (x[1] - x[0])
    high = (y[-1] - y[-2]) / (x[-1] - x[-2])

    def forward(p):
        lp = np.log(np.asarray(p, dtype="float64"))
        return np.where(lp < x[0], y[0] + low * (lp - x[0]),
                        np.where(lp > x[-1], y[-1] + high * (lp - x[-1]), np.interp(lp, x, y)))

    rank = np.argsort(Phi_k)
    return forward, (lambda Phi: np.exp(np.interp(Phi, Phi_k[rank], np.log(p_tab)[rank])))


# ---------------------------------------------------------------------------
# 1. Uniform gravity, as the specification states it, and what is exactly true
# ---------------------------------------------------------------------------
J_zero = np.zeros_like(J)
flat_constants = (0.0, GM, J_zero, DEGREES, R_NORM)
Phi_nodes = mesh.geopotential_m2s2
flat = lm.column(phi_a, Phi_nodes, r_a, (lambda Phi: 0.0), *flat_constants)
g0 = GM / r_a ** 2
g_top = GM / float(flat.radius_m[-1]) ** 2
exact = r_a / (1.0 - r_a * Phi_nodes / GM)
central = worst_relative(flat.radius_m, exact)
flat_form = worst_relative(flat.radius_m - r_a, Phi_nodes / g0)
record(1, "the central field of J = 0, Omega = 0, u = 0 returns the closed form "
          "r = r0 / (1 - r0 Phi / GM) to 1e-12 relative (SPEC_04 section 13 ruling 1)",
       central <= 1e-12,
       f"max relative departure {central:.3e}, bound 1e-12\n"
       f"dr/dPhi = r^2 / GM integrates to 1/r0 - 1/r = Phi / GM, which is the closed form above\n"
       f"the v0.10 form of this check, r - r0 = Phi / g, departs by {flat_form:.3e}: those "
       f"constants leave GM / r^2, not a uniform g, and it runs {g0:.6f} m/s2 at the gauge to "
       f"{g_top:.6f} at the top node, a change of {100 * abs(g_top - g0) / g0:.2f} percent")

g_flat, G_phi_flat, magnitude_flat, _ = g_eff_vector(0.0, r_a, phi_a, *flat_constants)
difference = flat.z_local_vertical_m - (flat.radius_m - r_a)
apart = np.abs(flat.radius_m - r_a) > 1.0
relative_z = float(np.max(np.abs(difference[apart] / (flat.radius_m - r_a)[apart])))
absolute_z = float(np.max(np.abs(difference)))
record("1b", "in that central field z_lv = r - r0 to 1e-7 m absolute (SPEC_04 section 13 ruling 2)",
       absolute_z <= 1e-7,
       f"largest absolute departure {absolute_z:.3e} m, bound 1e-7; in relative terms "
       f"{relative_z:.3e} over the {int(apart.sum())} nodes where |r - r0| > 1 m\n"
       f"G_phi is {float(G_phi_flat)!r} and |g_eff| - g is {float(magnitude_flat) - float(g_flat)!r} "
       f"exactly, so the two integrations are the same arithmetic and the difference is not "
       f"physical: r accumulates from {r_a:.0f} m and z from 0, so r - r0 cancels about seven "
       f"digits. One ulp of r is {np.spacing(r_a):.2e} m, and the departure is "
       f"{absolute_z / np.spacing(r_a):.0f} ulp, which is the bound the ruling sets")

# ---------------------------------------------------------------------------
# 2. The column at phi_c against the stated values
# ---------------------------------------------------------------------------
rank = np.argsort(Phi_closure)
forward_map, inverse_map = isobar_maps(Phi_closure)
r0_anchor = float(surface_on(np.array([phi_a]), closure_field)[0])
anchor_on_levels = lm.column(
    phi_a, Phi_closure[rank], r0_anchor,
    (lambda Phi: closure_field.wind_at(phi_a, inverse_map(Phi))), *CONSTANTS)
back = np.argsort(rank)
z_lv = anchor_on_levels.z_local_vertical_m[back]
dr = (anchor_on_levels.radius_m - r0_anchor)[back]
wanted = {"z_lv top": (float(z_lv[0]), 286715.2), "z_lv bottom": (float(z_lv[-1]), -104098.0),
          "r - r0 top": (float(dr[0]), 288051.3), "r - r0 bottom": (float(dr[-1]), -104576.8)}
worst = max(abs(a - b) for a, b in wanted.values())
record(2, "the column at phi_c reproduces the stated z_lv and r - r0 to 0.1 m",
       worst <= 0.1,
       "\n".join(f"{name:14s} {got:12.1f} m, expected {want:10.1f}, difference {got - want:+.3f}"
                 for name, (got, want) in wanted.items())
       + f"\nlargest difference {worst:.3f} m\n"
       f"the tabulated h - h_ref is {h_levels[0]:.1f} and {h_levels[-1]:.1f} m; the 15 m and 2 m "
       f"are the second-order drift terms\n"
       f"(r - r0) / z_lv at the top is {float(dr[0] / z_lv[0]):.6f} and 1 / cos psi is "
       f"{1.0 / math.cos(float(anchor_on_levels.psi_rad[back][0])):.6f} at psi = "
       f"{math.degrees(float(anchor_on_levels.psi_rad[back][0])):.4f} degrees")

# ---------------------------------------------------------------------------
# 3. Phi recovered from the column's own r by the trapezoid of g
# ---------------------------------------------------------------------------
mesh_column = lm.column(
    phi_a, Phi_nodes, r0_anchor,
    (lambda Phi: closure_field.wind_at(phi_a, inverse_map(Phi))), *CONSTANTS)
g_nodes, r_nodes = mesh_column.g_radial_ms2, mesh_column.radius_m
origin = mesh.gauge_geopotential_index
increments = 0.5 * (g_nodes[:-1] + g_nodes[1:]) * np.diff(r_nodes)
recovered = np.empty_like(Phi_nodes)
recovered[origin] = 0.0
recovered[origin + 1:] = np.cumsum(increments[origin:])
recovered[:origin] = -np.cumsum(increments[:origin][::-1])[::-1]
nonzero = Phi_nodes != 0.0
worst_phi = float(np.max(np.abs((recovered[nonzero] - Phi_nodes[nonzero]) / Phi_nodes[nonzero])))
record(3, "Phi recomputed from the column's own r by the trapezoid of g returns the nodes to "
          "1e-9 relative",
       worst_phi <= 1e-9,
       f"max relative departure {worst_phi:.3e} over {int(nonzero.sum())} nodes, bound 1e-9\n"
       f"{Phi_nodes.size} nodes at {GEOPOTENTIAL_SPACING:.0f} m2/s2, largest layer "
       f"{float(np.max(np.abs(np.diff(r_nodes)))):.0f} m; the gauge node is excluded because "
       f"a relative test there divides by zero")

# ---------------------------------------------------------------------------
# 4. The mesh refuses a spacing that is not positive
# ---------------------------------------------------------------------------
cases = [("latitude spacing zero", (0.0, GEOPOTENTIAL_SPACING)),
         ("latitude spacing negative", (-LATITUDE_SPACING, GEOPOTENTIAL_SPACING)),
         ("geopotential spacing zero", (LATITUDE_SPACING, 0.0)),
         ("geopotential spacing negative", (LATITUDE_SPACING, -GEOPOTENTIAL_SPACING)),
         ("latitude spacing not finite", (float("nan"), GEOPOTENTIAL_SPACING))]
lines, refused = [], 0
for name, (a, b) in cases:
    message = refusal(lm.build_mesh, [phi_a, target], Phi_closure, a, b)
    refused += message is not None
    lines.append(f"{name}: {'refused, ' + message[:90] if message else 'NOT refused'}")
record(4, "the mesh refuses a spacing that is not positive", refused == len(cases),
       "\n".join(lines))

# ---------------------------------------------------------------------------
# 5. extend adds whole spacings, and the columns match a mesh built that way
# ---------------------------------------------------------------------------
def columns_of(a_mesh):
    r0 = surface_on(a_mesh.latitude_rad, closure_field)
    return lm.build_columns(
        a_mesh, r0,
        (lambda i: (lambda Phi, la=float(a_mesh.latitude_rad[i]):
                    closure_field.wind_at(la, inverse_map(Phi)))), *CONSTANTS)


lines, ok = [], True
for side, amount in (("above", 12000.0), ("below", 7000.0)):
    extended = mesh.extend(side, amount)
    steps = int(np.ceil(amount / GEOPOTENTIAL_SPACING))
    nodes = extended.geopotential_m2s2
    grew = nodes.size - Phi_nodes.size == steps
    whole = bool(np.allclose(np.diff(nodes), GEOPOTENTIAL_SPACING, rtol=0, atol=1e-9))
    kept = (np.array_equal(nodes[steps:], Phi_nodes) if side == "below"
            else np.array_equal(nodes[:Phi_nodes.size], Phi_nodes))
    ok = ok and grew and whole and kept
    lines.append(f"{side} by {amount:.0f}: {steps} whole spacings added, node count "
                 f"{Phi_nodes.size} -> {nodes.size}, existing nodes kept exactly {kept}, "
                 f"run still uniform {whole}")

wide = mesh.extend("above", 12000.0)
# "A mesh built with them from the start" means one whose own build lands on the extended nodes.
# `build_mesh` always adds a spacing of margin each side, so it is given anchor levels one node
# inside the extended range; passing the extended nodes themselves would add the margin twice.
inside = np.array([wide.geopotential_m2s2[1], wide.geopotential_m2s2[-2]])
from_scratch = lm.build_mesh([phi_a, target], inside, LATITUDE_SPACING, GEOPOTENTIAL_SPACING)
same_nodes = np.array_equal(wide.geopotential_m2s2, from_scratch.geopotential_m2s2)
a_columns, b_columns = columns_of(wide), columns_of(from_scratch)
column_gap = max(worst_relative(a.radius_m, b.radius_m) for a, b in zip(a_columns, b_columns))
z_gap = max(worst_relative(a.z_local_vertical_m, b.z_local_vertical_m)
            for a, b in zip(a_columns, b_columns))
ok = ok and same_nodes and column_gap <= 1e-12 and z_gap <= 1e-12
lines.append(f"a mesh extended above by 12000 has the same nodes as one built with them from the "
             f"start: {same_nodes}; the columns agree to {column_gap:.3e} in r and {z_gap:.3e} in "
             f"z_lv over {len(a_columns)} latitudes, bound 1e-12")
record(5, "extend adds whole spacings on the named side and the columns on the new nodes match a "
          "mesh built with them from the start to 1e-12", ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# 6. Two anchors at different latitudes
# ---------------------------------------------------------------------------
second = math.radians(60.0)
two = lm.build_mesh([phi_a, second, target], Phi_closure, LATITUDE_SPACING, GEOPOTENTIAL_SPACING)
exact_nodes = {name: float(two.latitude_rad[two.latitude_index(value)]) == value
               for name, value in (("the Lindal anchor", phi_a), ("a second anchor at 60 N", second),
                                   ("the target at 10 N", target))}
span_ok = (two.latitude_rad[0] <= target - LATITUDE_SPACING + 1e-15
           and two.latitude_rad[-1] >= second + LATITUDE_SPACING - 1e-15)
record(6, "a namelist with two anchors at different latitudes gives a mesh with both as exact nodes",
       all(exact_nodes.values()) and span_ok,
       "\n".join(f"{name}: an exact node {value}" for name, value in exact_nodes.items())
       + f"\n{two.latitude_rad.size} latitude nodes from {math.degrees(two.latitude_rad[0]):.3f} to "
       f"{math.degrees(two.latitude_rad[-1]):.3f} degrees, one spacing of margin each side: "
       f"{span_ok}\n"
       f"the three inserted latitudes split their cells, so the run is uniform except there: "
       f"{int(np.sum(~np.isclose(np.diff(two.latitude_rad), LATITUDE_SPACING, rtol=1e-9)))} "
       f"cells differ from {math.degrees(LATITUDE_SPACING):.3f} degrees")

# ---------------------------------------------------------------------------
# 7. Every quantity on the mesh is finite
# ---------------------------------------------------------------------------
all_columns = columns_of(mesh)
fields = ("radius_m", "z_local_vertical_m", "g_radial_ms2", "G_phi_ms2", "g_magnitude_ms2",
          "psi_rad", "u_ms")
counts = {name: int(sum(int(np.sum(~np.isfinite(getattr(c, name)))) for c in all_columns))
          for name in fields}
record(7, "every quantity on the mesh is finite",
       all(v == 0 for v in counts.values()),
       f"{len(all_columns)} columns by {Phi_nodes.size} nodes = "
       f"{len(all_columns) * Phi_nodes.size} points\n"
       + "; ".join(f"{name} {value} non-finite" for name, value in counts.items()))

# ---------------------------------------------------------------------------
# 8 and 9. The cylinder-extended wind of decision P
# ---------------------------------------------------------------------------
# The inversion curve s_ref(phi) is sampled at the mesh's latitude spacing over the file's whole
# range, so that no latitude a level maps to is more than one mesh spacing from a sample (SPEC_04
# section 13 ruling 3). Sampling it on the file's 0.5 degree grid instead makes the inversion use
# that cell's chord, which cost up to 1e-2 m/s on the anchor's column and is what the first filing
# measured as a uniform offset. The file's own nodes and the run's exact latitudes are kept in the
# sample so that the reference-level identity stays exact on them.
fine = np.arange(float(lat_file[0]), float(lat_file[-1]) + 0.5 * LATITUDE_SPACING,
                 LATITUDE_SPACING)
lat_inversion = np.unique(np.clip(np.concatenate([fine, lat_file, [phi_a, target]]),
                                  float(lat_file[0]), float(lat_file[-1])))
on_file = np.searchsorted(lat_inversion, lat_file)
northern = lat_inversion >= 0.0

Phi_k, field_now, passes, history = Phi_closure.copy(), closure_field, 0, []
while True:
    forward_map, inverse_map = isobar_maps(Phi_k)
    Phi_of_pressure = forward_map(p_file)
    nodes = np.unique(np.concatenate([Phi_of_pressure, Phi_k, [0.0]]))
    at_pressure = np.searchsorted(nodes, Phi_of_pressure)
    at_level = np.searchsorted(nodes, Phi_k)
    r0 = surface_on(lat_inversion, field_now)
    radius = np.empty((lat_inversion.size, p_file.size))
    anchor_column = None
    for i, latitude in enumerate(lat_inversion):
        built = lm.column(latitude, nodes, float(r0[i]),
                          (lambda Phi, la=latitude: field_now.wind_at(la, inverse_map(Phi))),
                          *CONSTANTS)
        radius[i] = built.radius_m[at_pressure]
        if latitude == phi_a:
            anchor_column = built
    s_reference = radius[:, reference_column] * np.cos(lat_inversion)

    # The inversion is per hemisphere. A cylinder of radius s cuts the reference surface once in
    # each hemisphere, so one global U(s) would force u_north = u_south at equal s; the wind is
    # not symmetric (107.8 m/s at 60 N against 27.3 at 60 S), and a global inversion breaks the
    # pinning identity in the south by up to 112 m/s, measured. Each hemisphere carries its own
    # U(s), which keeps du/dz zero within a hemisphere, which is what S = 0 needs.
    def branch(mask):
        s, p = s_reference[mask], lat_inversion[mask]
        rising = np.argsort(s)
        return s[rising], p[rising]

    s_north, phi_north = branch(lat_inversion >= 0.0)
    s_south, phi_south = branch(lat_inversion <= 0.0)
    s_grid = radius[on_file] * np.cos(lat_file)[:, None]
    in_north = lat_file >= 0.0
    phi_star = np.empty_like(s_grid)
    phi_star[in_north] = np.interp(s_grid[in_north], s_north, phi_north)
    phi_star[~in_north] = np.interp(s_grid[~in_north], s_south, phi_south)
    u_total = np.interp(phi_star, lat_file, u_reference_file)
    u_total[0], u_total[-1] = 0.0, 0.0
    s_nodes, phi_nodes = s_north, phi_north
    change = float(np.max(np.abs(u_total - field_now.u_total_ms)))
    passes += 1
    history.append(change)
    # Assigning a tuple replaces the variable and drops its attributes with it, and kind W
    # requires units on every variable (SPEC_00 section 5), so the source file's attributes are
    # carried over and only the value_source is restated for what this field actually is.
    built = inputs.wind.copy(deep=True)
    dims = ("latitude_planetocentric", "pressure")
    total_attrs = dict(inputs.wind["u_total_ms"].attrs)
    total_attrs["value_source"] = ("the closure file's u_reference extended along cylinders "
                                   "(SPEC_04 decision P)")
    shear_attrs = dict(inputs.wind["u_shear_ms"].attrs)
    built["u_total_ms"] = (dims, u_total, total_attrs)
    built["u_shear_ms"] = (dims, u_total - u_reference_file[:, None], shear_attrs)
    field_now = wf.WindField(built)
    read_back = np.asarray(field_now.wind_at(np.full(p_tab.shape, phi_a), p_tab), dtype="float64")
    Phi_k = fp.produce(profile, inputs, gauge_level, p_b,
                       u_column=read_back).geopotential.geopotential_m2s2
    if change < 1.0e-6 or passes >= 10:
        break

built.attrs["title"] = "lindal cylinder-extended wind, kind W (SPEC_04 decision P)"
built.attrs["method"] = ("the closure file's reference-level wind extended along cylinders "
                         "coaxial with the rotation axis, pinned on the file's own reference "
                         "level, on the library's radial columns under the anchor's flat-isobar "
                         "map; a fixed point with the Step 1 reference surface")
built.attrs["vertical_structure"] = "constant on cylinders"
path = HERE / "lindal_cylinder_wind.nc"
written = cio.write(path, built, "wind", created_by=TOOL)
reread = wf.WindField(cio.read(written, "wind"))

sum_identity = float(np.max(np.abs(
    np.asarray(cio.read(written, "wind")["u_total_ms"].values, dtype="float64")
    - (u_reference_file[:, None] + np.asarray(cio.read(written, "wind")["u_shear_ms"].values,
                                              dtype="float64")))))
poles = (float(np.max(np.abs(u_total[0]))), float(np.max(np.abs(u_total[-1]))))
identity = float(np.max(np.abs(u_total[:, reference_column] - u_reference_file)))
record(8, "the cylinder-extended wind of decision P converges, holds the reference-level identity "
          "and the sum identity, and is written and read back through the schema",
       # The sum identity is bounded at round-off, not at bitwise equality: `u_shear` is formed
       # as `u_total - u_reference` in float64 and read back through netCDF, so `u_reference +
       # u_shear` returns `u_total` to an ulp and not always to the bit. The poles and the
       # reference-level identity are exact by construction and are held to that.
       change < 1.0e-6 and sum_identity <= 1e-12 and max(poles) == 0.0 and identity == 0.0,
       f"{passes} passes, changes " + ", ".join(f"{c:.2e}" for c in history) + " m/s\n"
       f"reference-level identity u_total(phi, p_ref) = u_reference(phi): largest departure "
       f"{identity:.2e} m/s over {lat_file.size} latitudes\n"
       f"sum identity u_total = u_reference + u_shear, computed here from the written file: "
       f"{sum_identity:.2e} m/s; both poles exactly zero: {max(poles) == 0.0}\n"
       f"written to {written} through lib.io.write under kind wind, which validates its schema; "
       f"tools.wind.build_wind is the published-curve path and does not apply here\n"
       f"the file carries the run's season, solar_longitude_deg "
       f"{float(cio.read(written, 'wind').attrs['solar_longitude_deg'])}")

construction = np.interp(
    np.interp(anchor_column.radius_m[at_level] * math.cos(phi_a), s_nodes, phi_nodes),
    lat_file, u_reference_file)
one_bar = int(np.argmin(np.abs(p_tab - 1.0e5)))
expected = {"top": (0, 9.490), "gauge": (gauge_level, 3.717), "1 bar": (one_bar, 2.168),
            "bottom": (p_tab.size - 1, 1.926)}
worst_u = max(abs(float(construction[i]) - w) for i, w in expected.values())
shift = fp.produce(profile, inputs, gauge_level, p_b,
                   u_column=construction).geopotential.geopotential_m2s2 - Phi_closure
worst_shift = max(abs(float(shift[0]) + 322.3), abs(float(shift[-1]) - 18.0))
record(9, "on the anchor's column the cylinder wind reproduces the stated values to 1e-2 m/s and "
          "the geopotential shift to 1 m2/s2",
       worst_u <= 1e-2 and worst_shift <= 1.0,
       "\n".join(f"{name:6s} {float(construction[i]):10.6f} m/s, expected {w:6.3f}, "
                 f"difference {float(construction[i]) - w:+.6f}"
                 for name, (i, w) in expected.items())
       + f"\nlargest departure {worst_u:.6f} m/s, bound 1e-2\n"
       f"Phi shift at the top {float(shift[0]):+.3f} m2/s2, expected -322.3; at the bottom "
       f"{float(shift[-1]):+.3f}, expected +18.0; largest departure {worst_shift:.3f}, bound 1\n"
       f"the level nearest 1 bar is {p_tab[one_bar]:.1f} Pa and u_reference at the anchor is "
       f"{float(np.interp(phi_a, lat_file, u_reference_file)):.8f} m/s, which the construction "
       f"returns there to {abs(float(construction[one_bar]) - float(np.interp(phi_a, lat_file, u_reference_file))):.2e}")

# ---------------------------------------------------------------------------
# 10. The column at 10 N
# ---------------------------------------------------------------------------
at_target = np.asarray(reread.wind_at(np.full(p_tab.shape, target), p_tab), dtype="float64")
target_expected = {"top": (0, 427.09), "gauge": (gauge_level, 345.98),
                   "bottom": (p_tab.size - 1, 326.40)}
worst_target = max(abs(float(at_target[i]) - w) for i, w in target_expected.values())
on_reference = float(reread.wind_at(target, p_ref))
record(10, "at 10 N the column carries the stated values to 0.05 m/s",
       worst_target <= 0.05 and abs(on_reference - 329.46) <= 0.05,
       "\n".join(f"{name:6s} {float(at_target[i]):9.2f} m/s, expected {w:7.2f}, difference "
                 f"{float(at_target[i]) - w:+.3f}" for name, (i, w) in target_expected.items())
       + f"\non the file's reference level {on_reference:.2f}, expected 329.46, difference "
       f"{on_reference - 329.46:+.3f}\nlargest departure "
       f"{max(worst_target, abs(on_reference - 329.46)):.3f} m/s, bound 0.05\n"
       f"10 N is a node of the file's 0.5 degree grid, which is why the read-back matches the "
       f"construction here and not at the anchor (check 11)")

# ---------------------------------------------------------------------------
# 11. Beyond the specification: the grid read-back against the construction
# ---------------------------------------------------------------------------
gaps = {name: abs(float(read_back[i]) - float(construction[i])) for name, (i, _) in expected.items()}
record(11, "beyond the specification: what the 0.5 degree file grid costs at the anchor, which is "
           "what Step 4 will read",
       True,
       "\n".join(f"{name:6s} construction {float(construction[i]):10.6f}, read back through "
                 f"wind_at {float(read_back[i]):10.6f}, difference {gaps[name]:+.6f}"
                 for name, (i, _) in expected.items())
       + f"\nlargest {max(gaps.values()):.3f} m/s, at the gauge level\n"
       f"the anchor at {math.degrees(phi_a):.4f} degrees falls between the file's "
       f"{math.degrees(lat_file[int(np.searchsorted(lat_file, phi_a)) - 1]):.1f} and "
       f"{math.degrees(lat_file[int(np.searchsorted(lat_file, phi_a))]):.1f} degree nodes, on the "
       f"steepest interval of the curve; decision P's expected values are the construction, and "
       f"this is the grid interpolation the decision says Steps 3 and 4 bound. A description, "
       f"not a bound")

# ---------------------------------------------------------------------------
# 12. Beyond the specification: the three-point difference on unequal spacing
# ---------------------------------------------------------------------------
quadratic = 3.0 * mesh.latitude_rad ** 2 - 2.0 * mesh.latitude_rad + 1.0
slope = mesh.d_dlatitude(quadratic)
analytic = 6.0 * mesh.latitude_rad - 2.0
uneven = int(np.sum(~np.isclose(np.diff(mesh.latitude_rad), LATITUDE_SPACING, rtol=1e-9)))
# The comparison is absolute against the scale of the derivative, not relative pointwise:
# 6 phi - 2 passes through zero at phi = 1/3 rad, which is inside this mesh, and a pointwise
# relative test there divides by nothing.
absolute = float(np.max(np.abs(slope - analytic)))
scale = float(np.max(np.abs(analytic)))
record(12, "beyond the specification: the centered difference in latitude is exact on a quadratic, "
           "including at the cells the inserted nodes split",
       absolute <= 1e-12 * scale,
       f"max absolute departure from 6 phi - 2: {absolute:.3e}, against {1e-12 * scale:.3e}, "
       f"which is 1e-12 of the derivative's range {scale:.4f}; the derivative passes through zero "
       f"inside this mesh, at phi = 1/3 radian, so a pointwise relative test is not the measure\n"
       f"{uneven} of {mesh.latitude_rad.size - 1} cells are not the declared spacing, the ones "
       f"the anchor and the target split; a formula assuming uniform spacing would be first order "
       f"there, and the three-point formula for unequal spacing is exact on a quadratic everywhere")

# ---------------------------------------------------------------------------
print()
passed = sum(1 for _, _, ok, _ in results if ok)
print(f"{passed} of {len(results)} checks pass")
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n"
              + "\n".join("        " + line for line in str(detail).splitlines())
              for n, d, ok, detail in results) + f"\n\n{passed} of {len(results)} checks pass\n",
    encoding="utf-8")
sys.exit(0 if passed == len(results) else 1)
