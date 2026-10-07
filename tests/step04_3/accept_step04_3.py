"""Acceptance checks for SPEC_04 v0.12 Step 3: `lib.kernel`, the shear kernel, its vertical
integral and the transfer kernel.

Every check prints its measured value; checks beyond the specification are labeled so. The anchor
and the run's inputs are the clean swept products of the Step 2 acceptance commit `562c769`, so
nothing is relaxed and no `-dirty` refusal is touched. This step registers nothing and writes
nothing inside the repository outside `reports/step04_3/`, which git ignores.

**How the synthetic wind states are posed.** `u = beta r sin(phi)` is `u = beta Z`, so
`(du/dZ)_R = beta` exactly, and `u = dOmega r cos(phi)` is `u = dOmega R`, so `(du/dZ)_R = 0`. Both
identities are kinematic: they hold for `u` as a function of position, so the map that defines `u`
and the geometry the kernel differentiates on have to be the same map. The synthetic files are
therefore built on the mesh's own columns, on a grid that covers the mesh and nothing more, and are
held in memory rather than written: a field that is nonzero at the poles is not a wind the model
runs on, and has no business on disk.

**The map.** Both directions of the anchor's flat-isobar map are continued at the slope of their
last interval outside the anchor's levels, as Step 2 decision 8 states. The first filing continued
`Phi(p)` and clamped its inverse, which put several mesh nodes at one pressure: a slope
discontinuity decision 8 does not have, two coinciding wind-file pressure nodes, and a kernel
residual an order of magnitude larger at the ends. SPEC_04 section 14 ruling 2 identifies that as
the script's, not the specification's, and it is corrected here.

**What check 9 measures.** The transfer applies the line integral of the kernel along an isobar,
not its maximum (decision Q). The residual of a gridded wind is spiky and alternates in sign inside
the cells holding a slope jump, so its maximum does not fall under refinement while its integral
does, and bounding `d ln N` by the maximum times the span overstates it by two orders of magnitude.
The first filing did exactly that.
"""

import math
import sys
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.forward import production as fp
from casspian.lib import control as ctl
from casspian.lib import geoid
from casspian.lib import kernel as lk
from casspian.lib import mesh as lm
from casspian.lib import windfield as wf
from casspian.lib.gravity import omega_abs
from casspian.lib.reduction import mean_over_species

HERE = Path("reports/step04_3")
CLOSURE = Path("forward/lindal_closure")
LATITUDE_SPACING = math.radians(0.05)
GEOPOTENTIAL_SPACING = 5.0e3

HERE.mkdir(parents=True, exist_ok=True)
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}", flush=True)
    for line in str(detail).splitlines():
        print(f"        {line}", flush=True)


def worst_absolute(got, want):
    return float(np.max(np.abs(np.asarray(got, dtype="float64")
                               - np.asarray(want, dtype="float64"))))


# ---------------------------------------------------------------------------
# The run, the anchor, the closure geometry
# ---------------------------------------------------------------------------
namelist = ctl.read_run_namelist(CLOSURE / "lindal_closure.toml")
inputs = ctl.load_run_inputs(namelist)
anchor = inputs.anchor
anchor_root = anchor.to_dataset(inherit=False)

phi_a = math.radians(float(anchor_root["latitude_planetocentric_deg"].values))
r_a = float(anchor_root["anchor_isobar_radius_m"].values)
p_tab = np.asarray(anchor["inputs"]["thermo"]["pressure_Pa"].values, dtype="float64")
gauge_level = int(np.flatnonzero(p_tab == namelist.gauge_isobar_Pa)[0])

Omega = float(inputs.rotation["angular_rate_rad_s"])
CONSTANTS = (Omega, float(inputs.gravity["GM_m3s2"]), np.asarray(inputs.gravity["J"].values),
             np.asarray(inputs.gravity["degree"].values),
             float(inputs.gravity["normalization_radius_m"]))

closure_field = wf.WindField(inputs.wind)
lat_file = np.radians(np.asarray(inputs.wind["latitude_planetocentric_deg"].values,
                                 dtype="float64"))
p_file = np.asarray(inputs.wind["pressure_Pa"].values, dtype="float64")
u_reference_file = np.asarray(inputs.wind["u_reference_ms"].values, dtype="float64")
p_ref = float(inputs.wind["reference_level_pressure_Pa"])

profile = fp.profile_from_anchor(anchor)
p_b = float(p_tab[0])
closure_column = np.asarray(
    closure_field.wind_at(np.full(p_tab.shape, phi_a), p_tab), dtype="float64")
Phi_closure = fp.produce(profile, inputs, gauge_level, p_b,
                         u_column=closure_column).geopotential.geopotential_m2s2

_rank = np.argsort(Phi_closure)
_order = np.argsort(np.log(p_tab))
_x, _y = np.log(p_tab)[_order], Phi_closure[_order]
_low = (_y[1] - _y[0]) / (_x[1] - _x[0])
_high = (_y[-1] - _y[-2]) / (_x[-1] - _x[-2])


def Phi_of_p(p):
    """The anchor's flat-isobar map, linear in `ln p`, continued at the last slope outside."""
    lp = np.log(np.asarray(p, dtype="float64"))
    return np.where(lp < _x[0], _y[0] + _low * (lp - _x[0]),
                    np.where(lp > _x[-1], _y[-1] + _high * (lp - _x[-1]), np.interp(lp, _x, _y)))


def p_of_Phi(Phi):
    """The inverse of the flat-isobar map, continued at the slope of its last interval.

    `np.interp` clamps, and the mesh reaches one geopotential spacing beyond the anchor's levels
    at each end, so a clamped inverse puts several nodes at one pressure: a slope discontinuity
    where decision 8 has none, two coinciding wind-file pressure nodes, and a spurious kernel
    residual an order of magnitude larger at the ends (SPEC_04 section 14 ruling 2). The
    continuation is the inverse of `Phi_of_p`, which continues rather than clamps, so the two are
    one map.
    """
    Phi = np.asarray(Phi, dtype="float64")
    y, x = Phi_closure[_rank], np.log(p_tab)[_rank]
    below = (x[1] - x[0]) / (y[1] - y[0])
    above = (x[-1] - x[-2]) / (y[-1] - y[-2])
    ln_p = np.where(Phi < y[0], x[0] + below * (Phi - y[0]),
                    np.where(Phi > y[-1], x[-1] + above * (Phi - y[-1]), np.interp(Phi, y, x)))
    return np.exp(ln_p)


def on_file_grid(Phi):
    """`p_of_Phi` held inside the wind file's own pressure range.

    The column node set for the cylinder construction spans the file's whole pressure range, and
    the round trip `Phi_of_p` then `p_of_Phi` is not exact to the bit at its ends, so the pressure
    can land one ulp outside the grid and `lib.windfield` refuses it, correctly. This is a clamp of
    one ulp and nothing else: the excursion is asserted below, and it is not the structural clamp
    of the first filing, which moved whole mesh rows onto one pressure.
    """
    p = p_of_Phi(Phi)
    low, high = float(p_file[0]), float(p_file[-1])
    excursion = max(float(np.max(low - np.minimum(p, low))) / low,
                    float(np.max(np.maximum(p, high) - high)) / high)
    assert excursion < 1e-12, f"p_of_Phi leaves the file's grid by {excursion:.3e} relative"
    return np.clip(p, low, high)


def mesh_for(target_deg):
    """The working mesh and its closure columns for a run to one target latitude."""
    target = math.radians(target_deg)
    mesh = lm.build_mesh([phi_a, target], Phi_closure, LATITUDE_SPACING, GEOPOTENTIAL_SPACING)
    r0 = geoid.through_anchor(
        mesh.latitude_rad, phi_a, r_a,
        lambda x: closure_field.wind_at(np.asarray(x, dtype="float64"),
                                        namelist.gauge_isobar_Pa), *CONSTANTS).radius
    columns = lm.build_columns(
        mesh, r0,
        (lambda i: (lambda Phi, la=float(mesh.latitude_rad[i]):
                    closure_field.wind_at(la, p_of_Phi(Phi)))), *CONSTANTS)
    radius = np.stack([c.radius_m for c in columns])
    g = np.stack([c.g_radial_ms2 for c in columns])
    pressure = p_of_Phi(mesh.geopotential_m2s2)[None, :] * np.ones_like(radius)
    return mesh, radius, g, pressure, target


mesh10, R10, G10, P10, target10 = mesh_for(10.0)
lat10 = mesh10.latitude_rad[:, None]

# ---------------------------------------------------------------------------
# 1. The closure wind: S / g is the closed form, no difference taken
# ---------------------------------------------------------------------------
U10 = np.asarray(closure_field.wind_at(np.broadcast_to(lat10, R10.shape), P10), dtype="float64")
S10, dudZ10 = lk.shear_kernel(mesh10, R10, P10, G10, U10, closure_field, Omega)
s_over_g10 = S10 / G10
u_prime, du_dlnp = closure_field.wind_derivatives(np.broadcast_to(lat10, R10.shape), P10)
closed = 2.0 * omega_abs(U10, R10, lat10, Omega) * np.cos(lat10) * u_prime / G10
record(1, "under the closure wind S / g is 2 Omega_abs cos(phi) u'(phi) / g at every node, to "
          "round-off, since no difference is taken",
       worst_absolute(s_over_g10, closed) <= 1e-15 and float(np.max(np.abs(du_dlnp))) == 0.0,
       f"the file carries no shear: max |du/dln p| over the mesh is "
       f"{float(np.max(np.abs(du_dlnp))):.1e}, exactly zero, so both conversion terms multiply "
       f"zero and no mesh difference reaches S\n"
       f"largest |S/g - 2 Omega_abs cos(phi) u'/g| over {s_over_g10.size} nodes: "
       f"{worst_absolute(s_over_g10, closed):.3e}\n"
       f"the residual is the r (cos phi / r) round trip in S = 2 Omega_abs r (cos phi / r) u' "
       f"against 2 Omega_abs cos(phi) u', not an approximation")

# ---------------------------------------------------------------------------
# 2. The largest |S/g| over each run's span
# ---------------------------------------------------------------------------
span10 = (mesh10.latitude_rad >= target10) & (mesh10.latitude_rad <= phi_a)
largest10 = float(np.max(np.abs(s_over_g10[span10])))

mesh60, R60, G60, P60, target60 = mesh_for(60.0)
lat60 = mesh60.latitude_rad[:, None]
U60 = np.asarray(closure_field.wind_at(np.broadcast_to(lat60, R60.shape), P60), dtype="float64")
S60, _ = lk.shear_kernel(mesh60, R60, P60, G60, U60, closure_field, Omega)
span60 = (mesh60.latitude_rad >= phi_a) & (mesh60.latitude_rad <= target60)
largest60 = float(np.max(np.abs((S60 / G60)[span60])))
record(2, "the largest |S/g| on the mesh is 0.0977 per radian from phi_c to 10 N and 0.0540 to 60 N",
       abs(largest10 - 0.0977) <= 5e-4 and abs(largest60 - 0.0540) <= 5e-4,
       f"to 10 N: {largest10:.6f} per radian, expected 0.0977, difference {largest10 - 0.0977:+.6f}\n"
       f"to 60 N: {largest60:.6f} per radian, expected 0.0540, difference {largest60 - 0.0540:+.6f}\n"
       f"the 10 N mesh has {mesh10.latitude_rad.size} latitudes and the 60 N mesh "
       f"{mesh60.latitude_rad.size}, both at 0.05 degrees")

# ---------------------------------------------------------------------------
# 3. I is zero at the Phi = 0 node of every column
# ---------------------------------------------------------------------------
I10 = lk.shear_integral(mesh10, s_over_g10)
at_gauge = I10[:, mesh10.gauge_geopotential_index]
record(3, "I = 0 at the Phi = 0 node of every column exactly",
       bool(np.all(at_gauge == 0.0)),
       f"{at_gauge.size} columns, every one exactly zero: {bool(np.all(at_gauge == 0.0))}; "
       f"largest |I| there {float(np.max(np.abs(at_gauge)))!r}\n"
       f"I runs {float(I10.min()):+.4e} to {float(I10.max()):+.4e} over the mesh; the sums are "
       f"taken outward from the gauge in both directions, as lib.geopotential takes the "
       f"geopotential")

# ---------------------------------------------------------------------------
# The synthetic wind states, posed on the mesh's own geometry
# ---------------------------------------------------------------------------
def synthetic_field(mesh, radius, u_of_radius_latitude):
    """A `WindField` carrying an analytic `u(r, phi)` on the mesh's own columns.

    The grid is the mesh's own latitudes and the anchor's own levels, so `r` at a file node is the
    mesh's `r` interpolated to that level and the identity the check rests on is posed on the
    geometry the kernel differentiates on. Held in memory: the field is nonzero at the poles and
    does not span them, so it is not a wind the model runs on and has no business on disk.

    The grid is the mesh's own geopotential nodes mapped through `p_of_Phi`, which is strictly
    monotonic now that the map is continued rather than clamped (SPEC_04 section 14 ruling 2), so
    every node of the mesh is a node of the file and none is clipped.
    """
    pressure = p_of_Phi(mesh.geopotential_m2s2)
    rising = np.argsort(pressure)
    u = u_of_radius_latitude(radius, mesh.latitude_rad[:, None])
    data = xr.Dataset(
        {
            "u_total_ms": (("latitude_planetocentric", "pressure"), u[:, rising]),
            "latitude_planetocentric_deg": (("latitude_planetocentric",),
                                            np.degrees(mesh.latitude_rad)),
            "pressure_Pa": (("pressure",), pressure[rising]),
            "reference_level_pressure_Pa": ((), float(pressure[rising][0])),
        }
    )
    return wf.WindField(data)


BETA = 1.0e-6
P_on = P10
beta_field = synthetic_field(mesh10, R10, lambda r, phi: BETA * r * np.sin(phi))
U_beta = np.asarray(beta_field.wind_at(np.broadcast_to(lat10, R10.shape), P_on), dtype="float64")
S_beta, dudZ_beta = lk.shear_kernel(mesh10, R10, P_on, G10, U_beta, beta_field, Omega)
want_beta = 2.0 * omega_abs(U_beta, R10, lat10, Omega) * R10 * BETA

# With the map continued there is no slope discontinuity at the anchor's ends, so no node set is
# set aside: the whole mesh is measured.
ratio = np.abs(dudZ_beta / BETA - 1.0)
record(4, "a synthetic u = beta r sin(phi) has (du/dZ)_R = beta and S = 2 Omega_abs r beta to the "
          "truncation error of the differences",
       float(np.max(ratio)) <= 5e-2,
       f"beta = {BETA:g} per second, posed on the mesh's own columns\n"
       f"(du/dZ)_R / beta - 1 over all {ratio.size} nodes: median "
       f"{float(np.median(ratio)):.3e}, 99th percentile {float(np.quantile(ratio, 0.99)):.3e}, "
       f"largest {float(np.max(ratio)):.3e}\n"
       f"the map is continued at the slope of its last interval (Step 2 decision 8), so there is "
       f"no slope discontinuity at the ends of the anchor's levels and no node set is set aside; "
       f"the first filing clamped it instead and measured 4.203e-01 at those ends\n"
       f"S against 2 Omega_abs r beta: largest relative "
       f"{float(np.max(np.abs((S_beta - want_beta) / want_beta))):.3e}, which is the truncation "
       f"of the file's linear interpolant and is reported, not bounded (SPEC_04 v0.12)")

DELTA_OMEGA = 1.0e-6
solid_field = synthetic_field(mesh10, R10, lambda r, phi: DELTA_OMEGA * r * np.cos(phi))
U_solid = np.asarray(solid_field.wind_at(np.broadcast_to(lat10, R10.shape), P_on), dtype="float64")
S_solid, dudZ_solid = lk.shear_kernel(mesh10, R10, P_on, G10, U_solid, solid_field, Omega)
solid_worst = float(np.max(np.abs(S_solid / G10)))
record(5, "a solid-body offset u = dOmega r cos(phi) has S = 0 analytically, and on the mesh "
          "|S/g| is below 1e-4 per radian (SPEC_04 v0.12)",
       solid_worst <= 1e-4,
       f"dOmega = {DELTA_OMEGA:g} per second, u = dOmega R, a function of the cylindrical radius "
       f"alone, so (du/dZ)_R vanishes for any r and phi\n"
       f"largest |S/g| over the {S_solid.size} nodes of the mesh: {solid_worst:.3e} per radian, "
       f"bound 1e-4\n"
       f"(du/dZ)_R largest {float(np.max(np.abs(dudZ_solid))):.3e} against dOmega "
       f"{DELTA_OMEGA:.1e}, a ratio of {float(np.max(np.abs(dudZ_solid))) / DELTA_OMEGA:.3e}; "
       f"that ratio is the truncation of the file's linear interpolant, not round-off, which is "
       f"what v0.12 restates")

# ---------------------------------------------------------------------------
# 6 and 7. The composition term
# ---------------------------------------------------------------------------
species = inputs.composition["species"].to_dataset(inherit=False)
names = [str(n) for n in species["species_name"].values]
R_i = dict(zip(names, np.asarray(species["refractivity_per_molecule_m3"].values,
                                 dtype="float64").tolist()))
M_i = dict(zip(names, np.asarray(species["molar_mass_kg_mol"].values, dtype="float64").tolist()))

comp_root = inputs.composition.to_dataset(inherit=False)
comp_p = np.asarray(comp_root["pressure_Pa"].values, dtype="float64")
fine_lat = np.arange(-90.0, 90.0 + 0.005, 0.01)
x_He = 0.06 + 0.02 * np.sin(np.radians(fine_lat))
x_H2 = 1.0 - x_He
synthetic_comp = xr.DataTree()
synthetic_comp.dataset = xr.Dataset({
    "x_H2": (("level", "latitude_planetocentric"),
             np.tile(x_H2, (comp_p.size, 1))),
    "x_He": (("level", "latitude_planetocentric"),
             np.tile(x_He, (comp_p.size, 1))),
    "pressure_Pa": (("level",), comp_p),
    "latitude_planetocentric_deg": (("latitude_planetocentric",), fine_lat),
})
synthetic_comp["species"] = xr.DataTree(xr.Dataset({
    "species_name": (("species",), np.array(["H2", "He"])),
    "refractivity_per_molecule_m3": (("species",), np.array([R_i["H2"], R_i["He"]])),
    "molar_mass_kg_mol": (("species",), np.array([M_i["H2"], M_i["He"]])),
}))

labels = np.array([p_tab[0], namelist.gauge_isobar_Pa, p_tab[-1]])
latitude_c, term = lk.composition_term(synthetic_comp, labels)
phi_c_grid = latitude_c
xh = 0.06 + 0.02 * np.sin(phi_c_grid)
dx = 0.02 * np.cos(phi_c_grid)
R_bar = (1.0 - xh) * R_i["H2"] + xh * R_i["He"]
m_bar = (1.0 - xh) * M_i["H2"] + xh * M_i["He"]
analytic = dx * (R_i["He"] - R_i["H2"]) / R_bar - dx * (M_i["He"] - M_i["H2"]) / m_bar
inner = slice(1, -1)
worst_comp = worst_absolute(term[:, inner], np.broadcast_to(analytic[inner], term[:, inner].shape))
record(6, "a synthetic composition with x_He = 0.06 + 0.02 sin(phi) and H2 as the closure gives "
          "the composition term equal to its closed form to 1e-8",
       worst_comp <= 1e-8,
       f"largest departure {worst_comp:.3e} over {term[:, inner].size} points, "
       f"{labels.size} isobar labels by {int(phi_c_grid.size) - 2} interior latitudes, bound 1e-8\n"
       f"the closed form is d ln(R_bar / m_bar) / dphi with dx_He/dphi = 0.02 cos(phi), evaluated "
       f"from the species table: R_H2 {R_i['H2']:.6e}, R_He {R_i['He']:.6e} m3, m_H2 "
       f"{M_i['H2']:.6e}, m_He {M_i['He']:.6e} kg/mol\n"
       f"the term runs {float(term.min()):+.6e} to {float(term.max()):+.6e} per radian; the "
       f"latitude grid is {float(np.degrees(np.mean(np.diff(phi_c_grid)))):.3f} degrees, chosen "
       f"fine enough that the centered difference meets the bound (the two edge latitudes are "
       f"one-sided and are excluded)")

latitude_t, term_t = lk.composition_term(inputs.composition, labels)
varies = [name for name, v in inputs.composition.to_dataset(inherit=False).data_vars.items()
          if "latitude_planetocentric" in v.dims and not np.all(v.values == v.values[..., :1])]
record(7, "the run's composition, whose columns are identical, gives a composition term of exactly "
          "zero",
       bool(np.all(term_t == 0.0)),
       f"largest |term| {float(np.max(np.abs(term_t)))!r} over {term_t.size} points, exactly "
       f"zero: {bool(np.all(term_t == 0.0))}\n"
       f"the file varies with latitude in {varies or 'nothing'}, so the difference is of equal "
       f"numbers and not of nearly equal ones")

# ---------------------------------------------------------------------------
# 8. Beyond the specification: transfer_kernel on the mesh's own nodes
# ---------------------------------------------------------------------------
probe_lat = mesh10.latitude_rad[[5, 100, 300]]
probe_Phi = mesh10.geopotential_m2s2[[7, 200, 600]]
on_node = lk.transfer_kernel(mesh10, s_over_g10, probe_lat[:, None], probe_Phi[None, :])
direct = s_over_g10[np.ix_([5, 100, 300], [7, 200, 600])]
with_comp = lk.transfer_kernel(mesh10, s_over_g10, probe_lat[:, None], probe_Phi[None, :],
                               latitude_c, term[1])
record(8, "beyond the specification: transfer_kernel returns S / g at a node and adds the "
          "composition term off it",
       worst_absolute(on_node, direct) == 0.0,
       f"at nine mesh nodes the bilinear form returns the node's own S / g exactly: "
       f"largest departure {worst_absolute(on_node, direct):.1e}\n"
       f"with the synthetic composition term added, K differs from S / g by "
       f"{float(np.max(np.abs(with_comp - on_node))):.6e} per radian, which is the composition "
       f"term there\n"
       f"K = S / g + (d ln(R_bar / m_bar) / dphi)_p is Eq. A27")

# ---------------------------------------------------------------------------
# 9. The cylinder-extended wind: S is zero analytically
# ---------------------------------------------------------------------------
# `u = f(R) ` is a function of the cylindrical radius alone, so the two terms of `(du/dZ)_R`
# cancel for any `f` and `S` vanishes identically. On the mesh the residual comes from the wind
# file's grid: `U(s)` is stored on latitude and pressure nodes and read back by the decision L
# interpolant, so the field the kernel sees is a piecewise-bilinear approximation of one that is
# constant on cylinders. The residual is reported as the largest `|S/g|` and the `|d ln N|` it
# would put on the run's span, and is measured at three pressure resolutions.
#
# The construction is decision P's, built here rather than read from Step 2's output so that this
# suite stands alone. The inversion curve and the columns do not depend on the file's pressure
# grid, so the fixed point is converged once and re-sampled onto each grid.
fine = np.arange(float(lat_file[0]), float(lat_file[-1]) + 0.5 * LATITUDE_SPACING,
                 LATITUDE_SPACING)
lat_inv = np.unique(np.clip(np.concatenate([fine, lat_file, [phi_a, target10]]),
                            float(lat_file[0]), float(lat_file[-1])))
on_file = np.searchsorted(lat_inv, lat_file)
ref_col = int(np.flatnonzero(p_file == p_ref)[0])

Phi_k, field_now, passes = Phi_closure.copy(), closure_field, 0
while True:
    Phi_p = Phi_of_p(p_file)
    nodes = np.unique(np.concatenate([Phi_p, Phi_k, [0.0]]))
    at_pressure = np.searchsorted(nodes, Phi_p)
    r0 = geoid.through_anchor(
        lat_inv, phi_a, r_a,
        lambda x: field_now.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
        *CONSTANTS).radius
    radius_inv = np.empty((lat_inv.size, p_file.size))
    for i, latitude in enumerate(lat_inv):
        radius_inv[i] = lm.column(
            latitude, nodes, float(r0[i]),
            (lambda Phi, la=latitude: field_now.wind_at(la, on_file_grid(Phi))),
            *CONSTANTS).radius_m[at_pressure]
    s_reference = radius_inv[:, ref_col] * np.cos(lat_inv)

    def branch(mask):
        s, q = s_reference[mask], lat_inv[mask]
        rising = np.argsort(s)
        return s[rising], q[rising]

    s_north, phi_north = branch(lat_inv >= 0.0)
    s_south, phi_south = branch(lat_inv <= 0.0)
    s_grid = radius_inv[on_file] * np.cos(lat_file)[:, None]
    in_north = lat_file >= 0.0
    phi_star = np.empty_like(s_grid)
    phi_star[in_north] = np.interp(s_grid[in_north], s_north, phi_north)
    phi_star[~in_north] = np.interp(s_grid[~in_north], s_south, phi_south)
    u_total = np.interp(phi_star, lat_file, u_reference_file)
    u_total[0], u_total[-1] = 0.0, 0.0
    change = float(np.max(np.abs(u_total - field_now.u_total_ms)))
    passes += 1
    built = inputs.wind.copy(deep=True)
    built["u_total_ms"] = (("latitude_planetocentric", "pressure"), u_total,
                           dict(inputs.wind["u_total_ms"].attrs))
    field_now = wf.WindField(built)
    Phi_k = fp.produce(profile, inputs, gauge_level, p_b, u_column=np.asarray(
        field_now.wind_at(np.full(p_tab.shape, phi_a), p_tab),
        dtype="float64")).geopotential.geopotential_m2s2
    if change < 1.0e-6 or passes >= 10:
        break


# The columns are integrated once, on a node set fine enough for every grid below, and the radius
# at any pressure is then read off by interpolation in Phi. The columns do not depend on the file's
# grids, so re-integrating them per grid would only repeat work.
def log_grid(per_decade):
    """A log spaced pressure grid over the file's own range, its ends pinned to the file's.

    `np.logspace` undershoots by an ulp, and the file's own grid is what the cylinder field is
    read from, so an endpoint a hair outside it is refused (correctly) by `lib.windfield`.
    """
    span = math.log10(float(p_file[-1]) / float(p_file[0]))
    grid = np.logspace(math.log10(float(p_file[0])), math.log10(float(p_file[-1])),
                       int(round(span * per_decade)) + 1)
    grid[0], grid[-1] = float(p_file[0]), float(p_file[-1])
    return grid


_finest = log_grid(40)
_nodes = np.unique(np.concatenate([Phi_of_p(_finest), Phi_k, [0.0]]))
_r0 = geoid.through_anchor(
    lat_inv, phi_a, r_a,
    lambda x: field_now.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
    *CONSTANTS).radius
_radius = np.empty((lat_inv.size, _nodes.size))
for _i, _la in enumerate(lat_inv):
    _radius[_i] = lm.column(
        _la, _nodes, float(_r0[_i]),
        (lambda Phi, la=_la: field_now.wind_at(la, on_file_grid(Phi))), *CONSTANTS).radius_m


def cylinder_on(pressures, latitudes):
    """The converged field re-sampled onto a latitude by pressure grid, as a `WindField`.

    `latitudes` is a subset of the inversion latitudes, so the radius at each of its nodes is one
    of the rows integrated above and the field is posed on one geometry throughout.
    """
    rows = np.searchsorted(lat_inv, latitudes)
    radius = np.stack([np.interp(Phi_of_p(pressures), _nodes, _radius[i]) for i in rows])
    s = radius * np.cos(latitudes)[:, None]
    north = latitudes >= 0.0
    star = np.empty_like(s)
    star[north] = np.interp(s[north], s_north, phi_north)
    star[~north] = np.interp(s[~north], s_south, phi_south)
    u = np.interp(star, lat_file, u_reference_file)
    u[0], u[-1] = 0.0, 0.0
    data = xr.Dataset({
        "u_total_ms": (("latitude_planetocentric", "pressure"), u),
        "latitude_planetocentric_deg": (("latitude_planetocentric",), np.degrees(latitudes)),
        "pressure_Pa": (("pressure",), pressures),
        "reference_level_pressure_Pa": ((), p_ref),
    })
    return wf.WindField(data)


decades = math.log10(float(p_file[-1]) / float(p_file[0]))


def latitudes_at(step_deg):
    """A subset of the inversion latitudes at the named spacing, spanning the file's range.

    Both poles are kept: the reference surface is marched to them, so a field that stops short of
    90 degrees is asked for a latitude it does not cover and refuses, which is `lib.windfield`
    doing its job.
    """
    want = np.radians(np.clip(np.arange(-90.0, 90.0 + 0.5 * step_deg, step_deg), -90.0, 90.0))
    nearest = np.abs(lat_inv[:, None] - want[None, :]).argmin(axis=0)
    index = np.unique(np.concatenate([nearest, [0, lat_inv.size - 1]]))
    return lat_inv[index]


lines, residuals = [], {}
for per_decade, step_deg in ((10, 0.5), (20, 0.5), (40, 0.5), (10, 0.25), (10, 0.1)):
    grid = log_grid(per_decade)
    file_lat = latitudes_at(step_deg)
    field_c = cylinder_on(grid, file_lat)
    # The mesh has to be the cylinder wind's own: `u = U(s)` is constant on cylinders of the
    # surface the cylinder field defines, and the closure columns put a different radius at the
    # same node, so reading this field on the closure mesh asks it for U at the wrong cylinder.
    # That is the same posing error the synthetic states carry, and it is worth 0.2 per radian.
    r0_c = geoid.through_anchor(
        mesh10.latitude_rad, phi_a, r_a,
        lambda x: field_c.wind_at(np.asarray(x, dtype="float64"), namelist.gauge_isobar_Pa),
        *CONSTANTS).radius
    columns_c = lm.build_columns(
        mesh10, r0_c,
        (lambda i: (lambda Phi, la=float(mesh10.latitude_rad[i]):
                    field_c.wind_at(la, on_file_grid(Phi)))), *CONSTANTS)
    R_c = np.stack([c.radius_m for c in columns_c])
    G_c = np.stack([c.g_radial_ms2 for c in columns_c])
    U_c = np.stack([c.u_ms for c in columns_c])
    on_grid = np.clip(P10, grid.min(), grid.max())
    S_c, _ = lk.shear_kernel(mesh10, R_c, on_grid, G_c, U_c, field_c, Omega)
    worst = float(np.max(np.abs((S_c / G_c)[span10])))
    # What the transfer applies is the line integral of the kernel along an isobar, not its
    # maximum (SPEC_04 v0.12 decision Q). The isobars are flat on this pass, so isobar k is the
    # line Phi = Phi_k, and `d ln N` on it is the integral of S/g in latitude from the anchor to
    # the target. The residual alternates in sign inside the cells holding a slope jump, so the
    # integral is far smaller than the maximum times the span, which is what the first filing
    # reported and what the review corrects.
    span_lat = mesh10.latitude_rad[span10]
    ascending = np.argsort(span_lat)
    on_isobar = lk.transfer_kernel(mesh10, S_c / G_c, span_lat[ascending, None],
                                   Phi_closure[None, :])
    I_c = lk.shear_integral(mesh10, S_c / G_c)
    I_isobar = lk.transfer_kernel(mesh10, I_c, span_lat[ascending, None], Phi_closure[None, :])
    # From the anchor to the target, which runs equatorward, so the integral is negated.
    d_ln_N_levels = -np.trapezoid(on_isobar, span_lat[ascending], axis=0)
    shift_levels = np.trapezoid(I_isobar, span_lat[ascending], axis=0)
    d_ln_N = float(np.max(np.abs(d_ln_N_levels)))
    shift = float(np.max(np.abs(shift_levels)))
    residuals[(per_decade, step_deg)] = (worst, d_ln_N, shift, d_ln_N_levels)
    lines.append(f"{step_deg:5.2f} deg by {per_decade:3d} per decade "
                 f"({file_lat.size:5d} latitudes, {grid.size:4d} pressures): largest |d ln N| "
                 f"{d_ln_N:.3e}, largest |shift| {shift:6.1f} m2/s2, largest |S/g| "
                 f"{worst:.3e} per radian")

base_worst, base_d_ln_N, base_shift, base_levels = residuals[(10, 0.5)]
one_bar = int(np.argmin(np.abs(p_tab - 1.0e5)))
ten_mbar = int(np.argmin(np.abs(p_tab - 1.0e3)))
stated = {"the top level": (0, 1.2e-4), "10 mbar": (ten_mbar, -5.2e-4),
          "the gauge": (gauge_level, -1.7e-4), "the bottom": (p_tab.size - 1, -1.1e-3)}
by_latitude = [residuals[(10, s)][1] for s in (0.5, 0.25, 0.1)]
record(9, "under the cylinder-extended wind the line integral of S / g along each of the anchor's "
          "isobars is below 2e-3 in d ln N and the isobar shift below 300 m2/s2 (SPEC_04 v0.12 "
          "decision Q)",
       base_d_ln_N <= 2e-3 and base_shift <= 300.0,
       f"the fixed point converged in {passes} passes\n" + "\n".join(lines)
       + f"\non the file's own 0.5 degree grid: largest |d ln N| {base_d_ln_N:.3e}, bound 2e-3; "
       f"largest |shift| {base_shift:.1f} m2/s2, bound 300\n"
       + "\n".join(f"    {name:14s} d ln N {float(base_levels[i]):+.2e}, stated {want:+.1e}"
                   for name, (i, want) in stated.items())
       + f"\nthe largest |S/g| is {base_worst:.3e} per radian, reported and not bounded: the "
       f"residual is spiky and alternates in sign inside the cells holding a slope jump, so its "
       f"maximum does not fall under refinement while its integral does\n"
       f"the line integral with the cylinder construction resampled at 0.5, 0.25 and 0.1 degrees, "
       f"the mesh unchanged: " + ", ".join(f"{v:.3e}" for v in by_latitude)
       + f", ratios {by_latitude[0] / by_latitude[1]:.2f} and {by_latitude[1] / by_latitude[2]:.2f}; "
       f"no order is claimed\n"
       f"this is the truncation floor of a wind hypothesis on a grid (decision Q), not a defect "
       f"and not a tension between decision L and A15: the report's first filing bounded d ln N "
       f"by the maximum times the span, which overstates a sign-alternating integrand")

# ---------------------------------------------------------------------------
print(flush=True)
passed = sum(1 for _, _, ok, _ in results if ok)
print(f"{passed} of {len(results)} checks pass", flush=True)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n"
              + "\n".join("        " + line for line in str(detail).splitlines())
              for n, d, ok, detail in results) + f"\n\n{passed} of {len(results)} checks pass\n",
    encoding="utf-8")
sys.exit(0 if passed == len(results) else 1)
