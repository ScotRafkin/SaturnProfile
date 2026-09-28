"""Acceptance checks for SPEC_02 v0.6 Step 3, `lib.reduction` and `refrac.reduce`.

Every check reports its measured value. The partials of the anchor terms are cross checked
against full reruns of the anchor fixed point on perturbed copies of the inputs, held in memory.
Nothing under `occul_data/` is written.
"""

import dataclasses
import math
import sys

import numpy as np

from casspian.lib import composition as comp
from casspian.lib import control as ctl
from casspian.lib import reduction as red
from casspian.lib.constants import BOLTZMANN_CONSTANT, CODATA_RELEASE
from casspian.refrac import reduce as rdc
from casspian.refrac.anchor import freeze_anchor

MANIFEST = "occul_data/lindal/lindal_reduction.toml"
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


manifest = ctl.read_reduction_manifest(MANIFEST)
# SPEC_02 v0.8 Step 6 may change the committed manifest's anchor rule. This suite tests Step 3
# as accepted, under mean polar anchoring, so it pins that rule and quantity in memory.
manifest = dataclasses.replace(manifest, anchor_rule="mean_polar_radius",
                               anchor_quantity="radius_polar_m")
inputs = ctl.load_reduction_inputs(manifest)
anchor = freeze_anchor(inputs, manifest)
r = rdc.reduce_profile(inputs, manifest, anchor)
p = r.pressure_Pa


def at(pressure):
    k = np.flatnonzero(p == pressure)
    if k.size != 1:
        raise RuntimeError(f"{pressure} Pa is not exactly one tabulated level")
    return int(k[0])


# SPEC_03 v0.5 Step 0: levels off the exact decades are on the declared grid, so they are found
# by the value Table I prints.
printed = np.asarray(inputs.thermo["pressure_printed_Pa"].values, dtype="float64")


def at_printed(pressure):
    k = np.flatnonzero(printed == pressure)
    if k.size != 1:
        raise RuntimeError(f"{pressure} Pa is not exactly one printed level")
    return int(k[0])


species = inputs.composition["species"].dataset
names = [str(v) for v in species["species_name"].values]
R_i = dict(zip(names, species["refractivity_per_molecule_m3"].values.tolist()))
# SPEC_04 decision O: kind C is a field on (level, latitude) for every use, and the reduction
# takes the column at the anchor's latitude by decision L. The checks below are on that column.
# The Lindal composition is uniform in latitude by construction, asserted here so that taking a
# column hides nothing.
_field = inputs.composition.to_dataset(inherit=False)
_varying = [name for name, v in _field.data_vars.items()
            if comp.LATITUDE_DIM in v.dims and not np.all(v.values == v.values[:, :1])]
assert not _varying, f"the composition varies with latitude in {_varying}; this suite reads "\
                     "the column the reduction reads and assumes the Lindal field is uniform"
composition_column = comp.column_at(inputs.composition,
                                    float(np.degrees(anchor.phi_c_rad)))
x_NH3 = np.asarray(composition_column.dataset["x_NH3"].values, dtype="float64")
R_DRY = 0.94 * R_i["H2"] + 0.06 * R_i["He"]

# ---------------------------------------------------------------------------
# 1. The 1 bar level
# ---------------------------------------------------------------------------
k1 = at(100000.0)
n1, N1 = float(r.number_density_m3[k1]), float(r.refractivity[k1])
expected_N1 = 5.3731235e25 * 4.836272e-30 * (1 - 10.9e-6)
record(1, "at 1 bar: n = 5.3731e25 m-3 and N = 2.59856e-4 with the tabulated 10.9 ppm of NH3",
       float(r.temperature_K[k1]) == 134.8 and abs(n1 - 5.373124e25) / 5.373124e25 < 1e-6
       and abs(N1 - 2.59856e-4) < 0.000005e-4,
       f"p = {p[k1]!r} Pa, T = {float(r.temperature_K[k1])!r} K, k_B = {BOLTZMANN_CONSTANT!r} "
       f"(CODATA {CODATA_RELEASE})\n"
       f"n = {n1:.7e} m-3 (5.3731e25; 5.373124e25 with CODATA 2018)\n"
       f"x_NH3 = {x_NH3[k1]!r}, R_bar = {float(r.mean_refractivity_m3[k1]):.7e} m3; the dry "
       f"value times (1 - 10.9e-6) is {4.836272e-30 * (1 - 10.9e-6):.7e}\n"
       f"N = {N1:.7e} (2.59856e-4; the spec's own product {expected_N1:.7e})")

# ---------------------------------------------------------------------------
# 2. The 794.33 mbar level, where NH3 is exactly zero
# ---------------------------------------------------------------------------
k2 = at_printed(79433.0)
ratio = float(r.refractivity[k2] / r.number_density_m3[k2])
record(2, "at 794.33 mbar, where NH3 is exactly zero, N / n equals the dry R_bar to round-off",
       x_NH3[k2] == 0.0 and abs(ratio - R_DRY) / R_DRY < 1e-15
       and abs(R_DRY - 4.836272e-30) < 5e-37,
       f"x_NH3 = {x_NH3[k2]!r}\n"
       f"N / n = {ratio!r} m3\n"
       f"dry R_bar = 0.94 R_H2 + 0.06 R_He from the species group = {R_DRY!r} (4.836272e-30)\n"
       f"relative difference {abs(ratio - R_DRY) / R_DRY:.2e}")

# ---------------------------------------------------------------------------
# 3. The top and bottom levels
# ---------------------------------------------------------------------------
kt, kb = at_printed(20.0), at(129848.0)
ks = at_printed(25.0)
nt, nb = float(r.number_density_m3[kt]), float(r.number_density_m3[kb])
ns = float(r.number_density_m3[ks])
factor = float(r.mean_refractivity_m3[kb]) / R_DRY
# SPEC_02 v0.10: on the grid, top n = 1.041934e22 and second level 1.270497e22.
record(3, "top n = 1.041934e22, second 1.270497e22, bottom n = 6.43287e25 with N from the NH3 "
          "reduced R_bar",
       float(r.temperature_K[kt]) == 138.7 and float(r.temperature_K[kb]) == 146.2
       and abs(nt - 1.041934e22) < 0.0000005e22 and abs(ns - 1.270497e22) < 0.0000005e22
       and abs(nb - 6.43287e25) < 0.000005e25
       and abs((1 - factor) - 79.3e-6) < 0.05e-6,
       f"top    p = {p[kt]!r} Pa, T = {float(r.temperature_K[kt])!r} K, n = {nt:.7e} (1.041934e22)\n"
       f"second p = {p[ks]!r} Pa, T = {float(r.temperature_K[ks])!r} K, n = {ns:.7e} (1.270497e22)\n"
       f"bottom p = {p[kb]!r} Pa, T = {float(r.temperature_K[kb])!r} K, n = {nb:.6e} (6.43287e25)\n"
       f"bottom x_NH3 = {x_NH3[kb]:.6e}, R_bar / dry = 1 - {1 - factor:.4e} (1 - 79.3e-6), "
       f"N = {float(r.refractivity[kb]):.7e}")

# ---------------------------------------------------------------------------
# 4. The radius at the anchor and at 1 bar
# ---------------------------------------------------------------------------
ka = at(manifest.anchor_isobar_Pa)
# SPEC_02 v0.10: Eq. B1 with the tilt, r = r0 + (h - h_ref) cos psi.
cos_psi = math.cos(anchor.psi_rad)
r1_expected = anchor.r0_m + (0.0 - r.h_ref_m) * cos_psi
rt_expected = anchor.r0_m + (float(r.height_m[0]) - r.h_ref_m) * cos_psi
record(4, "radius_m at the anchor level equals r0 exactly, at 1 bar r0 - h_ref cos psi, and at "
          "the top r0 + (h_top - h_ref) cos psi to 1e-6 m",
       r.radius_m[ka] == anchor.r0_m and r.radius_m[k1] == r1_expected
       and r.height_m[k1] == 0.0 and abs(r.radius_m[0] - rt_expected) < 1e-6,
       f"anchor level {ka} at {p[ka]:g} Pa: h_ref = {r.h_ref_m!r} m, cos psi = {cos_psi!r}\n"
       f"radius_m there {r.radius_m[ka]!r}, r0 {anchor.r0_m!r}, equal {r.radius_m[ka] == anchor.r0_m}\n"
       f"1 bar: h = {r.height_m[k1]!r}, radius_m {r.radius_m[k1]!r}, r0 - h_ref cos psi "
       f"{r1_expected!r}, equal {r.radius_m[k1] == r1_expected}\n"
       f"top: radius_m {r.radius_m[0]!r}, r0 + (h_top - h_ref) cos psi {rt_expected!r}, "
       f"difference {r.radius_m[0] - rt_expected:.3e} m (under the pinned mean_polar_radius rule)")

# ---------------------------------------------------------------------------
# 5. The fractional uncertainty of N
# ---------------------------------------------------------------------------
cN = r.companions["refractivity_uncertainty"]
frac = cN.value / r.refractivity
closed_form = 0.03 * (1 - x_NH3) * (R_i["H2"] - R_i["He"]) / r.mean_refractivity_m3
record(5, "fractional uncertainty of N is the composition term 0.0233 +- 0.0002 at every level, "
       "with the term lists",
       bool(np.all(np.abs(frac - 0.0233) <= 0.0002)) and cN.included == ("composition",)
       and cN.unstated == ("pressure", "temperature"),
       f"delta N / N over {frac.size} levels: min {frac.min():.10f}, max {frac.max():.10f}, "
       f"spread {frac.max() - frac.min():.2e}\n"
       f"closed form 0.03 (1 - x_NH3)(R_H2 - R_He) / R_bar: max departure "
       f"{float(np.max(np.abs(frac - closed_form))):.2e}; STP check 0.03 x (136 - 35) / 129.94 = "
       f"{0.03 * (136 - 35) / 129.94:.5f}\n"
       f"uncertainty_terms_included = {cN.included}\n"
       f"uncertainty_terms_unstated = {cN.unstated}\n"
       f"uncertainty_kind_conversions = {cN.conversions}\n"
       f"closure read from the declaration and checked against the values: rule "
       f"{r.closure['rule']!r}, share {r.closure['share']}, partner {r.closure['partner']}, "
       f"assigned first {r.closure['assigned_first']}\n"
       "inside the composition term: "
       f"{r.companions['mean_refractivity_uncertainty_m3'].included} entered, "
       f"{r.companions['mean_refractivity_uncertainty_m3'].unstated} left out (the per molecule "
       "uncertainties are NaN since the SPEC_01 v0.19 rebuild, and x_NH3 states none)")

# ---------------------------------------------------------------------------
# 6. The radius uncertainty, its two terms and their partials, cross checked
# ---------------------------------------------------------------------------
cr = r.companions["radius_uncertainty_m"]
P, Tm = r.partials, r.terms
quad = math.hypot(Tm["r0_anchor_radius_m"], Tm["r0_label_latitude_m"])


def rerun(thermo_attrs=None, anchor_radius_shift=0.0):
    thermo = inputs.thermo.copy()
    geodesy = inputs.geodesy.copy(deep=True)
    if thermo_attrs:
        thermo.attrs.update(thermo_attrs)
    if anchor_radius_shift:
        geodesy[manifest.anchor_quantity].values[:] += anchor_radius_shift
    return freeze_anchor(dataclasses.replace(inputs, thermo=thermo, geodesy=geodesy), manifest)


step_deg = 0.05
phi_g = float(inputs.thermo.attrs["latitude_planetographic_deg"])
north = rerun({"latitude_planetographic_deg": np.float64(phi_g + step_deg)})
south = rerun({"latitude_planetographic_deg": np.float64(phi_g - step_deg)})
fd_dphi = (north.phi_c_rad - south.phi_c_rad) / math.radians(2 * step_deg)
fd_dr0_dphi_g = (north.r0_m - south.r0_m) / math.radians(2 * step_deg)
up = rerun(anchor_radius_shift=1000.0)
down = rerun(anchor_radius_shift=-1000.0)
fd_dr0_dr = (up.r0_m - down.r0_m) / 2000.0
coupling = P["dr0_dphi_c_m_per_rad"] * P["dphi_c_dr_anchor_rad_per_m"]

cphi = r.companions["latitude_planetocentric_uncertainty_rad"]
range_text = "label_latitude: range taken as a uniform half width, divided by sqrt(3)"
ok6 = (bool(np.all(cr.value == cr.value[0])) and abs(cr.value[0] - quad) < 1e-6
       and cr.included == ("anchor_radius", "label_latitude")
       and cr.unstated == ("height", "anchor_height")
       and cr.conversions == (range_text,)
       and abs(fd_dphi - P["dphi_c_dphi_g"]) < 1e-3
       and abs(fd_dr0_dphi_g - P["dr0_dphi_c_m_per_rad"] * P["dphi_c_dphi_g"])
       < 1e-3 * abs(fd_dr0_dphi_g)
       and abs(fd_dr0_dr - P["dr0_dr_anchor"]) < 1e-4
       and abs(Tm["r0_anchor_radius_m"] / 1e3 - 12.37) < 0.005
       and abs(Tm["r0_label_latitude_m"] / 1e3 - 10.73) < 0.005
       and abs(cr.value[0] / 1e3 - 16.38) < 0.005
       and abs(math.degrees(float(cphi.value)) - 0.109) < 0.0005)
record(6, "radius_uncertainty_m is the same at every level and is the quadrature of the anchor "
       "term (12.37 km) and the label term (10.73 km), 16.38 km; phi_c carries 0.109 deg",
       ok6,
       f"radius_uncertainty_m: {cr.value[0]:.3f} m at all {cr.value.size} levels "
       f"({bool(np.all(cr.value == cr.value[0]))}); included {cr.included}, unstated {cr.unstated}\n"
       f"uncertainty_kind_conversions {cr.conversions}\n"
       f"height companions NaN at every level: "
       f"{bool(np.all(np.isnan(inputs.thermo['height_uncertainty_m'].values)))}\n"
       f"ANCHOR TERM  dr0/dr_anchor, total = {P['dr0_dr_anchor']:.6f}: latitude held "
       f"{P['dr0_dr_anchor_latitude_held']:.6f} (central difference, +- "
       f"{P['anchor_radius_step_m']:g} m) plus the coupling dr0/dphi_c x dphi_c/dr_anchor "
       f"{P['dr0_dr_anchor_latitude_coupling']:.6f}\n"
       f"             x {P['anchor_radius_uncertainty_1sigma_m']:g} m (declared "
       f"{P['declared_anchor_radius_uncertainty_m']:g} m {P['declared_anchor_radius_uncertainty_kind']}) "
       f"= {Tm['r0_anchor_radius_m']:.3f} m (12.37 km)\n"
       f"LABEL TERM   dr0/dphi_c = r0 G_phi / g = {P['dr0_dphi_c_m_per_rad'] / 1e3:.4f} km/rad "
       f"(march difference {P['dr0_dphi_c_difference_m_per_rad'] / 1e3:.4f})\n"
       f"             dphi_c/dphi_g = 1 / (1 + dpsi/dphi_c) = {P['dphi_c_dphi_g']:.6f}, with "
       f"dpsi/dphi_c along the surface {P['dpsi_dphi_c_along_surface']:.6f}\n"
       f"             x {math.degrees(P['label_latitude_uncertainty_1sigma_rad']):.6f} deg (declared "
       f"{math.degrees(P['declared_label_latitude_uncertainty_rad']):g} deg "
       f"{P['declared_label_latitude_uncertainty_kind']}, divided by sqrt(3)) = "
       f"{Tm['r0_label_latitude_m']:.3f} m (10.73 km)\n"
       f"quadrature {quad:.3f} m (16.38 km)\n"
       f"cross check by rerunning the whole anchor fixed point:\n"
       f"  label 36.3 +- {step_deg} deg: dphi_c/dphi_g = {fd_dphi:.6f} "
       f"(propagation {P['dphi_c_dphi_g']:.6f}); dr0/dphi_g = {fd_dr0_dphi_g / 1e3:.4f} km/rad "
       f"(propagation {P['dr0_dphi_c_m_per_rad'] * P['dphi_c_dphi_g'] / 1e3:.4f})\n"
       f"  anchor radius +- 1 km: total dr0/dr_anchor = {fd_dr0_dr:.6f}, propagation "
       f"{P['dr0_dr_anchor']:.6f}, difference {fd_dr0_dr - P['dr0_dr_anchor']:+.2e}\n"
       f"scalar companions: r0 {r.companions['anchor_isobar_radius_uncertainty_m'].value:.3f} m, "
       f"phi_c {math.degrees(float(cphi.value)):.6f} deg (0.109), "
       f"psi {math.degrees(r.companions['psi_uncertainty_rad'].value):.6f} deg; conversions on "
       f"phi_c {cphi.conversions}")

# ---------------------------------------------------------------------------
# 7. The inverse closes
# ---------------------------------------------------------------------------
T_back = red.temperature_from_refractivity(p, r.mean_refractivity_m3, r.refractivity)
rel = np.abs(T_back - r.temperature_K) / r.temperature_K
record(7, "T = p R_bar / (k_B N) returns the tabulated temperature to 1e-12 relative at every "
       "level",
       float(rel.max()) < 1e-12,
       f"max relative departure {float(rel.max()):.2e} over {rel.size} levels, at "
       f"{p[int(np.argmax(rel))]:g} Pa")

# ---------------------------------------------------------------------------
# 8. Beyond the specification: every companion's term lists, and the NaN rule
# ---------------------------------------------------------------------------
lines = []
for key, c in r.companions.items():
    v = np.atleast_1d(c.value)
    shown = "NaN" if np.all(np.isnan(v)) else f"{np.nanmin(v):.6g} to {np.nanmax(v):.6g}"
    lines.append(f"{key:42s} {shown:28s} included {c.included}, unstated {c.unstated}")
cn = r.companions["number_density_uncertainty_m3"]
total, inc, uns = red.quadrature({"a": np.array([3.0, np.nan]), "b": np.array([4.0, 4.0]),
                                  "c": np.array([np.nan, np.nan])})
ok8 = (bool(np.all(np.isnan(cn.value))) and cn.included == ()
       and cn.unstated == ("pressure", "temperature")
       and total.tolist() == [5.0, 4.0] and inc == ("a", "b") and uns == ("a", "c"))
record(8, "beyond the specification: the term lists of every companion, and unstated terms left "
       "out rather than zeroed",
       ok8,
       "\n".join(lines)
       + f"\nnumber_density_uncertainty_m3 has no stated term and is NaN: "
         f"{bool(np.all(np.isnan(cn.value)))}\n"
         f"quadrature({{a: [3, NaN], b: [4, 4], c: [NaN, NaN]}}) = {total.tolist()}, included "
         f"{inc}, unstated {uns}: a partly stated term is listed in both")

# ---------------------------------------------------------------------------
# 9. Beyond the specification: the closure is refused when the file does not show one
# ---------------------------------------------------------------------------
def closure_refusal(edit):
    tree = inputs.composition.copy()
    tree.dataset = edit(tree.to_dataset(inherit=False))
    try:
        rdc.composition_closure(tree)
        return False, "not refused"
    except rdc.ReductionError as exc:
        return True, str(exc)


def without_attr(ds, name):
    ds = ds.copy()
    ds.attrs.pop(name)
    return ds


def with_attr(ds, name, value):
    ds = ds.copy()
    ds.attrs[name] = value
    return ds


cases = [
    ("x_H2_uncertainty dropped", lambda ds: ds.drop_vars("x_H2_uncertainty")),
    ("closure_rule removed", lambda ds: without_attr(ds, "closure_rule")),
    ("closure_species declared as 'He H2'", lambda ds: with_attr(ds, "closure_species", "He H2")),
    ("closure_species declared as 'H2 NH3'", lambda ds: with_attr(ds, "closure_species", "H2 NH3")),
]
outcomes = [(label, *closure_refusal(edit)) for label, edit in cases]
record(9, "beyond the specification: a closure declaration that is missing or disagrees with "
       "the values is refused",
       all(refused for _, refused, _ in outcomes),
       "\n".join(f"{label}: refused {refused}\n  {message}" for label, refused, message in outcomes))

print()
failed = [x for x in results if not x[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
