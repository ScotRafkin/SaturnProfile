"""Acceptance for SPEC_17 v0.2 Step 1: what the Voyager-epoch observations require of the wind.

Part A: Lindal et al. (1985) Fig. 9's 100 mbar altitudes (`occul_data/lindal/isobar_100mbar.csv`,
written by `digitize_fig9.py`) against the delivered 100 mbar surface above the same reference geoid.
Part B: the linearized balance weight `lindal_wind.balance_weight`, the segment constraints from the
three Voyager profiles, and the linearization tested in a copy of the tree (`perturb.py`).

Reads the SPEC_15 runs as they are, built if absent; the perturbed runs are copies under
`reports/step17_1/perturbed/`, not registered runs. Measures only: no bounds. Writes under
`reports/step17_1/`.

Latitudes are given planetographic, with planetocentric in parentheses; Lindal's Fig. 9 and the
model's grids are planetocentric.

Run from the repository root.
"""

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import digitize_fig9  # noqa: E402
import perturb  # noqa: E402

from casspian.forward import transfer as forward_transfer  # noqa: E402
from casspian.lib import geoid as gd  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.lib.constants import BOLTZMANN_CONSTANT, MOLAR_GAS_CONSTANT  # noqa: E402
from casspian.lib.gravity import g_eff_radial, omega_abs  # noqa: E402
from casspian.tools.lindal import lindal_wind as lw  # noqa: E402
from casspian.tools.lindal import radio_profiles as rp  # noqa: E402
from casspian.tools.run import run_inputs  # noqa: E402

HERE = Path("reports/step17_1")
HERE.mkdir(parents=True, exist_ok=True)
ANCHOR_C = 30.80556842739218
#: Lindal's reference geoid (Fig. 9 caption): the polar radius it passes through.
REFERENCE_POLAR_RADIUS_m = 54438.0e3
HYPOTHESIS = ("hypothesis 1 (case ii)", "hypothesis 2 (case iii)")
#: Each Voyager target: its planetocentric latitude (SPEC_15), its runs (hypothesis 1, 2), its profile.
TARGETS = {
    "3.0 S": {"c": -2.3736176374048172, "runs": ("lindal_iris_v_sl_ii_110_03s", "lindal_iris_v_sl_iii_110_03s"),
              "observed": "occul_data/lindal/voyager1_egress.csv", "point": "voyager1", "epoch": "Voyager 1 egress"},
    "31.2 S": {"c": -26.06095353762054, "runs": ("lindal_iris_v_sl_ii_110", "lindal_iris_v_sl_iii_110"),
               "observed": "occul_data/lindal/voyager2_egress.csv", "point": "voyager2", "epoch": "Voyager 2 egress"},
    "74.2 S": {"c": -71.06804211320447, "runs": ("lindal_iris_v_sl_ii_110_74s", "lindal_iris_v_sl_iii_110_74s"),
               "observed": "occul_data/lindal/voyager1_ingress.csv", "point": "voyager1", "epoch": "Voyager 1 ingress"},
}
#: The 74.2 S runs carry the isobars from the anchor across all three segments; each hypothesis's w is
#: computed on its own run's (SPEC_17 v0.2 item 3).
PATH_RUNS = ("lindal_iris_v_sl_ii_110_74s", "lindal_iris_v_sl_iii_110_74s")
SEGMENTS = (("segment 1", ANCHOR_C, TARGETS["3.0 S"]["c"]),
            ("segment 2", TARGETS["3.0 S"]["c"], TARGETS["31.2 S"]["c"]),
            ("segment 3", TARGETS["31.2 S"]["c"], TARGETS["74.2 S"]["c"]))
TABLE_MBAR = (730.0, 290.0, 110.0, 30.0, 10.0, 4.5, 1.0)
#: SPEC_15: the equatorial profiles carry wave structure above 30 mbar. REVIEW_15_step1: Fig. 5's
#: drawn line is 1.5 to 3 K warm from 600 to 1200 mbar, which segment 3 carries.
WAVES_Pa = 3000.0
DEPTH_BIAS_Pa = (60000.0, 120000.0)
results, summary = [], {}


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def product(name):
    run = Path("forward") / name
    path = run / "output" / f"{name}_profile.nc"
    if not path.exists():
        run_inputs.build(run / f"{name}_build.toml")
        path = forward_transfer.run(run / f"{name}.toml").product
    return tree(path)


def run_inputs_of(name):
    run = Path("forward") / name / "inputs"
    return {kind: cio.read(run / f"{name}_{kind}.nc", kind) for kind in ("wind", "composition", "gravity", "rotation")}


def interp_ln_p(p, values, at):
    p, values = np.asarray(p, dtype="float64"), np.asarray(values, dtype="float64")
    order = np.argsort(p)
    return np.interp(np.log(at), np.log(p[order]), values[order], left=np.nan, right=np.nan)


def segment_integral(latitude_deg, values, start, end):
    """`int_start^end values dphi` (phi in radians) along the latitude grid, endpoints interpolated;
    oriented from `start` to `end`, so negative in extent when going south. `values` (..., latitude)."""
    lat = np.asarray(latitude_deg, dtype="float64")
    lo, hi = min(start, end), max(start, end)
    inside = (lat > lo) & (lat < hi)
    nodes = np.concatenate([[lo], lat[inside], [hi]])
    v = np.asarray(values, dtype="float64")
    at = np.stack([np.interp(nodes, lat, row) for row in v.reshape(-1, lat.size)])
    total = np.trapezoid(at, np.radians(nodes), axis=-1).reshape(v.shape[:-1])
    return total if end > start else -total


def temperature_on_wind_grid(prod, wind):
    """T on each isobar of the run, `(level, wind latitude)`, NaN off the run's path: the isobars'
    ln N with the target's mean refractivity (uniform in latitude, the run's composition) and the
    isobar's label."""
    root = prod.to_dataset(inherit=False)
    iso = prod["isobars"].to_dataset(inherit=False)
    lat_m = np.asarray(iso["latitude_planetocentric_deg"].values)
    labels = np.asarray(root["pressure_label_Pa"].values)
    r_bar = np.asarray(root["mean_refractivity_m3"].values)
    lat_w = np.asarray(wind["latitude_planetocentric_deg"].values)
    T = labels[:, None] * r_bar[:, None] / (BOLTZMANN_CONSTANT * np.exp(np.asarray(iso["ln_refractivity_lindal"].values)))
    # Held for one wind grid step beyond the path's ends, so that the anchor and the target, which lie
    # inside the step, are interpolated between two nodes.
    step = float(np.max(np.diff(lat_w)))
    off = (lat_w < lat_m.min() - step) | (lat_w > lat_m.max() + step)
    out = np.stack([np.interp(lat_w, lat_m, row) for row in T])
    out[:, off] = np.nan
    return labels, out


def full_balance(u_field, pressure_grid, inputs, levels, temperature):
    """`(d ln T / dphi)_p = -(2 Omega_abs r / g)(s D + (cos phi / r)(du/dphi)_p)`, the relation
    `construct` solves, evaluated whole for a wind field `u_field` `(latitude, pressure)`, at `levels`
    with `temperature` `(level, latitude)` held. Its change between two winds, less `w Delta s`, is the
    share of the wind's own change, which the linearization holds fixed."""
    wind, gravity, rotation = inputs["wind"], inputs["gravity"], inputs["rotation"]
    lat = np.asarray(wind["latitude_planetocentric_deg"].values, dtype="float64")
    phi = np.radians(lat)
    Omega = float(rotation["angular_rate_rad_s"])
    GM, R_norm = float(gravity["GM_m3s2"]), float(gravity["normalization_radius_m"])
    J, degrees = np.asarray(gravity["J"].values, dtype="float64"), np.asarray(gravity["degree"].values)
    m_bar = lw.mean_molar_mass(inputs["composition"], levels, lat)
    x = np.log(pressure_grid)
    order = np.argsort(x)
    shear_grid = np.gradient(u_field[:, order], x[order], axis=1)
    out = np.full((levels.size, lat.size), np.nan)
    for k, level in enumerate(levels):
        u_k = np.array([np.interp(np.log(level), x[order], row[order]) for row in u_field])
        s_k = np.array([np.interp(np.log(level), x[order], row) for row in shear_grid])
        r = np.asarray(gd.wind_geoid(phi, R_norm, "equatorial_radius", lambda y, u_k=u_k: np.interp(y, phi, u_k),
                                     Omega, GM, J, degrees, R_norm).radius, dtype="float64")
        g = g_eff_radial(u_k, r, phi, Omega, GM, J, degrees, R_norm)
        rate = omega_abs(u_k, r, phi, Omega)
        d = -g * m_bar[k] / (MOLAR_GAS_CONSTANT * temperature[k]) * (np.sin(phi) - np.cos(phi) / r * np.gradient(r, phi))
        out[k] = -(2.0 * rate * r / g) * (s_k * d + np.cos(phi) / r * np.gradient(u_k, phi))
    return out


products = {name: product(name) for target in TARGETS.values() for name in target["runs"]}

# ---------------------------------------------------------------------------
# Part A. The 100 mbar surface
# ---------------------------------------------------------------------------
fig9 = digitize_fig9.digitize(digitize_fig9.render())
table = np.genfromtxt([line for line in Path(digitize_fig9.OUT).read_text(encoding="utf-8").splitlines()
                       if not line.startswith("#")], delimiter=",", names=True, dtype=None, encoding="utf-8")
assert np.allclose(table["latitude_planetocentric_deg"], np.round(fig9.latitude_deg, 2)), "the data file is not the digitizer's"
assert np.allclose(table["altitude_km"], np.round(fig9.altitude_km, 1)), "the data file is not the digitizer's"
base_inputs = run_inputs_of(perturb.BASE)
gravity, rotation = base_inputs["gravity"], base_inputs["rotation"]
Omega = float(rotation["angular_rate_rad_s"])
GM, R_norm = float(gravity["GM_m3s2"]), float(gravity["normalization_radius_m"])
J, degrees = np.asarray(gravity["J"].values, dtype="float64"), np.asarray(gravity["degree"].values)
lat_ref = np.arange(-90.0, 90.0001, 0.05)
r_ref, _, residual_U = gd.reference_geoid(np.radians(lat_ref), REFERENCE_POLAR_RADIUS_m, Omega, GM, J, degrees, R_norm)
reference = lambda lat_deg: np.interp(lat_deg, lat_ref, r_ref)

anchor_group = products[TARGETS["31.2 S"]["runs"][0]]["anchors/lindal"].to_dataset(inherit=False)
assert float(anchor_group["anchor_isobar_pressure_Pa"]) == 1.0e4
anchor_alt_km = (float(anchor_group["anchor_isobar_radius_m"]) - reference(ANCHOR_C)) / 1e3


def surface(name):
    """The run's 100 mbar surface along its path: latitude, radius (the gauge isobar, Eq. B3)."""
    rs = products[name]["reference_surface"].to_dataset(inherit=False)
    return np.asarray(rs["latitude_planetocentric_deg"].values), np.asarray(rs["reference_surface_radius_m"].values)


lines_a = [f"Fig. 9 calibration: {fig9.altitude_ticks} altitude ticks, RMS {fig9.altitude_rms_km:.2f} km; "
           f"{fig9.latitude_ticks} latitude ticks, RMS {fig9.latitude_rms_deg:.2f} deg (the ticks are drafted unevenly)",
           f"reference geoid: lib.geoid.reference_geoid through {REFERENCE_POLAR_RADIUS_m / 1e3:.0f} km at the pole, the run's "
           f"gravity (GM {GM:.6e}, J {', '.join(f'{j:.4e}' for j in J)}) and rotation ({float(rotation['period_s']):.1f} s); "
           f"equatorial radius {reference(0.0) / 1e3:.1f} km, largest |U - U_ref| {float(np.max(residual_U)):.2e} m2/s2",
           "delivered minus observed 100 mbar altitude above the reference geoid, km, at the comparison latitude:",
           "  point              digitized  compared  difference  observed   hyp 1    hyp 2    hyp1-obs  hyp2-obs  slope (km/deg)"]
rows_a = []
for lat_obs, alt_obs, obs in zip(table["latitude_planetocentric_deg"], table["altitude_km"], table["observation"]):
    if obs == "pioneer11":
        lines_a.append(f"  pioneer11          {lat_obs:+8.2f}  reference only, not compared (two years from the anchor); "
                       f"observed {alt_obs:.1f} km")
        continue
    if abs(lat_obs - ANCHOR_C) < 2.0:
        where, at = "anchor 36.3 N", ANCHOR_C
        delivered = (anchor_alt_km, anchor_alt_km)
        slope = np.nan
    else:
        where = next(k for k, t in TARGETS.items() if abs(t["c"] - lat_obs) < 2.0)
        spec = TARGETS[where]
        delivered, slope = [], np.nan
        for name in spec["runs"]:
            lat_s, r_s = surface(name)
            covered = lat_s.min() <= lat_obs <= lat_s.max()
            at = lat_obs if covered else spec["c"]
            delivered.append((np.interp(at, lat_s, r_s) - reference(at)) / 1e3)
        # The slope of the delivered altitude there, from hypothesis 1's 74.2 S path, for the latitude difference.
        lat_p, r_p = surface(PATH_RUNS[0])
        alt_p = (r_p - reference(lat_p)) / 1e3
        slope = float(np.interp(at, lat_p, np.gradient(alt_p, lat_p)))
        where = f"{where} {spec['epoch'].split()[-1]}"
    d1, d2 = delivered[0] - alt_obs, delivered[1] - alt_obs
    rows_a.append({"point": where, "digitized_deg": float(lat_obs), "compared_deg": float(at), "observed_km": float(alt_obs),
                   "delivered_km": [float(delivered[0]), float(delivered[1])], "slope_km_per_deg": slope})
    lines_a.append(f"  {where:18s} {lat_obs:+8.2f}  {at:+8.2f}  {lat_obs - at:+8.2f}    {alt_obs:7.1f}  {delivered[0]:7.1f}  "
                   f"{delivered[1]:7.1f}  {d1:+8.1f}  {d2:+8.1f}  {slope:+8.2f}")
summary["isobar_100mbar"] = rows_a
lines_a.append("the comparison is at the point's own latitude where the run's path carries the surface there, else at "
               "the target's; the slope times the difference is the part of the difference that is the plotted latitude")

fig, ax = plt.subplots(figsize=(8, 9))
ax.plot(fig9.curve_altitude_km, fig9.curve_latitude_deg, ".", color="0.55", markersize=1.2,
        label="Lindal's best-fit geoid (stippled curve, traced)")
for k, name in enumerate(PATH_RUNS):
    lat_p, r_p = surface(name)
    ax.plot((r_p - reference(lat_p)) / 1e3, lat_p, color=("C2", "C4")[k], linewidth=1.3,
            label=f"delivered 100 mbar surface, {HYPOTHESIS[k]} (the 74.2 S run's path)")
for lat_obs, alt_obs, obs in zip(table["latitude_planetocentric_deg"], table["altitude_km"], table["observation"]):
    if obs == "pioneer11":
        ax.plot(alt_obs, lat_obs, "o", mfc="none", mec="0.3", markersize=8)
        ax.annotate("Pioneer 11 (reference only)", (alt_obs, lat_obs), xytext=(6, -12), textcoords="offset points", fontsize=8)
    else:
        ax.plot(alt_obs, lat_obs, "o", color="C3", markersize=7)
ax.plot([], [], "o", color="C3", label="Lindal Fig. 9, Voyager 100 mbar altitudes")
for row in rows_a:
    for k in range(2):
        ax.plot(row["delivered_km"][k], row["compared_deg"], ("s", "D")[k], mfc="none", mec=("C2", "C4")[k], markersize=8)
ax.axvline(0, color="0.6", linewidth=0.8)
ax.set_xlabel("altitude of the 100 mbar surface above the reference geoid (km)")
ax.set_ylabel("planetocentric latitude (deg)")
ax.set_ylim(-90, 90)
ax.grid(alpha=0.3)
ax.legend(fontsize=7, loc="lower right")
ax.set_title("The 100 mbar surface above Lindal's reference geoid (polar radius 54,438 km, no wind)", fontsize=9)
fig.tight_layout()
fig.savefig(HERE / "isobar_100mbar.png", dpi=120)
plt.close(fig)

# ---------------------------------------------------------------------------
# Part B. The weight and the segment constraints
# ---------------------------------------------------------------------------
weights = {}
for name in PATH_RUNS:
    inputs = run_inputs_of(name)
    labels, T_w = temperature_on_wind_grid(products[name], inputs["wind"])
    w = lw.balance_weight(inputs["wind"], inputs["composition"], inputs["gravity"], inputs["rotation"], labels, T_w)
    weights[name] = {"labels": labels, "T": T_w, "w": w, "inputs": inputs,
                     "lat": np.asarray(inputs["wind"]["latitude_planetocentric_deg"].values)}

misfit = {}
for where, spec in TARGETS.items():
    observed = rp.read_csv(spec["observed"])
    for k, name in enumerate(spec["runs"]):
        root = products[name].to_dataset(inherit=False)
        p, T = np.asarray(root["pressure_Pa"].values), np.asarray(root["temperature_K"].values)
        T_obs = interp_ln_p(np.asarray(observed["pressure_mbar"]) * 100.0, observed["temperature_K"], p)
        misfit[(where, k)] = np.log(T_obs / T)

constraints = {}
for k, name in enumerate(PATH_RUNS):
    W = weights[name]
    required = {"segment 1": misfit[("3.0 S", k)],
                "segment 2": misfit[("31.2 S", k)] - misfit[("3.0 S", k)],
                "segment 3": misfit[("74.2 S", k)] - misfit[("31.2 S", k)]}
    for segment, start, end in SEGMENTS:
        integral_w = segment_integral(W["lat"], W["w"], start, end)
        constraints[(segment, k)] = {"required": required[segment], "integral_w": integral_w,
                                     "mean_shear": required[segment] / integral_w}

table_b = ["segment-mean Delta s (m/s per unit ln p) = required integral of w Delta s / integral of w; at "
           + ", ".join(f"{p:g}" for p in TABLE_MBAR) + " mbar:"]
summary["segments"] = {}
for segment, start, end in SEGMENTS:
    table_b.append(f"{segment}, {start:+.2f} to {end:+.2f} deg planetocentric:")
    for k, name in enumerate(PATH_RUNS):
        c = constraints[(segment, k)]
        labels = weights[name]["labels"]
        at = lambda v: interp_ln_p(labels, v, np.array(TABLE_MBAR) * 100.0)
        table_b.append(f"  {HYPOTHESIS[k]:24s} mean Delta s " + " ".join(f"{v:+8.2f}" for v in at(c["mean_shear"])))
        table_b.append(f"  {'':24s} required     " + " ".join(f"{v:+8.4f}" for v in at(c["required"])))
        table_b.append(f"  {'':24s} int w dphi   " + " ".join(f"{v:+8.4f}" for v in at(c["integral_w"])))
        summary["segments"][f"{segment}|{k}"] = {"mean_shear": at(c["mean_shear"]).tolist(), "required": at(c["required"]).tolist(),
                                                 "integral_w": at(c["integral_w"]).tolist()}
table_b.append("required is ln(T_obs / T_delivered) differenced between the segment's ends; int w dphi is oriented from the "
               "segment's north end to its south end, in s/m. Segments 1 and 2 above 30 mbar carry the equatorial profiles' "
               "wave structure (SPEC_15); segment 3 from 600 to 1200 mbar carries Fig. 5's depth bias (REVIEW_15_step1), not "
               "corrected.")

fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
for k, name in enumerate(PATH_RUNS):
    W = weights[name]
    lim = np.nanpercentile(np.abs(W["w"]), 99)
    mesh = axes[k].pcolormesh(W["lat"], W["labels"] / 100.0, W["w"], cmap="RdBu_r", vmin=-lim, vmax=lim, shading="nearest")
    for _, start, end in SEGMENTS:
        axes[k].axvline(start, color="0.3", linewidth=0.7, linestyle=":")
    axes[k].axvline(SEGMENTS[-1][2], color="0.3", linewidth=0.7, linestyle=":")
    axes[k].set_yscale("log")
    axes[k].set_xlim(-75, 35)
    axes[k].set_xlabel("planetocentric latitude (deg)")
    axes[k].set_title(f"w, {HYPOTHESIS[k]} ({name})", fontsize=9)
    fig.colorbar(mesh, ax=axes[k], label="w (s m-1 rad-1)")
axes[0].invert_yaxis()
axes[0].set_ylabel("pressure (mbar)")
fig.suptitle("The linearized balance weight: Delta(d ln T / dphi) = w Delta s, on each run's isobars along its path "
             "(dotted: the targets)", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "weight.png", dpi=110)
plt.close(fig)

fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
for j, (segment, start, end) in enumerate(SEGMENTS):
    for k, name in enumerate(PATH_RUNS):
        axes[j].plot(constraints[(segment, k)]["mean_shear"], weights[name]["labels"] / 100.0, color=("C2", "C4")[k],
                     linewidth=1.3, label=HYPOTHESIS[k])
    if segment != "segment 3":
        axes[j].axhspan(0.1, WAVES_Pa / 100.0, color="C1", alpha=0.12, label="wave structure in the equatorial profile (SPEC_15)")
    else:
        axes[j].axhspan(DEPTH_BIAS_Pa[0] / 100.0, DEPTH_BIAS_Pa[1] / 100.0, color="C1", alpha=0.12,
                        label="Fig. 5 depth bias, 1.5 to 3 K warm (REVIEW_15_step1)")
    axes[j].axvline(0, color="0.5", linewidth=0.8)
    axes[j].set_yscale("log")
    axes[j].set_xlim(-60, 60)
    axes[j].set_title(f"{segment}: {start:+.2f} to {end:+.2f} deg planetocentric", fontsize=9)
    axes[j].set_xlabel("segment-mean Delta s (m/s per unit ln p)")
    axes[j].grid(alpha=0.3, which="both")
    axes[j].legend(fontsize=7)
axes[0].set_ylim(1300, 0.3)
axes[0].set_ylabel("pressure (mbar)")
fig.suptitle("The shear change each segment asks for, on its mean weight (no fit)", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "segments.png", dpi=110)
plt.close(fig)

# ---------------------------------------------------------------------------
# 1. The linearization, tested
# ---------------------------------------------------------------------------
W = weights[PATH_RUNS[0]]
base_wind = W["inputs"]["wind"]
p_grid = np.asarray(base_wind["pressure_Pa"].values, dtype="float64")
u_base = np.asarray(base_wind["u_total_ms"].values, dtype="float64")
predicted_unit = {where: segment_integral(W["lat"], W["w"] * perturb.delta_s(W["lat"], W["labels"], 1.0), ANCHOR_C, spec["c"])
                  for where, spec in TARGETS.items()}
band = (W["labels"] >= perturb.BAND_Pa[0]) & (W["labels"] <= perturb.BAND_Pa[1])
T_target = products[TARGETS["31.2 S"]["runs"][0]].to_dataset(inherit=False)["temperature_K"].values
# c moves the 31.2 S temperature by about 2 K in the band; rounded to 0.5 m/s per unit ln p.
c = float(np.round(2.0 * 2.0 / np.nanmean(np.abs(predicted_unit["31.2 S"][band] * T_target[band]))) / 2.0)
amplitudes = {"c": c, "2c": 2.0 * c}

delivered_change, full_change = {}, {}
# Below the band's lower taper Delta u is zero, so the whole balance's change is zero there exactly.
moved = W["labels"] < perturb.BAND_Pa[1] * np.exp(perturb.TAPER_LN_P) * 1.05
whole_base = full_balance(u_base, p_grid, W["inputs"], W["labels"][moved], W["T"][moved])
for label, amp in amplitudes.items():
    u_pert = u_base + perturb.delta_u(W["lat"], p_grid, amp)
    whole = np.zeros_like(W["w"])
    whole[moved] = full_balance(u_pert, p_grid, W["inputs"], W["labels"][moved], W["T"][moved]) - whole_base
    for where, spec in TARGETS.items():
        registered = spec["runs"][0]
        pert = tree(perturb.run_copy(registered, label, amp)).to_dataset(inherit=False)
        base = products[registered].to_dataset(inherit=False)
        delivered_change[(where, label)] = np.log(np.asarray(pert["temperature_K"].values) / np.asarray(base["temperature_K"].values))
        full_change[(where, label)] = segment_integral(W["lat"], whole, ANCHOR_C, spec["c"])

lines_1 = [f"Delta s = c A(phi) B(p): segment 2, tapered over 2 deg inside each end; 10 to 30 mbar, tapered over 0.5 "
           f"in ln p outside each end; c = {c:g} m/s per unit ln p, chosen to move 31.2 S by about 2 K; base "
           f"{perturb.BASE} (hypothesis 1); copies under {perturb.WORK}",
           "in the band (10 to 30 mbar, B = 1), delivered change in ln T over predicted (w Delta s), mean [range]; "
           "the whole balance's prediction (the wind's own change included) over w Delta s; linearity, Delta(2c) / (2 Delta(c)):"]
summary["linearization"] = {"c": c}
for where in TARGETS:
    pred = predicted_unit[where]
    for label, amp in amplitudes.items():
        d, f = delivered_change[(where, label)], full_change[(where, label)]
        if np.nanmax(np.abs(pred[band])) * amp < 1e-6:
            lines_1.append(f"  {where} {label:2s}: predicted zero (outside segment 2's path); delivered largest "
                           f"|Delta ln T| {np.nanmax(np.abs(d)):.2e}, whole balance {np.nanmax(np.abs(f)):.2e}")
            continue
        ratio = d[band] / (amp * pred[band])
        whole_ratio = f[band] / (amp * pred[band])
        lines_1.append(f"  {where} {label:2s}: delivered / predicted {np.mean(ratio):.4f} [{ratio.min():.4f}, {ratio.max():.4f}]; "
                       f"whole / w Delta s {np.mean(whole_ratio):.4f}; delivered / whole {np.mean(d[band] / f[band]):.4f}; "
                       f"predicted {np.mean(amp * pred[band]) * 100:+.3f} % ({np.mean(amp * pred[band] * T_target[band]):+.2f} K at 31.2 S's T)")
        summary["linearization"][f"{where}|{label}"] = {"ratio_mean": float(np.mean(ratio)), "ratio_min": float(ratio.min()),
                                                        "ratio_max": float(ratio.max()), "whole_over_linear": float(np.mean(whole_ratio))}
    if np.nanmax(np.abs(pred[band])) > 1e-6:
        lin = delivered_change[(where, "2c")][band] / (2.0 * delivered_change[(where, "c")][band])
        lines_1.append(f"  {where} linearity: {np.mean(lin):.4f} [{lin.min():.4f}, {lin.max():.4f}]")
        above = W["labels"] < perturb.BAND_Pa[0] * np.exp(-perturb.TAPER_LN_P)
        if above.any():
            lines_1.append(f"  {where} above the band (Delta s = 0, Delta u uniform): delivered Delta ln T at c {np.nanmean(delivered_change[(where, 'c')][above]):+.2e} "
                           f"mean, whole balance {np.nanmean(full_change[(where, 'c')][above]):+.2e}")

fig, axes = plt.subplots(1, 3, figsize=(16, 7), sharey=True)
for j, where in enumerate(TARGETS):
    for label, amp in amplitudes.items():
        style = "-" if label == "c" else "--"
        axes[j].plot(100 * delivered_change[(where, label)], W["labels"] / 100.0, style, color="C3", label=f"delivered, {label}")
        axes[j].plot(100 * amp * predicted_unit[where], W["labels"] / 100.0, style, color="C0", label=f"w Delta s, {label}")
        axes[j].plot(100 * full_change[(where, label)], W["labels"] / 100.0, style, color="C2", label=f"whole balance, {label}")
    axes[j].axhspan(perturb.BAND_Pa[0] / 100.0, perturb.BAND_Pa[1] / 100.0, color="0.85", alpha=0.5)
    axes[j].set_yscale("log")
    axes[j].set_title(where, fontsize=9)
    axes[j].set_xlabel("change in ln T x 100 (%)")
    axes[j].grid(alpha=0.3, which="both")
    axes[j].legend(fontsize=7)
axes[0].set_ylim(1300, 0.3)
axes[0].set_ylabel("pressure (mbar)")
fig.suptitle(f"The linearization: a known shear over segment 2 and 10 to 30 mbar (shaded), c = {c:g} m/s per unit ln p", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "linearization.png", dpi=110)
plt.close(fig)

record(1, "the linearization: delivered change over predicted in the perturbed band, at c and 2c", True, "\n".join(lines_1)
       + f"\nfigure {HERE / 'linearization.png'}")
record(2, "for the author: the 100 mbar surface, the weight, the segment constraints", True,
       "\n".join(lines_a + table_b) + f"\nfigures {HERE / 'isobar_100mbar.png'}, {HERE / 'weight.png'}, {HERE / 'segments.png'}")

(HERE / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8", newline="\n")
sys.exit(1 if failed else 0)
