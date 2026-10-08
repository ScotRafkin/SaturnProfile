"""Acceptance for SPEC_12 v0.2 Step 1: the Voyager 2 egress profile, and the first comparison.

Digitizes Lindal et al. (1985) Fig. 4 with `casspian.tools.lindal.egress_profile`, converts the
ingress (solid) curve the same way as the egress for the check against Table I, and rebuilds and
reruns `forward/lindal_iris_ii_110` at its target, Lindal's egress label converted as the anchor's
is (SPEC_12 section 1 deliverable 4). Writes under `reports/step12_1/`; the committed
`occul_data/lindal/voyager2_egress.csv` is read, never written. The run's inputs and product are
written in the run directory, as a run's are.

**The relaxation.** The tool is uncommitted, so the run's inputs carry `-dirty`; the one refusal
point, `casspian.lib.control._refuse_dirty_commit`, is replaced for this script's run only.

Run from the repository root.
"""

import json
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from PIL import Image, ImageDraw

from casspian.forward import transfer as forward_transfer
from casspian.lib import control as ctl
from casspian.lib.constants import AVOGADRO_CONSTANT, MOLAR_GAS_CONSTANT
from casspian.tools.lindal import egress_profile as ep
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.run import run_inputs

HERE = Path("reports/step12_1")
HERE.mkdir(parents=True, exist_ok=True)
PDF = Path("docs/Lindal_et_al_1985_AJ90_1136.pdf")
COMMITTED = Path("occul_data/lindal/voyager2_egress.csv")
TABLE = np.genfromtxt("occul_data/lindal/raw/lindal_table1.csv", delimiter=",", names=True)
RUN = Path("forward/lindal_iris_ii_110")
EGRESS_AMU = 2.135
LEVELS_mbar = (110.0, 290.0, 730.0)
IRIS = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
#: Lindal's labels, planetographic; the planetocentric latitudes are the tool's.
INGRESS_G, EGRESS_G = 36.3, -31.2
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def latitude(planetographic, planetocentric):
    """A latitude as the figures give it: planetographic, with planetocentric in parentheses."""
    side = lambda v: "N" if v >= 0 else "S"
    return (f"{abs(planetographic):.1f} {side(planetographic)} "
            f"({abs(planetocentric):.2f} {side(planetocentric)} planetocentric)")


INGRESS_LABEL = latitude(INGRESS_G, ep.INGRESS_LATITUDE_PLANETOCENTRIC_DEG)
EGRESS_LABEL = latitude(EGRESS_G, ep.EGRESS_LATITUDE_PLANETOCENTRIC_DEG)


def at_pressure(p, values, targets):
    order = np.argsort(np.log(p))
    return np.interp(np.log(targets), np.log(p)[order], np.asarray(values, dtype="float64")[order])


warnings.simplefilter("ignore")
result = ep.digitize(PDF)
inputs = ep.lindal_inputs()

# ---------------------------------------------------------------------------
# 0. The CSV is what the tool writes (beyond the specification)
# ---------------------------------------------------------------------------
egress = ep.to_pressure(result.dashed.altitude_km, result.dashed.temperature_K,
                        ep.EGRESS_LATITUDE_PLANETOCENTRIC_DEG,
                        sigma_K=ep.sigma_on(result.dashed, result.sigma), **inputs)
again = ep.write_csv(egress, HERE / "voyager2_egress.csv")
same = again.read_bytes() == COMMITTED.read_bytes()
zones = "\n".join(f"  zone {z['zone']['top_km']:g} to {z['zone']['bottom_km']:g} km by {z['zone']['scan']}: "
                  f"pieces (first, second) {z['segments']}, line width {z['line_width_px']:.0f} px, "
                  f"dashed the {z['dashed']} track, pieces agree {z['pieces_agree']}"
                  for z in result.zone_absence)
record(0, "running the tool reproduces the committed CSV byte for byte; the tracing's record", same,
       f"identical {same}; {egress.altitude_km.size} levels, {egress.altitude_km.min():g} to "
       f"{egress.altitude_km.max():g} km, {egress.pressure_Pa.max() / 100:.5g} to "
       f"{egress.pressure_Pa.min() / 100:.4g} mbar\n"
       f"altitude ticks: largest departure from one straight line {np.abs(result.altitude_residual_km).max():.2f} km "
       f"(the scale is piecewise between ticks); temperature ticks: largest residual "
       f"{np.abs(result.temperature_residual_K).max():.3f} K\n"
       f"bins filled across a gap or a spike: solid {int(result.solid.filled.sum())} of "
       f"{result.solid.altitude_km.size}, dashed {int(result.dashed.filled.sum())} of "
       f"{result.dashed.altitude_km.size}, sigma_T {int(result.sigma.filled.sum())} of "
       f"{result.sigma.altitude_km.size}\n" + zones + "\n"
       f"egress at {EGRESS_LABEL}: 1 bar radius "
       f"{egress.radius_1bar_m / 1e3:.2f} km, g {egress.gravity_ms2.min():.3f} to {egress.gravity_ms2.max():.3f} m/s2, "
       f"mean molar mass {egress.mean_molar_mass_kg_mol.min() * 1e3:.4f} to "
       f"{egress.mean_molar_mass_kg_mol.max() * 1e3:.4f} amu")

# ---------------------------------------------------------------------------
# 1. The ingress check against Table I
# ---------------------------------------------------------------------------
ingress = ep.to_pressure(result.solid.altitude_km, result.solid.temperature_K,
                         ep.INGRESS_LATITUDE_PLANETOCENTRIC_DEG, **inputs)
h_t, T_t, p_t = TABLE["altitude_km"], TABLE["temperature_K"], TABLE["pressure_mbar"] * 100.0
inside = (h_t >= ingress.altitude_km.min()) & (h_t <= ingress.altitude_km.max())
dT = np.interp(h_t[inside], ingress.altitude_km, ingress.temperature_K) - T_t[inside]
p_at = np.exp(np.interp(h_t[inside], ingress.altitude_km, np.log(ingress.pressure_Pa)))
dlnp = np.log(p_at / p_t[inside])
worst_T, worst_p = int(np.argmax(np.abs(dT))), int(np.argmax(np.abs(dlnp)))
fig, axes = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
axes[0].plot(dT, h_t[inside], ".-", color="C0")
axes[1].plot(100 * dlnp, h_t[inside], ".-", color="C3")
axes[0].set_xlabel("traced - Table I temperature (K)")
axes[1].set_xlabel("traced - Table I pressure (percent, as ln p)")
axes[0].set_ylabel("altitude above 1 bar (km)")
for ax in axes:
    ax.axvline(0, color="0.6", linewidth=0.8)
    ax.grid(alpha=0.3)
fig.suptitle(f"The ingress check: Fig. 4's solid curve at {INGRESS_LABEL}, traced and converted, "
             "against Lindal's Table I", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "ingress_check.png", dpi=120)
plt.close(fig)
record(1, f"the ingress check: the solid curve, traced and converted the same way at {INGRESS_LABEL}, against "
          "Lindal's Table I at the table's altitudes (reported)", True,
       f"{int(inside.sum())} Table I levels, {h_t[inside].min():g} to {h_t[inside].max():g} km\n"
       f"temperature: RMS {np.sqrt(np.mean(dT ** 2)):.2f} K, largest {dT[worst_T]:+.2f} K at "
       f"{h_t[inside][worst_T]:g} km\n"
       f"pressure: RMS {100 * np.sqrt(np.mean(dlnp ** 2)):.2f} percent, largest {100 * dlnp[worst_p]:+.2f} "
       f"percent at {h_t[inside][worst_p]:g} km ({p_t[inside][worst_p] / 100:g} mbar)\n"
       f"ingress at {INGRESS_LABEL}: 1 bar radius {ingress.radius_1bar_m / 1e3:.2f} km; "
       f"mean molar mass at 1 bar {float(np.interp(0.0, ingress.altitude_km, ingress.mean_molar_mass_kg_mol)) * 1e3:.4f} amu "
       f"(Lindal 2.135)\nfigure {HERE / 'ingress_check.png'}")

# ---------------------------------------------------------------------------
# 2. The overlay
# ---------------------------------------------------------------------------
sc, t_fit, s_fit = result.altitude_scale, result.temperature_fit, result.sigma_fit
row = sc.row_of
col = lambda T: (np.asarray(T) - t_fit[1]) / t_fit[0]
scol = lambda s: (np.asarray(s) - s_fit[1]) / s_fit[0]
image = Image.fromarray(np.where(result.ink, 0, 255).astype(np.uint8)).convert("RGB")
draw = ImageDraw.Draw(image)
for trace, colour, to_x in ((result.solid, (0, 160, 0), col), (result.dashed, (220, 0, 0), col),
                            (result.sigma, (0, 0, 220), scol)):
    for h, value, filled in zip(trace.altitude_km, trace.temperature_K, trace.filled):
        x, y = to_x(value), row(h)
        draw.ellipse([x - 4, y - 4, x + 4, y + 4], fill=(255, 170, 0) if filled else colour)
overlay = HERE / "overlay.png"
image.crop(ep.FIGURE_BOX).save(overlay)
record(2, "the traced curves over the scan, for the author to view", overlay.exists(),
       f"{overlay}: green the solid curve, {INGRESS_LABEL}; red the dashed, {EGRESS_LABEL}; blue sigma_T, each 1 km "
       f"bin as a disk; orange a bin filled across a gap or a spike")

# ---------------------------------------------------------------------------
# 3. The comparison: the run rebuilt and rerun at its target
# ---------------------------------------------------------------------------
_original = ctl._refuse_dirty_commit
ctl._refuse_dirty_commit = lambda path, attrs, consumer, what=None: str(attrs.get("casspian_git_commit", ""))
run_inputs.build(RUN / f"{RUN.name}_build.toml")
state = forward_transfer.run(RUN / f"{RUN.name}.toml")
ctl._refuse_dirty_commit = _original
product = xr.open_datatree(state.product, engine="netcdf4").load()
product.close()
root = product.to_dataset(inherit=False)
anchor = product["anchors/lindal"].to_dataset(inherit=False)
thermo = product["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
p = np.asarray(root["pressure_Pa"].values)
T = np.asarray(root["temperature_K"].values)
rho = np.asarray(root["number_density_m3"].values) * np.asarray(root["mean_molar_mass_kg_mol"].values) / AVOGADRO_CONSTANT
p_a = np.asarray(thermo["pressure_Pa"].values)
T_a = np.asarray(thermo["temperature_K"].values)
rho_a = (np.asarray(anchor["number_density_m3"].values) * np.asarray(anchor["mean_molar_mass_kg_mol"].values)
         / AVOGADRO_CONSTANT)
p_e, T_e, s_e = egress.pressure_Pa, egress.temperature_K, egress.temperature_uncertainty_K
rho_e = p_e * EGRESS_AMU * 1e-3 / (MOLAR_GAS_CONSTANT * T_e)
identity_K = T * np.asarray(root["pressure_identity_residual"].values)
passes = int(state.state.record["passes"])
T_i = ingress.temperature_K
lines, change_lines, rows, iris_at = [], [], {}, {}
for level in LEVELS_mbar:
    pa = level * 100.0
    delivered, observed = float(at_pressure(p, T, [pa])[0]), float(at_pressure(p_e, T_e, [pa])[0])
    sigma = float(at_pressure(p_e, s_e, [pa])[0])
    anchor_T = float(at_pressure(p_a, T_a, [pa])[0])
    traced_in = float(at_pressure(ingress.pressure_Pa, T_i, [pa])[0])
    table_in = float(at_pressure(p_t, T_t, [pa])[0])
    d_rho = float(at_pressure(p, rho, [pa])[0] / at_pressure(p_e, rho_e, [pa])[0] - 1.0)
    use = IRIS["pressure_mbar"] == level
    f = iris.fit(IRIS["latitude_planetographic_deg"][use], IRIS["temperature_K"][use],
                 grid=np.array([INGRESS_G, EGRESS_G]))
    iris_n, iris_s = float(f.value_K[0]), float(f.value_K[1])
    iris_at[level] = (iris_n, iris_s)
    radio_change, transfer_change = observed - traced_in, delivered - anchor_T
    rows[level] = {"delivered_K": delivered, "egress_K": observed, "sigma_K": sigma, "anchor_K": anchor_T,
                   "density_relative": d_rho, "iris_ingress_K": iris_n, "iris_egress_K": iris_s,
                   "traced_ingress_K": traced_in, "table_ingress_K": table_in,
                   "radio_change_K": radio_change, "transfer_change_K": transfer_change,
                   "iris_change_K": iris_s - iris_n}
    lines.append(f"{level:g} mbar: delivered {delivered:.2f} K, egress {observed:.2f} K (sigma_T {sigma:.2f}), "
                 f"delivered - egress {delivered - observed:+.2f} K; anchor {anchor_T:.2f} K, anchor - egress "
                 f"{anchor_T - observed:+.2f} K; IRIS fit {iris_n:.2f} K at {INGRESS_G:g} N and {iris_s:.2f} K at "
                 f"{abs(EGRESS_G):g} S planetographic, IRIS at {abs(EGRESS_G):g} S - egress "
                 f"{iris_s - observed:+.2f} K; density delivered / egress - 1 {100 * d_rho:+.2f} percent")
    change_lines.append(f"{level:g} mbar: Fig. 4 egress - ingress {radio_change:+.2f} K; transfer delivered - "
                        f"anchor {transfer_change:+.2f} K; difference {transfer_change - radio_change:+.2f} K; "
                        f"IRIS fit {abs(EGRESS_G):g} S - {INGRESS_G:g} N {iris_s - iris_n:+.2f} K; traced "
                        f"ingress - Table I {traced_in - table_in:+.2f} K")
levels_Pa = np.array(LEVELS_mbar) * 100.0
iris_n_all = np.array([iris_at[v][0] for v in LEVELS_mbar])
iris_s_all = np.array([iris_at[v][1] for v in LEVELS_mbar])
common = (p_e <= min(p.max(), p_a.max())) & (p_e >= max(p.min(), p_a.min()))
fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
axes[0].fill_betweenx(p_e / 100, T_e - s_e, T_e + s_e, color="C3", alpha=0.18, linewidth=0,
                      label="egress +- sigma_T")
axes[0].plot(T_e, p_e / 100, color="C3", linewidth=1.0, label=f"Voyager 2 egress (Lindal Fig. 4), {EGRESS_LABEL}")
axes[0].plot(T, p / 100, color="k", linewidth=1.3, label="delivered, lindal_iris_ii_110")
axes[0].plot(T_a, p_a / 100, color="0.55", linewidth=0.9, linestyle="--", label=f"anchor, {INGRESS_LABEL}")
axes[0].plot(iris_n_all, levels_Pa / 100, "o", markerfacecolor="none", markeredgecolor="C0", markersize=8,
             label=f"IRIS fit, {INGRESS_LABEL}")
axes[0].plot(iris_s_all, levels_Pa / 100, "s", color="C0", markersize=7, label=f"IRIS fit, {EGRESS_LABEL}")
axes[0].set_xlabel("temperature (K)")
axes[0].set_ylabel("pressure (mbar)")
axes[1].fill_betweenx(p_e / 100, -s_e, s_e, color="C3", alpha=0.18, linewidth=0, label="egress sigma_T")
axes[1].plot(at_pressure(p, T, p_e[common]) - T_e[common], p_e[common] / 100, color="k", linewidth=1.3,
             label="delivered - egress")
axes[1].plot(at_pressure(p_a, T_a, p_e[common]) - T_e[common], p_e[common] / 100, color="0.55", linewidth=0.9,
             linestyle="--", label="anchor - egress")
axes[1].plot(iris_s_all - at_pressure(p_e, T_e, levels_Pa), levels_Pa / 100, "s", color="C0", markersize=7,
             label=f"IRIS fit at {EGRESS_LABEL} - egress")
axes[1].axvline(0, color="0.6", linewidth=0.8)
axes[1].set_xlabel("temperature minus the egress (K)")
axes[2].plot(rho_e, p_e / 100, color="C3", linewidth=1.0, label="egress, 2.135 amu")
axes[2].plot(rho, p / 100, color="k", linewidth=1.3, label="delivered")
axes[2].plot(rho_a, p_a / 100, color="0.55", linewidth=0.9, linestyle="--", label="anchor")
axes[2].set_xscale("log")
axes[2].set_xlabel("mass density (kg/m3)")
for ax in axes:
    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
axes[0].set_ylim(float(max(p.max(), p_e.max())) / 100 * 1.1, float(min(p.min(), p_e.min())) / 100 / 1.1)
fig.suptitle(f"The first comparison: the Lindal anchor at {INGRESS_LABEL} transferred to {EGRESS_LABEL} "
             "under the Lindal-only wind, against the Voyager 2 egress", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "comparison.png", dpi=130)
plt.close(fig)
record(3, "the comparison: the delivered temperature and density against pressure, with the egress "
          "profile and its sigma_T, the anchor's and the IRIS fit's, as differences from the egress, and as "
          "changes between the two latitudes, for the author", True,
       f"the run rebuilt and rerun at {EGRESS_LABEL}, {float(root['latitude_planetocentric_deg'])!r} deg "
       f"planetocentric; "
       f"outer loop {passes} passes; pressure identity at most {np.abs(identity_K).max():.3f} K\n"
       + "\n".join(lines) + f"\nas changes from {INGRESS_LABEL} to {EGRESS_LABEL}:\n"
       + "\n".join(change_lines) + f"\nfigure {HERE / 'comparison.png'}")

(HERE / "summary.json").write_text(json.dumps({"levels": rows, "passes": passes}, indent=2), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
