"""Acceptance for SPEC_15 v0.2 Step 1: three more radio profiles, and the two hypotheses at four targets.

Traces Lindal et al. (1985) Figs. 5 and 6 with `tools/lindal/radio_profiles.py`, compares the
calibration traces and the Pioneer 11 check trace, and reads the committed data files. Builds the six
new runs with `casspian-run-inputs` and runs them with the transfer driver where their products are
absent (31.2 S is SPEC_14's two runs, read as they are, built if absent). Run inputs and products are
written in the run directories, as a run's are; this script writes under `reports/step15_1/`.

Latitudes are given planetographic, with planetocentric in parentheses.

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
from casspian.lib.constants import AVOGADRO_CONSTANT, MOLAR_GAS_CONSTANT
from casspian.tools.lindal import iris_temperatures as iris
from casspian.tools.lindal import radio_profiles as rp
from casspian.tools.run import run_inputs

HERE = Path("reports/step15_1")
HERE.mkdir(parents=True, exist_ok=True)
PDF = "docs/Lindal_et_al_1985_AJ90_1136.pdf"
TABLE_I = Path("occul_data/lindal/raw/lindal_table1.csv")
EGRESS_V2 = Path("occul_data/lindal/voyager2_egress.csv")
IRIS_TABLE = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
EGRESS_AMU = 2.135
LEVELS_MBAR = (730.0, 290.0, 110.0, 10.0, 4.5, 1.0, 0.3)
IRIS_LEVELS = {110: 11000.0, 290: 29000.0, 730: 73000.0}
#: Each target: its label, its planetocentric latitude (target_latitudes.py), its observed profile,
#: and its two runs, hypothesis 1 then hypothesis 2.
TARGETS = {
    "31.2 S": {"label_g": -31.2, "observed": "Voyager 2 egress (Lindal Fig. 4, SPEC_12)",
               "runs": ("lindal_iris_v_sl_ii_110", "lindal_iris_v_sl_iii_110")},
    "74.2 S": {"label_g": -74.2, "observed": "Voyager 1 ingress (Lindal Fig. 5)",
               "runs": ("lindal_iris_v_sl_ii_110_74s", "lindal_iris_v_sl_iii_110_74s")},
    "3.0 S": {"label_g": -3.0, "observed": "Voyager 1 egress (Lindal Fig. 6)",
              "runs": ("lindal_iris_v_sl_ii_110_03s", "lindal_iris_v_sl_iii_110_03s")},
    "11.9 S": {"label_g": -11.9, "observed": "Pioneer 11 exit (Kliore et al. 1980 Table 1, rescaled)",
               "runs": ("lindal_iris_v_sl_ii_110_12s", "lindal_iris_v_sl_iii_110_12s")},
}
HYPOTHESIS = ("hypothesis 1 (case ii)", "hypothesis 2 (case iii)")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def at_pressure(p, values, targets):
    p, values = np.asarray(p, dtype="float64"), np.asarray(values, dtype="float64")
    order = np.argsort(np.log(p))
    out = np.interp(np.log(targets), np.log(p)[order], values[order], left=np.nan, right=np.nan)
    return out


def compare(trace_p, trace_T, ref_p, ref_T):
    """The trace against a reference on the reference's levels inside the trace: RMS, largest and its
    pressure, in K."""
    d = at_pressure(trace_p, trace_T, ref_p) - ref_T
    ok = np.isfinite(d)
    k = int(np.nanargmax(np.abs(d)))
    return float(np.sqrt(np.mean(d[ok] ** 2))), float(d[k]), float(ref_p[k]), int(ok.sum()), d


warnings.simplefilter("ignore")
summary = {}

# ---------------------------------------------------------------------------
# 1. The calibration traces
# ---------------------------------------------------------------------------
ink = rp.render(PDF)
calibrations, traces = rp.trace_all(ink)
table_i = np.genfromtxt(TABLE_I, delimiter=",", names=True)
egress_v2 = np.genfromtxt(EGRESS_V2, delimiter=",", names=True)
lines = []
for name, ref, label in (("36.3N", (table_i["pressure_mbar"], table_i["temperature_K"]), "Lindal Table I"),
                         ("31.2S", (egress_v2["pressure_mbar"], egress_v2["temperature_K"]),
                          "voyager2_egress.csv (Fig. 4, SPEC_12)")):
    t = traces[name]
    rms, worst, p_worst, n, d = compare(t.pressure_mbar, t.temperature_K, *ref)
    summary[f"calibration_{name}"] = {"rms_K": rms, "largest_K": worst, "at_mbar": p_worst, "n": n}
    big = np.flatnonzero(np.abs(d) > 1.0)
    lines.append(f"Fig. 5 {name} trace against {label}: {n} levels, RMS {rms:.2f} K, largest {worst:+.2f} K "
                 f"at {p_worst:.4g} mbar; {big.size} levels over 1 K"
                 + (f", at {', '.join(f'{ref[0][i]:.3g}' for i in big[:12])} mbar" if big.size else ""))
for name, cal in calibrations.items():
    lines.append(f"{name} calibration: ticks within {np.abs(cal.p_residual_dex).max():.4f} dex in pressure and "
                 f"{np.abs(cal.t_residual_K).max():.3f} K in temperature of the fits")
coincident = rp.coincident_ranges(traces["36.3N"], traces["31.2S"], calibrations["fig5"])
lines.append("36.3 N and 31.2 S within one line width (both on the shared ink): "
             + ("; ".join(f"{a:.3g} to {b:.3g} mbar" for a, b in coincident) or "none"))
for pair in (("74.2S", "36.3N"), ("74.2S", "31.2S")):
    ranges = rp.coincident_ranges(traces[pair[0]], traces[pair[1]], calibrations["fig5"])
    lines.append(f"{pair[0]} and {pair[1]} within one line width: "
                 + ("; ".join(f"{a:.3g} to {b:.3g} mbar" for a, b in ranges) or "none"))
ranges = rp.coincident_ranges(traces["3.0S"], traces["11.9S"], calibrations["fig6"])
lines.append("3.0S and 11.9S within one line width: " + ("; ".join(f"{a:.3g} to {b:.3g} mbar" for a, b in ranges) or "none"))
for name, t in traces.items():
    lines.append(f"{name}: {t.pressure_mbar.size} levels, {t.pressure_mbar.min():.3g} to {t.pressure_mbar.max():.4g} mbar, "
                 f"{int(t.filled.sum())} bins filled across a gap; bars "
                 + ("; ".join(f"{p:.3g} mbar +-{h:.2f} K (centre {c:.1f} K)" for p, h, c in t.bars) or "none"))
record(1, "the calibration traces: Fig. 5's 36.3 N against Lindal's Table I and 31.2 S against the Fig. 4 egress "
          "(reported; over 1 K explained in the report)", True, "\n".join(lines))

# ---------------------------------------------------------------------------
# 2. Pioneer 11, two ways
# ---------------------------------------------------------------------------
kliore = rp.read_pioneer11()
T_r, p_r, factors = rp.rescale_composition(kliore["temperature_K"], kliore["pressure_mbar"])
t = traces["11.9S"]
rms_pub, worst_pub, p_pub, n_pub, d_pub = compare(t.pressure_mbar, t.temperature_K, kliore["pressure_mbar"], kliore["temperature_K"])
rms_res, worst_res, p_res, n_res, d_res = compare(t.pressure_mbar, t.temperature_K, p_r, T_r)
summary["pioneer11"] = {"published": {"rms_K": rms_pub, "largest_K": worst_pub, "at_mbar": p_pub, "mean_K": float(np.nanmean(d_pub))},
                        "rescaled": {"rms_K": rms_res, "largest_K": worst_res, "at_mbar": p_res, "mean_K": float(np.nanmean(d_res))},
                        "factors": factors}
record(2, "Pioneer 11 two ways: the Fig. 6 stippled trace against Kliore's Table 1 as published (90/10) and "
          "rescaled to 94/6 by Lindal's rule (reported)", True,
       f"rescaling factors: temperature x{factors['temperature']:.4f} ({100 * (factors['temperature'] - 1):+.2f} %), "
       f"pressure x{factors['pressure']:.4f} (refractivity x{factors['refractivity']:.4f})\n"
       f"trace minus table as published: {n_pub} levels, mean {np.nanmean(d_pub):+.2f} K, RMS {rms_pub:.2f} K, "
       f"largest {worst_pub:+.2f} K at {p_pub:.4g} mbar\n"
       f"trace minus table rescaled: {n_res} levels, mean {np.nanmean(d_res):+.2f} K, RMS {rms_res:.2f} K, "
       f"largest {worst_res:+.2f} K at {p_res:.4g} mbar\n"
       f"Lindal rescaled the table if the second agrees better than the first.")

# ---------------------------------------------------------------------------
# 3. The tracing overlays
# ---------------------------------------------------------------------------
colours = {"36.3N": (220, 0, 0), "31.2S": (0, 150, 0), "74.2S": (0, 0, 230), "3.0S": (220, 0, 0), "11.9S": (0, 150, 0)}
page = Image.fromarray(np.where(ink, 0, 255).astype(np.uint8)).convert("RGB")
draw = ImageDraw.Draw(page)
for name, trace in traces.items():
    cal = calibrations[trace.figure]
    xy = list(zip(cal.column(trace.temperature_K), cal.row(trace.pressure_mbar)))
    draw.line(xy, fill=colours[name], width=4)
for name, cal in calibrations.items():
    f = cal.frame
    page.crop((f["left"] - 60, f["top"] - 60, f["right"] + 60, f["bottom"] + 60)).save(HERE / f"overlay_{name}.png")
record(3, "the tracing overlays, each figure with its traces drawn over it, on the binned values (for the author)", True,
       f"{HERE / 'overlay_fig5.png'}, {HERE / 'overlay_fig6.png'}; 36.3 N red, 31.2 S green, 74.2 S blue on Fig. 5; "
       "3.0 S red, 11.9 S green on Fig. 6")

# ---------------------------------------------------------------------------
# The runs
# ---------------------------------------------------------------------------
products = {}
built = []
for target in TARGETS.values():
    for name in target["runs"]:
        run = Path("forward") / name
        path = run / "output" / f"{name}_profile.nc"
        if not path.exists():
            run_inputs.build(run / f"{name}_build.toml")
            path = forward_transfer.run(run / f"{name}.toml").product
            built.append(name)
        products[name] = tree(path)

# ---------------------------------------------------------------------------
# 4. The pressure identity of each run
# ---------------------------------------------------------------------------
lines = []
for name, product in products.items():
    root = product.to_dataset(inherit=False)
    p, T = np.asarray(root["pressure_Pa"].values), np.asarray(root["temperature_K"].values)
    e = T * np.asarray(root["pressure_identity_residual"].values)
    k = int(np.argmax(np.abs(e)))
    with open(HERE / f"pressure_identity_{name}.csv", "w", encoding="utf-8", newline="\n") as handle:
        handle.write("pressure_mbar,temperature_K,identity_K\n")
        handle.writelines(f"{pp / 100:.6g},{tt:.4f},{ee:+.5f}\n" for pp, tt, ee in zip(p, T, e))
    summary.setdefault("identity_K", {})[name] = float(np.abs(e).max())
    lines.append(f"{name}: largest {abs(e[k]):.4f} K at {p[k] / 100:.4g} mbar, RMS {np.sqrt(np.mean(e ** 2)):.4f} K over {e.size} levels"
                 + (" (built here)" if name in built else ""))
record(4, "the pressure identity of each run, in kelvins, level by level (reported)", True,
       "\n".join(lines) + f"\nevery level in {HERE}/pressure_identity_<run>.csv")

# ---------------------------------------------------------------------------
# 5. For the author, at each of the four targets
# ---------------------------------------------------------------------------
v1_in = rp.read_csv("occul_data/lindal/voyager1_ingress.csv")
v1_out = rp.read_csv("occul_data/lindal/voyager1_egress.csv")
sigma_r = 0.5 * (kliore["temperature_uncertainty_plus_K"] - kliore["temperature_uncertainty_minus_K"]) * factors["temperature"]
observed = {
    "31.2 S": (egress_v2["pressure_mbar"], egress_v2["temperature_K"], egress_v2["temperature_uncertainty_K"]),
    "74.2 S": (v1_in["pressure_mbar"], v1_in["temperature_K"], v1_in["temperature_uncertainty_K"]),
    "3.0 S": (v1_out["pressure_mbar"], v1_out["temperature_K"], v1_out["temperature_uncertainty_K"]),
    "11.9 S": (p_r, T_r, sigma_r),
}
table_lines = ["delivered minus observed, K, at " + ", ".join(f"{l:g}" for l in LEVELS_MBAR) + " mbar (NaN where the profile has no data):"]
fig_s, axes_s = plt.subplots(2, 4, figsize=(20, 10), sharey=True)
for col, (where, spec) in enumerate(TARGETS.items()):
    p_o, T_o, s_o = (np.asarray(v, dtype="float64") for v in observed[where])
    order = np.argsort(p_o)
    p_o, T_o, s_o = p_o[order], T_o[order], s_o[order]
    rho_o = p_o * 100.0 * EGRESS_AMU * 1e-3 / (MOLAR_GAS_CONSTANT * T_o)
    first = products[spec["runs"][0]]
    anchor_thermo = first["anchors/lindal/inputs/thermo"].to_dataset(inherit=False)
    p_a, T_a = np.asarray(anchor_thermo["pressure_Pa"].values) / 100, np.asarray(anchor_thermo["temperature_K"].values)
    target_c = float(first.to_dataset(inherit=False)["latitude_planetocentric_deg"])
    title = f"{where} ({abs(target_c):.2f} {'S' if target_c < 0 else 'N'} planetocentric)"
    fig, axes = plt.subplots(1, 3, figsize=(17, 7), sharey=True)
    bars = np.isfinite(s_o)
    axes[0].plot(T_o, p_o, color="C3", linewidth=1.0, label=f"observed: {spec['observed']}")
    axes[0].errorbar(T_o[bars], p_o[bars], xerr=s_o[bars], fmt="none", ecolor="C3", capsize=3)
    axes[0].plot(T_a, p_a, color="0.3", linewidth=1.0, linestyle="--", label="anchor, Lindal ingress, 36.3 N (30.81 N planetocentric)")
    axes[1].errorbar(np.zeros(int(bars.sum())), p_o[bars], xerr=s_o[bars], fmt="none", ecolor="C3", capsize=3,
                     label="observed uncertainty")
    axes[2].plot(rho_o, p_o, color="C3", linewidth=1.0, label="observed, 2.135 amu")
    rows = {}
    for k, name in enumerate(spec["runs"]):
        root = products[name].to_dataset(inherit=False)
        p, T = np.asarray(root["pressure_Pa"].values) / 100, np.asarray(root["temperature_K"].values)
        rho = np.asarray(root["number_density_m3"].values) * np.asarray(root["mean_molar_mass_kg_mol"].values) / AVOGADRO_CONSTANT
        colour = ("C2", "C4")[k]
        axes[0].plot(T, p, color=colour, linewidth=1.2, label=f"{name}, {HYPOTHESIS[k]}")
        inside = (p_o >= p.min()) & (p_o <= p.max())
        diff = at_pressure(p, T, p_o[inside]) - T_o[inside]
        axes[1].plot(diff, p_o[inside], color=colour, linewidth=1.2, label=f"{HYPOTHESIS[k]} - observed")
        axes_s[k, col].plot(diff, p_o[inside], color=colour, linewidth=1.2)
        axes[2].plot(rho, p, color=colour, linewidth=1.2, label=name)
        rows[HYPOTHESIS[k]] = at_pressure(p, T, LEVELS_MBAR) - at_pressure(p_o, T_o, LEVELS_MBAR)
    a_inside = (p_o >= p_a.min()) & (p_o <= p_a.max())
    anchor_diff = at_pressure(p_a, T_a, p_o[a_inside]) - T_o[a_inside]
    axes[1].plot(anchor_diff, p_o[a_inside], color="0.3", linestyle="--", linewidth=1.0, label="anchor - observed (untransferred)")
    rows["anchor, untransferred"] = at_pressure(p_a, T_a, LEVELS_MBAR) - at_pressure(p_o, T_o, LEVELS_MBAR)
    iris_at = []
    for level, pa in IRIS_LEVELS.items():
        r = IRIS_TABLE["pressure_mbar"] == level
        v = iris.fit(IRIS_TABLE["latitude_planetographic_deg"][r], IRIS_TABLE["temperature_K"][r],
                     grid=np.array([spec["label_g"]]), window="smooth").value_K[0]
        iris_at.append(v)
    iris_at = np.array(iris_at)
    if np.isfinite(iris_at).any():
        axes[0].plot(iris_at, np.array(list(IRIS_LEVELS.values())) / 100, "s", color="C9", markersize=7,
                     label=f"IRIS fit (smooth window) at {where}")
        axes[1].plot(iris_at - at_pressure(p_o, T_o, np.array(list(IRIS_LEVELS.keys()), dtype=float)),
                     np.array(list(IRIS_LEVELS.keys())), "s", color="C9", markersize=7, label="IRIS fit - observed")
    for ax in axes:
        ax.set_yscale("log")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    axes[0].invert_yaxis()
    axes[1].axvline(0, color="0.6", linewidth=0.8)
    axes[0].set_xlabel("temperature (K)")
    axes[0].set_ylabel("pressure (mbar)")
    axes[1].set_xlabel("temperature minus observed (K)")
    axes[2].set_xscale("log")
    axes[2].set_xlabel("mass density (kg/m3)")
    fig.suptitle(f"The two hypotheses at {title} against the {spec['observed']}", fontsize=10)
    fig.tight_layout()
    tag = where.replace(" ", "").replace(".", "")
    fig.savefig(HERE / f"profiles_{tag}.png", dpi=120)
    plt.close(fig)
    sigma_at = at_pressure(p_o, np.where(np.isfinite(s_o), s_o, np.nan), LEVELS_MBAR)
    table_lines.append(f"  {title}, against the {spec['observed']}:")
    for label, values in rows.items():
        table_lines.append(f"    {label:28s}" + "  ".join(f"{v:+7.2f}" for v in values))
    table_lines.append(f"    {'IRIS fit - observed':28s}" + ", ".join(
        f"{lev} mbar {iris_at[i] - at_pressure(p_o, T_o, [float(lev)])[0]:+.2f}" for i, lev in enumerate(IRIS_LEVELS)))
    summary.setdefault("minus_observed", {})[where] = {k: [float(x) for x in v] for k, v in rows.items()}
    for k in range(2):
        axes_s[k, col].axvline(0, color="0.6", linewidth=0.8)
        axes_s[k, col].set_yscale("log")
        axes_s[k, col].grid(alpha=0.3)
        axes_s[k, col].set_title(f"{HYPOTHESIS[k]}, {title}", fontsize=8)
        axes_s[k, col].set_xlim(-30, 30)
axes_s[0, 0].invert_yaxis()
for ax in axes_s[1]:
    ax.set_xlabel("delivered minus observed (K)")
for ax in axes_s[:, 0]:
    ax.set_ylabel("pressure (mbar)")
fig_s.suptitle("Delivered minus observed against pressure, four targets by two hypotheses (hypothesis 2's first guesses were scaled from 31.2 S)", fontsize=10)
fig_s.tight_layout()
fig_s.savefig(HERE / "summary.png", dpi=110)
plt.close(fig_s)
record(5, "for the author: at each of the four targets, the profiles of both hypotheses with the anchor, the observed "
          "profile and its uncertainties and the IRIS fit; the table of delivered minus observed; the summary figure", True,
       "\n".join(table_lines) + f"\nfigures {HERE}/profiles_<target>.png, {HERE / 'summary.png'}")

(HERE / "summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8", newline="\n")
sys.exit(1 if failed else 0)
