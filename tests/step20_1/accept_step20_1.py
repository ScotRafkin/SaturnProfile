"""SPEC_20 Step 1 acceptance: the Voyager anchors, the sections and the anchor matrix.

Reads the anchors under `occul_data/lindal/`, their closure runs, the five sections under `sections/` and the
anchor matrix runs (`tests/step20_1/build_matrix.py`); runs one section check against SPEC_19's runs and one
linearization of `voyager_fit_1` for the shear's posterior sd. Writes under `reports/step20_1/`. Absolute
profiles come first in every figure set (author, 9 November).

Both dry adiabats, normal and equilibrium hydrogen, are drawn on every temperature profile panel and on the
lapse-rate figure, with `c_p(T)` from each anchor's kind C (`lib.heat_capacity`, SPEC_20 v0.5). Section figures
carry the zonal wind as dashed contours in m/s and no shear hatching (REVIEW_20_step1 ruling 4, standing).

Run from the repository root after the sections and the matrix; about five minutes.
"""

import csv
import json
import subprocess
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import xarray as xr  # noqa: E402

sys.path.insert(0, str(Path(__file__).parent))
import build_matrix as bm  # noqa: E402
from casspian.lib import geoid as gd  # noqa: E402
from casspian.lib import heat_capacity as hc  # noqa: E402
from casspian.lib import io as cio  # noqa: E402
from casspian.tools.lindal import radio_profiles as rp  # noqa: E402
from casspian.tools.section import section as sec  # noqa: E402
from casspian.tools.wind.estimation import driver  # noqa: E402
from casspian.tools.wind.estimation.radio import interp_ln_p, sigma_curve  # noqa: E402

warnings.filterwarnings("ignore")
HERE = Path("reports/step20_1")
FIG = Path("reports/figures")
LINDAL = Path("occul_data/lindal")
ANCHORS = ("lindal", "voyager2_egress", "voyager1_ingress", "voyager1_egress")
LABEL = {"lindal": "36.3 N (30.8 N), Lindal", "voyager2_egress": "31.2 S (26.1 S), V2 egress",
         "voyager1_ingress": "74.2 S (71.1 S), V1 ingress", "voyager1_egress": "3.0 S (2.4 S), V1 egress"}
COLOR = {"lindal": "C0", "voyager2_egress": "C3", "voyager1_ingress": "C2", "voyager1_egress": "C4"}
SECTIONS = {"lindal": "lindal_voyager_fit_1", "voyager2_egress": "voyager2_egress_voyager_fit_1",
            "voyager1_ingress": "voyager1_ingress_voyager_fit_1", "voyager1_egress": "voyager1_egress_voyager_fit_1"}
HYP1 = "lindal_iris_v_sl_ii_110"
TABLE_MBAR = (730.0, 290.0, 110.0, 30.0, 10.0, 1.0)
#: Pioneer 11's exit, planetocentric (SPEC_18), marked only.
PIONEER11 = -9.5
#: The common pressure grid of the section figures, Pa: 1.3 bar to 0.2 mbar.
GRID_PA = np.geomspace(1.3e5, 20.0, 160)
#: SPEC_15's calibration traces against Lindal's own data, RMS 0.95 and 1.51 K (REPORT_15_step1): the tracing's
#: precision, the bound of acceptance 2.
TRACING_K = 1.0
summary, results = {}, []


def record(number, description, passed, detail):
    results.append(passed)
    print(f"[{'REPORTED' if passed is None else 'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")
    sys.stdout.flush()


def tree(path):
    return driver.tree(path)


def observed(slug):
    """The profile the anchor was reduced from: Table I for Lindal, the trace for the others. mbar, K."""
    if slug == "lindal":
        th = xr.open_dataset(LINDAL / "lindal_thermo.nc")
        return th["pressure_Pa"].values / 100.0, th["temperature_K"].values
    t = rp.read_csv(LINDAL / f"{slug}.csv")
    good = np.isfinite(t["temperature_K"])
    return t["pressure_mbar"][good], t["temperature_K"][good]


def latitude(slug):
    return float(tree(LINDAL / f"{slug}_refractivity.nc").to_dataset(inherit=False)["latitude_planetocentric_deg"])


def on_grid(ds, name, grid=GRID_PA):
    """A section variable on the common pressure grid, `(latitude, grid)`, NaN outside the anchor's levels."""
    out = []
    for i in range(ds.sizes["latitude"]):
        p = ds["pressure_Pa"].values[i]
        out.append(interp_ln_p(p, ds[name].values[i], grid))
    return np.array(out)


ADIABAT_STYLE = {"normal": ("0.35", "--"), "equilibrium": ("0.35", ":")}
#: The adiabats are drawn from the observed profile's deepest point up to this pressure (SPEC_20 v0.4).
ADIABAT_TOP_Pa = 1.0e4


def composition(slug):
    """The anchor's kind C, the mixture its runs carry (94/6 with Lindal's ammonia, on its levels)."""
    handle = cio.read(LINDAL / f"{slug}_composition.nc", "composition")
    try:
        return handle.load()
    finally:
        handle.close()


def draw_adiabats(ax, slug, p_mbar, T, latitude_deg):
    deepest = int(np.argmax(p_mbar))
    for mode, (color, style) in ADIABAT_STYLE.items():
        p, t = hc.adiabat(float(T[deepest]), float(p_mbar[deepest]) * 100.0, ADIABAT_TOP_Pa, composition(slug), latitude_deg, mode)
        ax.plot(t, p / 100.0, color=color, linestyle=style, linewidth=1.1, label=f"dry adiabat, {mode} H2")


def main():
    HERE.mkdir(parents=True, exist_ok=True)
    sigma = sigma_curve(str(LINDAL / "voyager2_egress.csv"))
    lats = {s: latitude(s) for s in ANCHORS}

    # ---- 1. Lindal's answer unchanged --------------------------------------------------------------
    # Lindal's own files by name. The traced profiles beside them are SPEC_15's, and voyager2_egress.csv is
    # rebuilt by REVIEW_20_step1 ruling 2.
    paths = ["occul_data/lindal/lindal_*", "occul_data/lindal/raw/lindal_*", "forward/lindal_closure", "forward/lindal_transfer"]
    status = subprocess.run(["git", "status", "--porcelain", "--", *paths], capture_output=True, text=True).stdout.splitlines()
    touched = [line for line in status if not line.startswith("??")]
    record(1, "Lindal's files unchanged: no tracked lindal_* file under occul_data/lindal, nor lindal_closure or lindal_transfer, differs",
           not touched, "tracked changes: " + (", ".join(touched) if touched else "none") + f"; new files beside them: {len(status) - len(touched)}")

    # ---- 2. closures ------------------------------------------------------------------------------
    lines, ok = [f"delivered at its own latitude minus its own profile, on its levels; bound {TRACING_K} K, the tracing's precision"], True
    summary["closure"] = {}
    for slug in ANCHORS:
        r = tree(f"forward/{slug}_closure/output/{slug}_closure_profile.nc").to_dataset(inherit=False)
        dT = r["temperature_K"].values - r["temperature_tabulated_K"].values
        dp = r["pressure_Pa"].values / r["pressure_tabulated_Pa"].values - 1.0
        identity = float(np.nanmax(np.abs(dp)))
        summary["closure"][slug] = {"max_dT": float(np.nanmax(np.abs(dT))), "rms_dT": float(np.sqrt(np.nanmean(dT ** 2))), "max_dp": identity}
        ok &= slug == "lindal" or (np.nanmax(np.abs(dT)) <= TRACING_K and identity < 1e-2)
        lines.append(f"  {slug:17s} levels {dT.size:3d}  largest |dT| {np.nanmax(np.abs(dT)):.4f} K, RMS {np.sqrt(np.nanmean(dT ** 2)):.4f} K; "
                     f"largest |p / p_tabulated - 1| {identity:.2e}")
    lines.append("Lindal's own closure is shown for scale and is not new.")
    record(2, "each new anchor's closure within the tracing's precision, with the pressure identity", ok, "\n".join(lines))

    # ---- 3. the section against single runs ----------------------------------------------------------
    control = sec.read_control("sections/lindal_voyager_fit_1/lindal_voyager_fit_1.toml")
    runs = {lats[t]: bm.LINDAL_RUNS[t] for t in bm.LINDAL_RUNS}
    ds = sec.section(control, latitudes_deg=sorted(runs), work=control["here"] / "work_check")
    lines, worst = [], 0.0
    for lat, name in runs.items():
        run = tree(f"forward/{name}/output/{name}_profile.nc").to_dataset(inherit=False)
        row = ds.isel(latitude=int(np.flatnonzero(ds["latitude_planetocentric_deg"].values == lat)[0]))
        d = interp_ln_p(row["pressure_label_Pa"].values, row["temperature_K"].values, run["pressure_label_Pa"].values) - run["temperature_K"].values
        worst = max(worst, float(np.nanmax(np.abs(d))))
        lines.append(f"  {name}: largest |section - run| {np.nanmax(np.abs(d)):.2e} K, RMS {np.sqrt(np.nanmean(d ** 2)):.2e} K")
    summary["section_check_K"] = worst
    record(3, "the Lindal section under hypothesis 3 at 31.2 S, 74.2 S and 3.0 S against SPEC_19's runs, about 0.1 K", worst < 0.1, "\n".join(lines))

    # ---- 4. for the author, absolute first ----------------------------------------------------------------
    lines = []
    matrix = bm.runs()
    delivered = {pair: tree(f"forward/{run}/output/{run}_profile.nc").to_dataset(inherit=False) for pair, run in matrix.items()}

    def matrix_panel(ax, target, sources):
        p_obs, T_obs = observed(target)
        s = sigma(p_obs * 100.0)
        ax.fill_betweenx(p_obs, T_obs - s, T_obs + s, color="0.85", label="observed, sigma_T (Lindal Fig. 4)")
        ax.plot(T_obs, p_obs, color="k", linewidth=1.6, label=f"observed, {LABEL[target]}")
        draw_adiabats(ax, target, p_obs, T_obs, lats[target])
        for anchor in sources:
            r = delivered[(anchor, target)]
            ax.plot(r["temperature_K"].values, r["pressure_Pa"].values / 100.0, color=COLOR[anchor], linewidth=1.2,
                    label=f"from {LABEL[anchor]}")
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
        ax.set_xlabel("temperature (K)")
        ax.set_title(f"at {LABEL[target]}", fontsize=10)
        ax.grid(True, color="0.92", linewidth=0.6)
        ax.legend(fontsize=7, loc="lower left")

    fig, axes = plt.subplots(1, 4, figsize=(20, 8), sharey=True)
    for ax, target in zip(axes, ANCHORS):
        matrix_panel(ax, target, [a for a in ANCHORS if a != target])
    axes[0].set_ylabel("pressure (mbar)")
    fig.suptitle("The anchor matrix under hypothesis 3: absolute profiles", fontsize=12)
    fig.tight_layout()
    fig.savefig(FIG / "step20_1_matrix_absolute.png", dpi=100)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 8))
    matrix_panel(ax, "voyager1_ingress", ["voyager1_egress"])
    ax.set_ylabel("pressure (mbar)")
    ax.set_title("Voyager 1 egress to Voyager 1 ingress, hypothesis 3", fontsize=10)
    fig.tight_layout()
    fig.savefig(FIG / "step20_1_v1_egress_to_ingress.png", dpi=110)
    plt.close(fig)

    lines.append(f"delivered minus observed (K) at {', '.join(f'{v:g}' for v in TABLE_MBAR)} mbar, every pair; sigma_T last:")
    at = np.array(TABLE_MBAR) * 100.0
    summary["matrix"] = {}
    for target in ANCHORS:
        p_obs, T_obs = observed(target)
        obs = interp_ln_p(p_obs * 100.0, T_obs, at)
        lines.append(f"  at {LABEL[target]}:")
        for anchor in ANCHORS:
            if anchor == target:
                continue
            r = delivered[(anchor, target)]
            d = interp_ln_p(r["pressure_Pa"].values, r["temperature_K"].values, at) - obs
            summary["matrix"][f"{anchor} to {target}"] = d.tolist()
            lines.append(f"    from {anchor:17s}" + "".join(f"{v:+8.2f}" for v in d))
        lines.append(f"    {'sigma_T':22s}" + "".join(f"{v:8.2f}" for v in sigma(at)))
    lines.append("  Voyager 1 egress to Voyager 1 ingress, on its own: "
                 + " ".join(f"{v:+.2f}" for v in summary["matrix"]["voyager1_egress to voyager1_ingress"]))

    # The sections.
    sections = {slug: xr.open_dataset(f"sections/{name}/{name}_section.nc") for slug, name in SECTIONS.items()}
    sections["hypothesis 1"] = xr.open_dataset(f"sections/{HYP1}/{HYP1}_section.nc")
    T_grid = {k: on_grid(v, "temperature_K") for k, v in sections.items()}
    lat_s = sections["lindal"]["latitude_planetocentric_deg"].values

    wind = driver._read("fits/voyager_fit_1/voyager_fit_1_wind.nc", "wind")
    u_grid = {k: on_grid(v, "u_ms") for k, v in sections.items()}
    levels_u = [-50, 0, 100, 200, 300, 400, 500]

    iris = np.genfromtxt(LINDAL / "iris_temperatures.csv", delimiter=",", names=True)
    lat_c = np.asarray(wind["latitude_planetocentric_deg"].values)
    lat_g = np.asarray(wind["latitude_planetographic_deg"].values)
    iris_c = np.interp(iris["latitude_planetographic_deg"], lat_g, lat_c)
    levels_T = np.arange(70, 182, 4)

    def decorate(ax):
        for slug in ANCHORS:
            ax.axvline(lats[slug], color="k", linewidth=0.7, linestyle="--")
        ax.axvline(PIONEER11, color="k", linewidth=0.7, linestyle=":")
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
        ax.set_xlim(-90, 90)

    names = ["lindal", "voyager2_egress", "voyager1_ingress", "voyager1_egress", "hypothesis 1"]
    fig, axes = plt.subplots(2, 3, figsize=(19, 11), sharex=True, sharey=True)
    for ax, key in zip(axes.ravel(), names):
        cs = ax.contourf(lat_s, GRID_PA / 100, T_grid[key].T, levels=levels_T, cmap="viridis", extend="both")
        ax.contour(lat_s, GRID_PA / 100, T_grid[key].T, levels=levels_T[::2], colors="w", linewidths=0.4)
        ax.scatter(iris_c, iris["pressure_mbar"], c=iris["temperature_K"], cmap="viridis", vmin=levels_T[0], vmax=levels_T[-1],
                   s=14, edgecolors="w", linewidths=0.4)
        cu = ax.contour(lat_s, GRID_PA / 100, u_grid[key].T, levels=levels_u, colors="k", linestyles="--", linewidths=0.7)
        ax.clabel(cu, fmt="%d", fontsize=7, inline_spacing=2)
        decorate(ax)
        title = f"{LABEL[key]} anchor, hypothesis 3" if key in LABEL else "Lindal anchor, hypothesis 1"
        ax.set_title(title, fontsize=10)
    axes[1, 2].axis("off")
    fig.colorbar(cs, ax=axes, label="temperature (K)", shrink=0.8)
    for ax in axes[1]:
        ax.set_xlabel("planetocentric latitude (deg)")
    for ax in axes[:, 0]:
        ax.set_ylabel("pressure (mbar)")
    fig.suptitle("Temperature sections, the run's zonal wind in dashed black contours (m/s). Dots: IRIS at 110, 290, "
                 "730 mbar, same colors. Vertical lines: the four Voyager anchors (dashed), Pioneer 11 (dotted).", fontsize=10)
    fig.savefig(FIG / "step20_1_sections_temperature.png", dpi=95)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(19, 6.5), sharey=True)
    for ax, key in zip(axes, ANCHORS[1:]):
        cs = ax.contourf(lat_s, GRID_PA / 100, (T_grid[key] - T_grid["lindal"]).T, levels=np.arange(-20, 21, 2), cmap="RdBu_r", extend="both")
        cu = ax.contour(lat_s, GRID_PA / 100, u_grid[key].T, levels=levels_u, colors="k", linestyles="--", linewidths=0.7)
        ax.clabel(cu, fmt="%d", fontsize=7, inline_spacing=2)
        decorate(ax)
        ax.set_title(f"{LABEL[key]} section minus Lindal's", fontsize=10)
        ax.set_xlabel("planetocentric latitude (deg)")
    axes[0].set_ylabel("pressure (mbar)")
    fig.colorbar(cs, ax=axes, label="K")
    fig.savefig(FIG / "step20_1_sections_difference.png", dpi=95)
    plt.close(fig)
    summary["section_difference_K"] = {}
    for key in ANCHORS[1:]:
        d = T_grid[key] - T_grid["lindal"]
        rows = []
        for lo, hi in ((1.3e5, 2e4), (2e4, 3e3), (3e3, 20.0)):
            band = (GRID_PA <= lo) & (GRID_PA > hi)
            rows.append((float(np.nanmean(np.abs(d[:, band]))), float(np.nanmax(np.abs(d[:, band])))))
        summary["section_difference_K"][key] = rows
        lines.append(f"section minus Lindal's, {key}: mean |d| and largest, K, 1300 to 200, 200 to 30, 30 to 0.2 mbar: "
                     + "; ".join(f"{m:.1f}, {x:.1f}" for m, x in rows))

    p_w = np.asarray(wind["pressure_Pa"].values)
    keep = (p_w <= 1.3e5) & (p_w >= 20.0)
    fig, axes = plt.subplots(1, 2, figsize=(15, 6.5), sharey=True)
    cs = axes[0].contourf(lat_c, p_w[keep] / 100, wind["u_total_ms"].values[:, keep].T, levels=np.arange(-200, 601, 25), cmap="RdBu_r")
    fig.colorbar(cs, ax=axes[0], label="u (m/s)")
    cs = axes[1].contourf(lat_c, p_w[keep] / 100, wind["u_total_uncertainty_ms"].values[:, keep].T, levels=np.arange(0, 210, 10), cmap="magma_r")
    fig.colorbar(cs, ax=axes[1], label="posterior sd (m/s)")
    for ax, title in zip(axes, ("hypothesis 3 zonal wind", "its posterior sd")):
        decorate(ax)
        ax.set_title(title)
        ax.set_xlabel("planetocentric latitude (deg)")
    axes[0].set_ylabel("pressure (mbar)")
    fig.savefig(FIG / "step20_1_wind.png", dpi=95)
    plt.close(fig)

    # The lapse rate at each anchor's own latitude, against both dry adiabatic lapse rates, g m_bar / c_p, with g
    # from the section itself, g = -(dp/dz) / rho.
    fig, axes = plt.subplots(1, 4, figsize=(21, 7), sharey=True)
    summary["lapse"] = {}
    lines.append("lapse rate at each anchor's latitude against the dry adiabatic lapse rate, below 500 mbar: median, then the "
                 "levels above each adiabat (superadiabatic) and the largest excess, K/km:")
    for ax, key in zip(axes, ANCHORS):
        ds = sections[key]
        i = int(np.argmin(np.abs(ds["latitude_planetocentric_deg"].values - lats[key])))
        p, T = ds["pressure_Pa"].values[i], ds["temperature_K"].values[i]
        z, rho = ds["height_above_anchor_isobar_m"].values[i], ds["density_kg_m3"].values[i]
        m = ds["mean_molar_mass_kg_mol"].values[i]
        g = -np.gradient(p, z) / rho
        lapse = ds["lapse_rate_K_km"].values[i]
        ax.plot(lapse, p / 100, color=COLOR[key], linewidth=1.0, label="-dT/dz, delivered")
        deep = p >= 5e4
        text = [f"  {key} at {ds['latitude_planetocentric_deg'].values[i]:+.1f}: median {np.median(lapse[deep]):.3f}"]
        for mode, (color, style) in ADIABAT_STYLE.items():
            gamma = g * m / hc.mixture_cp(composition(key), p, lats[key], T, mode) * 1e3
            ax.plot(gamma, p / 100, color="k", linestyle=style, linewidth=1.2, label=f"dry adiabatic, {mode} H2")
            over = deep & (lapse > gamma)
            summary["lapse"].setdefault(key, {})[mode] = {"median_adiabat": float(np.median(gamma[deep])), "levels_over": int(over.sum()),
                                                          "levels": int(deep.sum()), "largest_excess": float(np.max((lapse - gamma)[deep]))}
            text.append(f"{mode} adiabat {np.median(gamma[deep]):.3f}, over it {over.sum()} of {deep.sum()} levels, largest excess "
                        f"{np.max((lapse - gamma)[deep]):+.3f}")
        lines.append("; ".join(text))
        ax.set_yscale("log")
        ax.set_ylim(1300, 0.2)
        ax.set_xlim(-3, 2.5)
        ax.axvline(0, color="0.7", linewidth=0.6)
        ax.set_title(f"lapse rate at {LABEL[key]}, its own anchor", fontsize=10)
        ax.set_xlabel("-dT/dz (K/km)")
        ax.grid(True, color="0.92", linewidth=0.6)
        ax.legend(fontsize=8, loc="upper left")
    axes[0].set_ylabel("pressure (mbar)")
    fig.tight_layout()
    fig.savefig(FIG / "step20_1_lapse_rate.png", dpi=95)
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(17, 5), sharex=True)
    summary["density"] = {}
    for ax, p_at in zip(axes, (1e5, 1e4, 1e3)):
        for key in ANCHORS:
            rho = on_grid(sections[key], "density_kg_m3", np.array([p_at]))[:, 0]
            ax.plot(lat_s, rho, color=COLOR[key], label=LABEL[key])
            summary["density"].setdefault(key, {})[f"{p_at:g}"] = [float(np.nanmin(rho)), float(np.nanmax(rho))]
        for slug in ANCHORS:
            ax.axvline(lats[slug], color="0.6", linewidth=0.6, linestyle="--")
        ax.set_title(f"density at {p_at / 100:g} mbar, hypothesis 3")
        ax.set_xlabel("planetocentric latitude (deg)")
        ax.set_ylabel("kg/m3")
    axes[0].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "step20_1_density.png", dpi=100)
    plt.close(fig)

    # Fig. 9 as a check of the Voyager 1 radii (REVIEW_20_step0 ruling 1): each anchor's 100 mbar radius above
    # Lindal's reference geoid against his point.
    g, r = driver._read(LINDAL / "lindal_gravity.nc", "gravity"), driver._read(LINDAL / "lindal_rotation.nc", "rotation")
    lat_ref = np.arange(-90.0, 90.0001, 0.05)
    r_ref = gd.reference_geoid(np.radians(lat_ref), 54438.0e3, float(r["angular_rate_rad_s"]), float(g["GM_m3s2"]),
                               np.asarray(g["J"].values), np.asarray(g["degree"].values), float(g["normalization_radius_m"]))[0]
    points = [row for row in csv.DictReader(line for line in open(LINDAL / "isobar_100mbar.csv", encoding="utf-8") if not line.startswith("#"))]
    lines.append("Fig. 9 check, each anchor's 100 mbar altitude above Lindal's reference geoid against his nearest point, km:")
    for slug in ANCHORS:
        n = tree(LINDAL / f"{slug}_refractivity.nc").to_dataset(inherit=False)
        alt = (float(n["anchor_isobar_radius_m"]) - np.interp(lats[slug], lat_ref, r_ref)) / 1e3
        near = min(points, key=lambda row: abs(float(row["latitude_planetocentric_deg"]) - lats[slug]))
        lines.append(f"  {slug:17s} at {lats[slug]:+.2f}: ours {alt:6.1f}; Fig. 9 {float(near['altitude_km']):6.1f} at "
                     f"{float(near['latitude_planetocentric_deg']):+.2f}; difference {alt - float(near['altitude_km']):+.1f}")

    lines.append("figures, absolute first: reports/figures/step20_1_matrix_absolute.png, step20_1_v1_egress_to_ingress.png, "
                 "step20_1_sections_temperature.png, step20_1_sections_difference.png, step20_1_wind.png, step20_1_lapse_rate.png, "
                 "step20_1_density.png")
    record(4, "for the author, absolute first: the matrix, the sections, their differences, the wind, the lapse rate, density", None, "\n".join(lines))
    (HERE / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
    bounded = [r for r in results if r is not None]
    print(f"{sum(bounded)} of {len(bounded)} bounded checks pass; {results.count(None)} reported without a bound")
    return 0 if all(bounded) else 1


if __name__ == "__main__":
    raise SystemExit(main())
