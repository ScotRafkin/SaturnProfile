"""Acceptance for SPEC_18 v0.3 Step 1: Kliore et al. (1980) as a second registered profile.

Reads the built chain under `occul_data/kliore1980/` and the closure under
`forward/kliore1980_closure/`, which the runbook builds, and rebuilds Lindal's files in a copy under
`reports/step18_1/lindal_copy/`. Writes under `reports/step18_1/` only.

Latitudes are given planetographic, with planetocentric in parentheses: 11.9 S (9.5 S).

Run from the repository root.
"""

import json
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.constants import MOLAR_GAS_CONSTANT
from casspian.lib.gravity import g_eff_vector, validated_degrees
from casspian.lib.schema import WRITER_FILLED_GLOBALS
from casspian.tools.composition import build_composition
from casspian.tools.lindal import build_inputs

warnings.filterwarnings("ignore", message=".*input_hashes.*")

HERE = Path("reports/step18_1")
HERE.mkdir(parents=True, exist_ok=True)
K = Path("occul_data/kliore1980")
L = Path("occul_data/lindal")
RUN = Path("forward/kliore1980_closure")
#: Half a unit of what Table 1 prints: 0.1 mbar and 0.1 km. SPEC_16: 0.1 n-units.
HALF_P_PA = 5.0
HALF_R_M = 50.0
HALF_N = 0.05
ANCHOR_PA = 1.0e4
#: SPEC_16's fitted offset, for the report's comparison.
SPEC_16_OFFSET_KM = 9.2
results = []
summary = {}


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


def tree(path):
    handle = xr.open_datatree(path, engine="netcdf4").load()
    handle.close()
    return handle


def nodes(path):
    return {n.path: n.to_dataset(inherit=False) for n in tree(path).subtree}


def same(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape:
        return False
    return np.array_equal(a, b, equal_nan=True) if a.dtype.kind == "f" else np.array_equal(a, b)


#: What a rebuild need not reproduce: what the writer fills, and the lineage.
NOT_VALUE = set(WRITER_FILLED_GLOBALS) | {"history", "control_file", "input_hashes", "raw_bundle",
                                          "latitude_conversion_inputs"}


def compare(registered, rebuilt):
    """Variables array-equal; every attribute the rebuild writes equal. A global only the registered
    file has is one SPEC_R1 stopped writing, and is counted, not compared."""
    old, new = nodes(registered), nodes(rebuilt)
    count, value_diffs, attr_diffs, retired = 0, [], [], 0
    if set(old) != set(new):
        value_diffs.append(f"groups {set(old) ^ set(new)}")
    for group in set(old) & set(new):
        o, n = old[group], new[group]
        for name in set(o.variables) | set(n.variables):
            count += 1
            if name not in o.variables or name not in n.variables:
                value_diffs.append(f"{group} {name}: in one only")
            elif not same(o[name].values, n[name].values):
                value_diffs.append(f"{group} {name}: values")
            elif {k: str(v) for k, v in o[name].attrs.items()} != {k: str(v) for k, v in n[name].attrs.items()}:
                attr_diffs.append(f"{group} {name}: attributes")
        for key in set(n.attrs) - NOT_VALUE:
            if str(o.attrs.get(key)) != str(n.attrs[key]):
                attr_diffs.append(f"{group} global {key}")
        retired += len(set(o.attrs) - set(n.attrs) - NOT_VALUE)
    return count, value_diffs, attr_diffs, retired


# ---------------------------------------------------------------------------
# 1. Lindal's answer unchanged
# ---------------------------------------------------------------------------
status = subprocess.run(["git", "status", "--porcelain", "--", "occul_data/lindal", "forward/lindal_*"],
                        capture_output=True, text=True, check=True).stdout.strip()
WORK = HERE / "lindal_copy"
if WORK.exists():
    shutil.rmtree(WORK)
(WORK / "occul_data").mkdir(parents=True)
shutil.copytree(L, WORK / L, ignore=shutil.ignore_patterns("figures"))
shutil.copytree("data_static", WORK / "data_static")
copy = WORK / L
# Deleted first: `write_unless_unchanged` keeps what is there, and the comparison would prove nothing.
for name in ("lindal_composition.nc", "lindal_thermo.nc", "lindal_geodesy.nc", "lindal_reduction.toml"):
    (copy / name).unlink()
build_composition.build(copy / "lindal_build.toml")
build_inputs.build(copy / "lindal_build.toml")
lines, bad = [], 0
for name in ("lindal_composition.nc", "lindal_thermo.nc", "lindal_geodesy.nc"):
    count, values, attributes, retired = compare(L / name, copy / name)
    bad += len(values) + len(attributes)
    lines.append(f"{name}: {count} variables, {len(values)} differ in value, {len(attributes)} in "
                 f"attributes written, {retired} retired globals not compared {values + attributes}")
manifest_equal = (L / "lindal_reduction.toml").read_bytes() == (copy / "lindal_reduction.toml").read_bytes()
runs = run_variables = run_bad = 0
for build_file in sorted(Path("forward").glob("*/*_build.toml")):
    run = build_file.parent
    if run == RUN:
        continue
    target = WORK / run
    (target / "inputs").mkdir(parents=True)
    copied = target / build_file.name
    shutil.copyfile(build_file, copied)
    registered = run / "inputs" / f"{run.name}_composition.nc"
    build_composition.build(copied)
    count, values, attributes, _ = compare(registered, target / "inputs" / registered.name)
    runs += 1
    run_variables += count
    run_bad += len(values) + len(attributes)
    if values or attributes:
        lines.append(f"{registered}: {values + attributes}")
lines.append(f"manifest byte-equal: {manifest_equal}")
lines.append(f"{runs} runs' kind C: {run_variables} variables, {run_bad} differences")
lines.append(f"git status of occul_data/lindal and forward/lindal_*: {status or 'clean'}")
summary["lindal"] = {"runs": runs, "run_variables": run_variables, "differences": bad + run_bad,
                     "manifest_equal": manifest_equal, "git_status": status}
record(1, "Lindal's answer unchanged: no registered file changes, and the changed tools rebuild his "
          "files equal in value in a copy",
       not status and bad == 0 and run_bad == 0 and manifest_equal, "\n".join(lines))
shutil.rmtree(WORK)

# ---------------------------------------------------------------------------
# The profile, increasing pressure
# ---------------------------------------------------------------------------
product = tree(K / "kliore1980_refractivity.nc")
root = product.to_dataset(inherit=False)
raw = tree(K / "raw/kliore1980_raw.nc")
table = raw["table1"].to_dataset(inherit=False)
p_all = table["pressure_Pa"].values
tabulated = p_all != ANCHOR_PA
inserted = int(np.flatnonzero(~tabulated)[0])
p_tab_kind_n = product["inputs/thermo"].to_dataset(inherit=False)["pressure_Pa"].values
order = np.argsort(p_tab_kind_n)

# ---------------------------------------------------------------------------
# 2. Kliore's refractivity
# ---------------------------------------------------------------------------
N = np.asarray(root["refractivity"].values, dtype="float64")[order] * 1e6
N_table = table["refractivity_n_units"].values
d = N[tabulated] - N_table[tabulated]
rms_N, largest_N = float(np.sqrt(np.mean(d ** 2))), float(d[np.argmax(np.abs(d))])
inside_half = int(np.sum(np.abs(d) <= HALF_N))
summary["refractivity"] = {"rms_n": rms_N, "largest_n": largest_N, "mean_n": float(d.mean()),
                           "inside_half_step": inside_half, "inserted_n": float(N[inserted])}
record(2, "kind N's refractivity against Table 1, the inserted level excluded and reported",
       rms_N <= HALF_N,
       f"{d.size} levels: RMS {rms_N:.4f} n-units, largest {largest_N:+.4f} at "
       f"{p_all[tabulated][np.argmax(np.abs(d))] / 100:g} mbar, mean {d.mean():+.4f}, against the table's "
       f"{HALF_N} n-unit half step; {inside_half} of {d.size} inside it\n"
       f"inserted 100 mbar level: {N[inserted]:.3f} n-units, between the table's "
       f"{N_table[inserted - 1]} and {N_table[inserted + 1]}")

# ---------------------------------------------------------------------------
# 3. Kliore's radius
# ---------------------------------------------------------------------------
radius = np.asarray(root["radius_m"].values, dtype="float64")[order]
r_local = table["radius_local_m"].values
dr = radius[tabulated] - r_local[tabulated]
summary["radius"] = {"mean_km": float(dr.mean() / 1e3), "std_m": float(dr.std()),
                     "range_m": float(dr.max() - dr.min())}
record(3, "kind N's radius against Kliore's local radius: a mean offset and no spread",
       float(dr.max() - dr.min()) < 1.0,
       f"kind N minus Kliore: mean {dr.mean() / 1e3:+.3f} km, spread {dr.std():.3f} m RMS, range "
       f"{dr.max() - dr.min():.3f} m\n"
       f"SPEC_16's +{SPEC_16_OFFSET_KM} km is the delivered radius under the Sanchez-Lavega and IRIS "
       "wind, not this reduction's own geoid; the two are not the same quantity")

# ---------------------------------------------------------------------------
# 4. The closure, against Kliore's printed precision
# ---------------------------------------------------------------------------
closure = tree(RUN / "output/kliore1980_closure_profile.nc")
c_root = closure.to_dataset(inherit=False)
rec = closure["production_record"].to_dataset(inherit=False)
scalar = {**{k: v for k, v in rec.attrs.items() if np.ndim(v) == 0},
          **{k: float(rec[k]) for k in rec.variables if rec[k].ndim == 0}}
p_c, p_tab = c_root["pressure_Pa"].values, c_root["pressure_tabulated_Pa"].values
T_c, T_tab = c_root["temperature_K"].values, c_root["temperature_tabulated_K"].values
M_c = c_root["mean_molar_mass_kg_mol"].values
Phi = c_root["geopotential_m2s2"].values
h = c_root["height_above_anchor_isobar_m"].values
res = p_c / p_tab - 1.0
g_c = np.abs(np.gradient(Phi, h))
H = MOLAR_GAS_CONSTANT * T_tab / (M_c * g_c)
#: The level's printed pressure, the boundary's (it is every level's starting value), and the level's
#: printed radius over a scale height. F6 is unchanged (SPEC_18 section 5 ruling 6).
envelope = 2.0 * HALF_P_PA / p_tab + HALF_R_M / H
inside = np.abs(res) <= envelope
P2 = 10.0 ** 2.3
bins = {"above 2 mbar": p_tab < P2, "2 to 10 mbar": (p_tab >= P2) & (p_tab < 1e3),
        "10 to 100 mbar": (p_tab >= 1e3) & (p_tab < 1e4),
        "below 100 mbar, bottom row excluded": (p_tab >= 1e4) & (p_tab != p_tab.max())}
lines, table_closure = [], {}
for name, mask in bins.items():
    worst = float(np.max(np.abs(res[mask])))
    worst_T = float(np.max(np.abs((T_c - T_tab)[mask])))
    table_closure[name] = {"levels": int(mask.sum()), "max_abs_r": worst, "mean_r": float(res[mask].mean()),
                           "max_abs_dT_K": worst_T, "inside": int(inside[mask].sum()),
                           "envelope_median": float(np.median(envelope[mask]))}
    lines.append(f"{name} ({int(mask.sum())} levels): max |r| {worst:.3e}, mean {res[mask].mean():+.3e}, max "
                 f"|dT| {worst_T:.3f} K; envelope median {np.median(envelope[mask]):.2e}, {int(inside[mask].sum())} "
                 "inside")
identity = float(np.max(np.abs((T_c / T_tab - 1.0) - res)))
bottom = int(np.argmax(p_tab))
mean_10 = float(res[p_tab >= 1e3].mean())
lines.append(f"bottom row ({p_tab[bottom] / 100:g} mbar): {res[bottom]:+.3e}; mean below 10 mbar {mean_10:+.3e}")
lines.append(f"T and p residuals agree to {identity:.2e}; {int(inside.sum())} of {res.size} levels inside the "
             "envelope (the drift of check 5 is why the others are not)")
lines.append("production_record agrees: "
             + ", ".join(f"{k} {scalar[k]:.3e}" for k in sorted(scalar) if k.startswith("residual_max_abs")
                         and "temperature" not in k))
summary["closure"] = {"bins": table_closure, "bottom_r": float(res[bottom]), "mean_below_10": mean_10,
                      "identity": identity, "inside": int(inside.sum()), "levels": int(res.size),
                      "boundary_Pa": float(c_root["boundary_pressure_Pa"])}
record(4, "the closure's residuals in the bins of SPEC_03 Step 4, against Kliore's printed precision; "
          "Lindal's bounds do not apply",
       identity <= 1e-12 and float(c_root["boundary_pressure_Pa"]) == float(p_tab.min()), "\n".join(lines))

# ---------------------------------------------------------------------------
# 5. Kliore's table, hydrostatically
# ---------------------------------------------------------------------------
p = p_all[tabulated]
T = table["temperature_K"].values[tabulated]
r = r_local[tabulated]
gravity = cio.read(K / "kliore1980_gravity.nc", "gravity").load()
rotation = cio.read(K / "kliore1980_rotation.nc", "rotation").load()
constants = (float(rotation["angular_rate_rad_s"]), float(gravity["GM_m3s2"]),
             np.asarray(gravity["J"].values, dtype="float64"), validated_degrees(gravity["degree"].values),
             float(gravity["normalization_radius_m"]))
M = float(np.unique(M_c)[0])
phi = np.radians(float(root["latitude_planetocentric_deg"]))
u = float(scalar["u_at_phi_c_ms"])
at = np.full(r.shape, phi)
_, _, mag_w, psi_w = g_eff_vector(np.full(r.shape, u), r, at, *constants)
g_r0, _, _, _ = g_eff_vector(np.zeros(r.shape), r, at, *constants)
# The gravity the closure uses: the magnitude along the local vertical, dh = dr / cos psi.
g_ours = mag_w / np.cos(psi_w)


def by_density(g, mass=M):
    """`dp = rho g dr` down from the top, rho from the table's own p and T: the closure's form."""
    q = mass * p / (MOLAR_GAS_CONSTANT * T) * g
    return p[0] + np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * (-np.diff(r)))])


def by_log(g, start, mass=M):
    """`d ln p = g M dr / (R T)` from level `start`: his temperatures alone, the pressures integrated."""
    q = mass * g / (MOLAR_GAS_CONSTANT * T)
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * (-np.diff(r)))])
    return p[start] * np.exp(cumulative - cumulative[start])


def binned(pressure):
    out = {}
    q = pressure / p - 1.0
    for name, lo, hi in (("2 to 10", P2, 1e3), ("10 to 30", 1e3, 3e3), ("30 to 100", 3e3, 1e4), ("100 to 182", 1e4, 2e4)):
        out[name] = float(q[(p >= lo) & (p < hi)].mean())
    out["bottom"] = float(q[-1])
    return out


def show(name, values):
    return f"{name:44s}" + " ".join(f"{values[k]:+.2e}" for k in ("2 to 10", "10 to 30", "30 to 100", "100 to 182", "bottom"))


top10 = int(np.argmin(np.abs(p - 1000.0)))
M_22 = 2.2e-3
cases = {
    "density, ours (the closure's form)": binned(by_density(g_ours)),
    "density, ours, M = 2.2 amu": binned(by_density(g_ours, M_22)),
    "density, no wind, radial": binned(by_density(g_r0)),
    "density, no wind, radial, M = 2.2 amu": binned(by_density(g_r0, M_22)),
    "ln p from 10 mbar, ours": binned(by_log(g_ours, top10)),
    "ln p from 10 mbar, ours, M = 2.2 amu": binned(by_log(g_ours, top10, M_22)),
    "ln p from 1.8 mbar, ours": binned(by_log(g_ours, 0)),
}
# The closure's own bins, on the tabulated levels, for the comparison.
keep_c = p_tab != ANCHOR_PA
closure_binned = {}
for name, lo, hi in (("2 to 10", P2, 1e3), ("10 to 30", 1e3, 3e3), ("30 to 100", 3e3, 1e4), ("100 to 182", 1e4, 2e4)):
    closure_binned[name] = float(res[keep_c & (p_tab >= lo) & (p_tab < hi)].mean())
closure_binned["bottom"] = float(res[bottom])
agreement = max(abs(closure_binned[k] - cases["density, ours (the closure's form)"][k]) for k in closure_binned)

# The residual as a constant and an offset: r = e + c / p, from 10 mbar down.
deep = p >= 1e3
direct = by_density(g_ours) / p - 1.0
design = np.vstack([np.ones(p.size), 1.0 / p]).T
(e_fit, c_fit), *_ = np.linalg.lstsq(design[deep], direct[deep], rcond=None)
left = float(np.std(direct[deep] - design[deep] @ np.array([e_fit, c_fit])))
closure_tab = res[keep_c]
(e_close, c_close), *_ = np.linalg.lstsq(design[deep], closure_tab[deep], rcond=None)
# The product M g the table itself implies: the slope of ln p against the integral of dr / (R T).
q = 1.0 / (MOLAR_GAS_CONSTANT * T)
integral = np.concatenate([[0.0], np.cumsum(0.5 * (q[1:] + q[:-1]) * (-np.diff(r)))])
Mg_table = float(np.polyfit(integral[deep], np.log(p[deep]), 1)[0])
Mg_ours = float(M * np.mean(g_ours[deep]))
Mg_no_wind = float(M * np.mean(g_r0[deep]))
reproduces = agreement <= 1e-3
explained = (abs(e_fit - (Mg_ours / Mg_table - 1.0)) <= 2e-3
             and abs(binned(by_density(g_ours, M_22))["100 to 182"]) <= 1.5e-3)
header = f"{'mean of p / p_table - 1':44s}" + " ".join(f"{k:>9s}" for k in ("2-10", "10-30", "30-100", "100-182", "bottom"))
lines = [header, show("the closure product", closure_binned)] + [show(k, v) for k, v in cases.items()]
inside_22 = int(np.sum(np.abs(by_density(g_ours, M_22) / p - 1.0 - c_fit / p)[deep]
                       <= (2.0 * HALF_P_PA / p + HALF_R_M * M_22 * g_ours / (MOLAR_GAS_CONSTANT * T))[deep]))
inside_e = int(np.sum(np.abs(direct - e_fit - c_fit / p)[deep]
                      <= (2.0 * HALF_P_PA / p + HALF_R_M * M * g_ours / (MOLAR_GAS_CONSTANT * T))[deep]))
lines += [
    f"inside the envelope from 10 mbar down, of {int(deep.sum())}: the fit removed {inside_e}; M = 2.2 amu and the offset removed {inside_22}",
    f"the density form reproduces the closure's bins to {agreement:.2e}",
    f"r = e + c / p from 10 mbar down: direct e = {e_fit:+.2e}, c = {c_fit:+.1f} Pa, {left:.2e} RMS left; "
    f"closure e = {e_close:+.2e}, c = {c_close:+.1f} Pa",
    f"M g implied by the table {Mg_table * 1e3:.4f} (g m/mol m/s2); ours {Mg_ours * 1e3:.4f} ({Mg_ours / Mg_table - 1:+.2%}); "
    f"no wind, radial {Mg_no_wind * 1e3:.4f} ({Mg_no_wind / Mg_table - 1:+.2%})",
    f"as a gravity with M = {M * 1e3:.4f} amu: {Mg_table / M:.3f} m/s2 against ours {Mg_ours / M:.3f}; "
    f"as a mass with our gravity: {Mg_table / (Mg_ours / M) * 1e3:.4f} amu",
]
summary["hydrostatic"] = {"cases": cases, "closure": closure_binned, "agreement": agreement, "e": float(e_fit),
                          "c_Pa": float(c_fit), "rms_left": left, "e_closure": float(e_close),
                          "c_closure_Pa": float(c_close), "Mg_table": Mg_table, "Mg_ours": Mg_ours,
                          "Mg_no_wind": Mg_no_wind, "implied_g": Mg_table / M, "ours_g": Mg_ours / M,
                          "implied_mass_amu": Mg_table / (Mg_ours / M) * 1e3}
record(5, "Kliore's table integrated hydrostatically outside the model reproduces the closure's drift, "
          "and the drift is a constant in M g plus an offset from the top",
       reproduces and explained, "\n".join(lines))

# ---------------------------------------------------------------------------
# 6. Beyond the specification: section 5's rulings as built
# ---------------------------------------------------------------------------
lines, ok = [], True
for kind in ("geodesy", "gravity", "rotation", "wind"):
    count, values, _, _ = compare(L / f"lindal_{kind}.nc", K / f"kliore1980_{kind}.nc")
    ok &= not values
    lines.append(f"kliore1980_{kind}.nc against Lindal's: {count} variables, {len(values)} differ in value {values}")
reach = []
for path in sorted(K.glob("*.nc")) + [K / "raw/kliore1980_raw.nc"]:
    text = str(tree(path).attrs.get("input_hashes", ""))
    if "lindal_" in text or "lindal/" in text:
        reach.append(path.name)
ok &= not reach
lines.append(f"files whose input_hashes reach Lindal's: {reach or 'none'}")
kind_t = cio.read(K / "kliore1980_thermo.nc", "thermo").load()
phi_c = float(root["latitude_planetocentric_deg"])
label = float(kind_t.attrs["latitude_planetographic_deg"])
ok &= abs(phi_c + 9.5) <= 1e-6 and p_all.size == 129 and int(tabulated.sum()) == 128
ok &= bool(np.all(np.isnan(kind_t["temperature_uncertainty_K"].values))) and "inserted_level" not in table.variables
transcribed = np.isfinite(table["temperature_range_minus_K"].values).sum()
lines.append(f"label {label:.5f} deg returns phi_c {phi_c:.8f} deg; {p_all.size} levels, {int(tabulated.sum())} "
             f"tabulated; inserted T {table['temperature_K'].values[inserted]:.4f} K")
lines.append(f"kind T's uncertainties NaN; {int(transcribed)} rows of Kliore's ranges in the raw bundle")
summary["built"] = {"label_deg": label, "phi_c_deg": phi_c, "levels": int(p_all.size),
                    "inserted_T_K": float(table["temperature_K"].values[inserted])}
record(6, "beyond the specification: kinds D, G, R and W equal Lindal's in value, nothing reaches his "
          "files, and the label returns 9.5 S",
       ok, "\n".join(lines))

# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7.5, 5.5))
ax.fill_betweenx(p_tab / 100, -envelope, envelope, color="0.85", label="printed precision")
ax.plot(res, p_tab / 100, color="tab:blue", marker="o", markersize=2.5, linewidth=0.8, label="closure, p / p_table - 1")
ax.plot(direct, p / 100, color="tab:orange", linewidth=1.0, label="the table integrated outside the model")
ax.plot(e_fit + c_fit / p[deep], p[deep] / 100, color="k", linestyle="--", linewidth=1.0,
        label=f"e + c / p: e = {e_fit:+.2e}, c = {c_fit:+.0f} Pa")
ax.plot(by_density(g_ours, M_22) / p - 1.0, p / 100, color="tab:green", linewidth=0.8, label="the same, M = 2.2 amu")
ax.axvline(0.0, color="0.5", linewidth=0.6)
ax.set_yscale("log")
ax.invert_yaxis()
ax.set_xlim(-0.04, 0.02)
ax.set_xlabel("fractional difference from the printed pressure")
ax.set_ylabel("pressure (mbar)")
ax.set_title("Kliore et al. (1980) Table 1, 11.9 S (9.5 S planetocentric)")
ax.legend(loc="lower left", fontsize=8)
fig.tight_layout()
fig.savefig(HERE / "closure_residual.png", dpi=150)

(HERE / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")
failed = [r for r in results if not r[2]]
print(f"{len(results) - len(failed)} of {len(results)} checks pass")
sys.exit(1 if failed else 0)
