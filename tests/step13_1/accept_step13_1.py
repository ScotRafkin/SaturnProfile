"""Acceptance for SPEC_13 v0.3 Step 1: the gradient from the fitted temperatures, and the
Sanchez-Lavega cloud wind.

Builds the Sanchez-Lavega cloud wind from `tests/step13_1/sl_wind_build.toml` and constructs the
`lindal_iris` wind from the Ingersoll and Pollard cloud wind of `forward/lindal_iris_ii_110` under
both gradients. Writes under `reports/step13_1/`; reads the run's inputs and writes nothing there.

Run from the repository root.
"""

import json
import sys
import tomllib
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr

from casspian.tools.lindal import iris_temperatures, lindal_wind
from casspian.tools.wind import build_wind, curve as curve_module
from casspian.tools.wind.shear import _read

HERE = Path("reports/step13_1")
HERE.mkdir(parents=True, exist_ok=True)
RUN = Path("forward/lindal_iris_ii_110")
INPUTS = RUN / "inputs"
TABLE = Path("data_static/winds/vasavada_saturn_winds.txt")
results = []


def record(number, description, passed, detail):
    results.append((number, description, passed, detail))
    print(f"[{'PASS' if passed else 'FAIL'}] {number}. {description}")
    for line in str(detail).splitlines():
        print(f"        {line}")


#: The attributes a wind no longer carries (SPEC_13 v0.5 section 2a item 3, SPEC_00 v0.24): present in
#: a file built before, absent from one built now, and allowed to be.
DROPPED = {"epoch_note", "observation_level_justification", "solar_longitude_source",
           "source_latitude_convention_source", "reference_level_pressure_Pa.value_source"}


def identical(a, b, global_attributes=True, dropped=frozenset()):
    """Same values and attributes, NaN equal to NaN, apart from the build's stamps and the
    `dropped` attributes of `a` absent from `b` (a variable's as `name.attribute`). Without
    `global_attributes`, the file's global attributes, which `lib.io.write` stamps, are not compared."""
    stamps = {"created_at", "casspian_git_commit", "input_hashes", "latitude_conversion_inputs", "history"}
    stamps |= {k for k in dropped if "." not in k and k not in b.attrs}
    differ = []
    for name in sorted(set(a.variables) | set(b.variables)):
        if name not in a.variables or name not in b.variables:
            differ.append(name)
            continue
        x, y = a[name].values, b[name].values
        if not np.array_equal(x, y, equal_nan=x.dtype.kind == "f"):
            differ.append(name)
        elif (set(a[name].attrs) - {k for k in a[name].attrs if f"{name}.{k}" in dropped}) != set(b[name].attrs) or not all(
                np.array_equal(np.asarray(a[name].attrs[k]), np.asarray(b[name].attrs[k])) for k in b[name].attrs):
            differ.append(f"{name} (attributes)")
    for key in sorted((set(a.attrs) | set(b.attrs)) - stamps) if global_attributes else ():
        if not np.array_equal(np.asarray(a.attrs.get(key)), np.asarray(b.attrs.get(key))):
            differ.append(f"attribute {key}")
    return differ


def toml_value(v):
    if isinstance(v, str):
        return '"' + v.replace("\\", "/").replace('"', '\\"') + '"'
    if isinstance(v, list):
        return "[" + ", ".join(toml_value(x) for x in v) + "]"
    return repr(v)


warnings.simplefilter("ignore")
with open(RUN / f"{RUN.name}_build.toml", "rb") as handle:
    run_build = tomllib.load(handle)

# ---------------------------------------------------------------------------
# 0. The existing paths are unchanged (beyond the specification)
# ---------------------------------------------------------------------------
section = dict(run_build["wind"])
for key in build_wind.PATH_KEYS:
    if key in section:
        section[key] = str((RUN / section[key]).resolve())
section["output"] = str((HERE / "ip_wind_source_rebuilt.nc").resolve())
control = HERE / "ip_wind_build.toml"
control.write_text("[wind]\n" + "".join(f"{k} = {toml_value(v)}\n" for k, v in section.items()),
                   encoding="utf-8", newline="\n")
rebuilt = build_wind.build(control)
cloud_ip = _read(INPUTS / f"{RUN.name}_wind_source.nc", "wind")
again = _read(rebuilt, "wind")
cloud_differs = identical(cloud_ip, again, dropped=DROPPED)
gone = sorted(k for k in DROPPED if ("." in k and k.split(".")[1] in cloud_ip[k.split(".")[0]].attrs
                                     and k.split(".")[1] not in again[k.split(".")[0]].attrs)
              or ("." not in k and k in cloud_ip.attrs and k not in again.attrs))

shear = run_build["shear"]
inputs = {key: _read(INPUTS / Path(shear[key]).name, key) for key in ("composition", "gravity", "rotation")}
iris_table = np.genfromtxt("occul_data/lindal/iris_temperatures.csv", delimiter=",", names=True)
args = (cloud_ip, iris_table, inputs["composition"], inputs["gravity"], inputs["rotation"],
        float(shear["top_level_Pa"]), shear["above_top"])
slope = lindal_wind.construct(*args)
value = lindal_wind.construct(*args, gradient="value")
smooth = lindal_wind.construct(*args, gradient="value", window="smooth")
committed = _read(INPUTS / f"{RUN.name}_wind.nc", "wind")
wind_differs = identical(committed, slope.dataset, global_attributes=False)
record(0, "the existing paths are unchanged: the Ingersoll and Pollard cloud wind rebuilt by the edited "
          "tool, and the lindal_iris wind at the default gradient, against the run's files", not cloud_differs
       and not wind_differs,
       f"cloud wind: {'identical in every variable and attribute but the build stamps and the attributes SPEC_00 v0.24 drops: ' + ', '.join(gone) if not cloud_differs else cloud_differs}\n"
       f"lindal_iris wind, gradient 'slope' (the default): "
       f"{'identical in every variable and its attributes' if not wind_differs else wind_differs}")

# ---------------------------------------------------------------------------
# 1. The gradient and the shear: the slope, and the derivative with each window
# ---------------------------------------------------------------------------
lat_g = np.asarray(cloud_ip["latitude_planetographic_deg"].values, dtype="float64")
order = np.argsort(lat_g)
labels = [f"{p / 100:g} mbar" for p in slope.level_Pa]
printed = (lindal_wind.TOP_LEVEL_MBAR, 290, 730)
cases = ((slope, "slope, adaptive window (SPEC_11)", "-", "0.45"),
         (value, "derivative, adaptive window", "--", "C1"),
         (smooth, "derivative, smooth window", "-", "C0"))
fig, axes = plt.subplots(3, 3, figsize=(16, 11), sharex=True)
lines = [f"iterations: slope {slope.iterations}, derivative adaptive {value.iterations}, derivative smooth "
         f"{smooth.iterations}"]
rows = {}
for k, label in enumerate(labels):
    use = iris_table["pressure_mbar"] == printed[k]
    lat_d, T_d = iris_table["latitude_planetographic_deg"][use], iris_table["temperature_K"][use]
    fine, w_smooth = iris_temperatures.smooth_width(lat_d, 4.0, 3)
    w_adaptive = np.array([iris_temperatures._adaptive_width(lat_d - phi, 4.0, 3) for phi in fine])
    axes[0, k].plot(fine, w_adaptive, color="C1", label="adaptive (SPEC_10)")
    axes[0, k].plot(fine, w_smooth, color="C0", label="smooth")
    axes[0, k].set_ylabel("window FWHM (deg)")
    axes[0, k].set_title(label)
    for result, name, style, colour in cases:
        axes[1, k].plot(lat_g[order], result.gradient_K_per_rad[k][order] * np.pi / 180.0, style, color=colour,
                        linewidth=1.0, label=name)
        axes[2, k].plot(lat_g[order], result.shear_ms[k][order], style, color=colour, linewidth=1.0, label=name)
    axes[1, k].set_ylabel("dT/dphi (K per planetocentric deg)")
    axes[2, k].set_ylabel("shear du/dln p (m/s)")
    axes[2, k].set_xlabel("planetographic latitude (deg)")
    for ax in axes[:, k]:
        ax.grid(alpha=0.3)
    for ax in axes[1:, k]:
        ax.axhline(0, color="0.7", linewidth=0.8)
    rms = {}
    for window in ("adaptive", "smooth"):
        at_points = iris_temperatures.fit(lat_d, T_d, grid=lat_d, window=window).value_K
        rms[window] = float(np.sqrt(np.mean((T_d - at_points) ** 2)))
    entry = {"fit_rms_adaptive_K": rms["adaptive"], "fit_rms_smooth_K": rms["smooth"],
             "smooth_minus_adaptive_width_min_deg": float((w_smooth - w_adaptive).min()),
             "width_adaptive_max_deg": float(w_adaptive.max()), "width_smooth_max_deg": float(w_smooth.max())}
    for result, name in ((value, "adaptive"), (smooth, "smooth")):
        both = np.isfinite(slope.gradient_K_per_rad[k]) & np.isfinite(result.gradient_K_per_rad[k])
        d_grad = (result.gradient_K_per_rad[k] - slope.gradient_K_per_rad[k])[both] * np.pi / 180.0
        d_shear = result.shear_ms[k] - slope.shear_ms[k]
        worst = int(np.nanargmax(np.abs(d_shear)))
        entry[name] = {"gradient_rms_K_per_deg": float(np.sqrt(np.mean(d_grad ** 2))),
                       "gradient_max_K_per_deg": float(np.max(np.abs(d_grad))),
                       "shear_rms_ms": float(np.sqrt(np.mean(d_shear ** 2))),
                       "shear_max_ms": float(d_shear[worst]), "shear_max_at_deg": float(lat_g[worst])}
    rows[label] = entry
    lines.append(f"{label}: the fit's RMS about its {int(use.sum())} points, adaptive {rms['adaptive']:.3f} K, smooth "
                 f"{rms['smooth']:.3f} K; window width adaptive up to {w_adaptive.max():.2f} deg, smooth up to "
                 f"{w_smooth.max():.2f} deg, smooth minus adaptive at least {(w_smooth - w_adaptive).min():.2e} deg")
    for name in ("adaptive", "smooth"):
        e = entry[name]
        lines.append(f"    derivative ({name} window) - slope: dT/dphi RMS {e['gradient_rms_K_per_deg']:.4f}, largest "
                     f"{e['gradient_max_K_per_deg']:.4f} K per planetocentric deg; shear RMS {e['shear_rms_ms']:.2f} m/s, "
                     f"largest {e['shear_max_ms']:+.2f} m/s at {e['shear_max_at_deg']:.1f} deg planetographic")
axes[0, 0].legend(fontsize=8)
axes[1, 0].legend(fontsize=8)
fig.suptitle("The fit's window, the IRIS gradient and the shear at the three levels: the fit's local slope, and the "
             "derivative of its value with the adaptive and the smooth window (Ingersoll and Pollard cloud wind)",
             fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "shear_slope_value.png", dpi=120)
plt.close(fig)
record(1, "the window, the gradient and the shear at the three levels from the slope (SPEC_11) and from the "
          "derivative with the adaptive and the smooth window, against latitude, with the fit's RMS about the "
          "points for each window, for the author", True, "\n".join(lines) + f"\nfigure {HERE / 'shear_slope_value.png'}")

# ---------------------------------------------------------------------------
# 2. The two cloud winds
# ---------------------------------------------------------------------------
sl_path = build_wind.build("tests/step13_1/sl_wind_build.toml")
cloud_sl = _read(sl_path, "wind")
table = curve_module.read_table(TABLE)


def reference(wind):
    p = np.asarray(wind["pressure_Pa"].values, dtype="float64")
    column = int(np.flatnonzero(p == float(wind["reference_level_pressure_Pa"]))[0])
    g = np.asarray(wind["latitude_planetographic_deg"].values, dtype="float64")
    o = np.argsort(g)
    return (g[o], np.asarray(wind["u_total_ms"].values)[o, column],
            np.asarray(wind["u_total_uncertainty_ms"].values)[o, column],
            np.asarray(wind["value_provenance"].values)[o, column])


g_ip, u_ip, s_ip, f_ip = reference(cloud_ip)
g_sl, u_sl, s_sl, f_sl = reference(cloud_sl)
fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
ok = np.isfinite(s_ip)
axes[0].fill_between(g_ip[ok], (u_ip - s_ip)[ok], (u_ip + s_ip)[ok], color="C0", alpha=0.15, linewidth=0,
                     label="Ingersoll and Pollard +- 1 sigma (Smith points about the curve)")
axes[0].plot(g_ip, u_ip, color="C0", linewidth=1.2, label="Ingersoll and Pollard (1982), the current cloud wind")
axes[0].plot(g_sl, u_sl, color="C3", linewidth=1.2, label="Sanchez-Lavega et al. (2000), as built")
axes[0].errorbar(table.latitude_deg, table.u_ms, yerr=table.u_rms_ms, fmt=".", color="k", markersize=2.5,
                 elinewidth=0.5, label="Table II rows +- u_rms")
ring = (g_sl > -10.7) & (g_sl < -2.1)
axes[0].axvspan(-10.7, -2.1, color="0.85", zorder=0, label="ring gap, -10.7 to -2.1 (uncertainty NaN)")
axes[0].set_ylabel("cloud-top zonal wind (m/s)")
axes[0].legend(fontsize=8, loc="upper right")
common = np.interp(g_sl, g_ip, u_ip)
axes[1].plot(g_sl, u_sl - common, color="k", linewidth=1.0)
axes[1].axhline(0, color="0.7", linewidth=0.8)
axes[1].axvspan(-10.7, -2.1, color="0.85", zorder=0)
axes[1].set_ylabel("Sanchez-Lavega - Ingersoll and Pollard (m/s)")
axes[1].set_xlabel("planetographic latitude (deg)")
for ax in axes:
    ax.grid(alpha=0.3)
    ax.set_xlim(-90, 90)
fig.suptitle("The two cloud winds against planetographic latitude", fontsize=10)
fig.tight_layout()
fig.savefig(HERE / "cloud_winds.png", dpi=120)
plt.close(fig)
meanings = str(cloud_sl["value_provenance"].attrs["flag_meanings"]).split()
counts = ", ".join(f"{meanings[c]} {int((f_sl == c).sum())}" for c in np.unique(f_sl))
rows_g = (g_sl >= table.latitude_deg.min()) & (g_sl <= table.latitude_deg.max())
diff = (u_sl - common)[rows_g & ~ring]
curve_at_rows = curve_module.TableCurve(table, gap_rule="pchip_bridge", polar_rule="pchip_to_zero",
                                        ring_gap_deg=(-10.7, -2.1))(table.latitude_deg)
peak = int(np.argmax(u_sl))
record(2, "the two cloud winds against latitude, with the table's points and u_rms, for the author", True,
       f"{sl_path}: {table.latitude_deg.size} rows read, {table.dropped} dropped as not data; through every row "
       f"to {np.max(np.abs(curve_at_rows - table.u_ms)):.1e} m/s\n"
       f"grid {g_sl.size} latitudes, provenance at the reference level: {counts}; uncertainty NaN at "
       f"{int((~np.isfinite(s_sl)).sum())} (ring gap {int((ring & ~np.isfinite(s_sl)).sum())} of {int(ring.sum())}, "
       f"the rest poleward of the rows)\n"
       f"in the ring gap the wind runs {u_sl[ring].max():.1f} to {u_sl[ring].min():.1f} m/s; the current wind there "
       f"{common[ring].max():.1f} to {common[ring].min():.1f}\n"
       f"Sanchez-Lavega - Ingersoll and Pollard, within the rows and outside the ring gap: mean {diff.mean():+.1f}, "
       f"RMS {np.sqrt(np.mean(diff ** 2)):.1f} m/s; peak {u_sl[peak]:.1f} m/s at {g_sl[peak]:.1f} deg "
       f"(current {u_ip.max():.1f} at {g_ip[np.argmax(u_ip)]:.1f})\n"
       f"figure {HERE / 'cloud_winds.png'}")

for d in (cloud_ip, again, committed, cloud_sl, *inputs.values()):
    d.close()
(HERE / "summary.json").write_text(json.dumps({"shear": rows}, indent=2), encoding="utf-8")
failed = [r for r in results if not r[2]]
line = f"{len(results) - len(failed)} of {len(results)} checks pass"
print(line)
(HERE / "output.txt").write_text(
    "\n".join(f"[{'PASS' if ok else 'FAIL'}] {n}. {d}\n        " + str(detail).replace("\n", "\n        ")
              for n, d, ok, detail in results) + f"\n\n{line}\n", encoding="utf-8")
sys.exit(1 if failed else 0)
