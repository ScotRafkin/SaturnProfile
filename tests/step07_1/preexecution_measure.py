"""SPEC_07 v0.2 §4, the coding agent's reading: measurements on every wind file.

A prototype of decision L2 built from `scipy.interpolate.PchipInterpolator` exactly as §1
describes it (latitude first on every pressure row, then `ln p` at the point), set beside decision
L (bilinear in latitude and `ln p`, as `lib.windfield` reads now). Nothing in `src` is used for L2.

1. The turnover loci: where, along latitude, two adjacent vertical secants of the latitude
   interpolated rows change from the same sign to opposite signs (the harmonic mean switching to
   zero), and the jump of `(du/dphi)` measured across each.
2. The predictor: the largest `|u_L2 - u_L|` on a sample of 8 points per cell along each axis,
   with its place, and at the anchor's `phi_c` on the reference level of every reduction's wind.
3. The experiments' grids near `p_s` and the stop pressure: the wind node interval there in
   `ln p` and the anchor layers it holds.

Writes `reports/step07_1/preexecution.json` and prints a summary. Run from the repository root.
"""

import json
import re
import warnings
from pathlib import Path

import numpy as np
import xarray as xr
from scipy.interpolate import PchipInterpolator

OUT = Path("reports/step07_1")
SAMPLES = 8
PHI_C_DEG = 30.80556842739218


def wind_files():
    found = sorted(set(Path("occul_data").rglob("*_wind*.nc")) | set(Path("forward").glob("*/inputs/*_wind*.nc")))
    return [p for p in found if "spacing" not in p.as_posix()]


def load(path):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        tree = xr.open_datatree(path, engine="netcdf4").load()
    tree.close()
    ds = tree.to_dataset(inherit=False)
    phi = np.radians(np.asarray(ds["latitude_planetocentric_deg"].values, dtype="float64"))
    p = np.asarray(ds["pressure_Pa"].values, dtype="float64")
    u = np.asarray(ds["u_total_ms"].values, dtype="float64")
    i, k = np.argsort(phi), np.argsort(p)
    return phi[i], np.log(p[k]), u[np.ix_(i, k)], float(ds["reference_level_pressure_Pa"])


def sample_axis(nodes, per_cell):
    t = (np.arange(per_cell) + 0.5) / per_cell
    return (nodes[:-1, None] + np.diff(nodes)[:, None] * t[None, :]).ravel()


def linear(phi, x, u, phis, xs):
    rows = np.stack([np.interp(phis, phi, u[:, k]) for k in range(x.size)], axis=1)
    return np.stack([np.interp(xs, x, rows[j]) for j in range(phis.size)])


def pchip(phi, x, u, phis, xs):
    rows = PchipInterpolator(phi, u, axis=0)(phis)
    return PchipInterpolator(x, rows, axis=1)(xs)


def loci(phi, x, u):
    """Turnover loci: latitudes where adjacent vertical secants of the PCHIP rows change from
    agreeing in sign to not, at an interior pressure node, and the jump of du/dphi there."""
    fine = np.linspace(phi[0], phi[-1], 20 * (phi.size - 1) + 1)
    rows = PchipInterpolator(phi, u, axis=0)(fine)
    secants = np.diff(rows, axis=1)
    if not np.any(secants):
        return []
    agree = secants[:, :-1] * secants[:, 1:] > 0
    found = []
    for k in range(agree.shape[1]):
        flips = np.flatnonzero(agree[1:, k] != agree[:-1, k])
        for j in flips:
            lo, hi = fine[j], fine[j + 1]
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                r = PchipInterpolator(phi, u, axis=0)(np.array([lo, mid]))
                s = np.diff(r[:, k:k + 3], axis=1)
                same_lo, same_mid = s[0, 0] * s[0, 1] > 0, s[1, 0] * s[1, 1] > 0
                lo, hi = (mid, hi) if same_lo == same_mid else (lo, mid)
            at = 0.5 * (lo + hi)
            xs = np.concatenate([np.linspace(x[k], x[k + 1], 9)[1:-1], np.linspace(x[k + 1], x[k + 2], 9)[1:-1]])
            eps, d = 1e-6, 1e-9

            def dudphi(where):
                return (pchip(phi, x, u, np.array([where + d]), xs) - pchip(phi, x, u, np.array([where - d]), xs))[0] / (2 * d)

            jump = float(np.max(np.abs(dudphi(at + eps) - dudphi(at - eps))))
            found.append({"latitude_deg": float(np.degrees(at)), "node_Pa": float(np.exp(x[k + 1])),
                          "jump_dudphi_ms_per_rad": jump})
    return found


def linear_jump(phi, x, u):
    """Decision L's largest jump of du/dphi across a latitude node, on the node rows."""
    slopes = np.diff(u, axis=0) / np.diff(phi)[:, None]
    return float(np.max(np.abs(np.diff(slopes, axis=0))))


results = {"files": {}, "phi_c": {}, "experiments": {}}
for path in wind_files():
    phi, x, u, reference = load(path)
    phis, xs = sample_axis(phi, SAMPLES), sample_axis(x, SAMPLES)
    diff = np.abs(pchip(phi, x, u, phis, xs) - linear(phi, x, u, phis, xs))
    j, k = np.unravel_index(int(np.argmax(diff)), diff.shape)
    found = loci(phi, x, u)
    entry = {
        "grid": [int(phi.size), int(x.size)],
        "largest_abs_u": float(np.max(np.abs(u))),
        "largest_diff_ms": float(diff[j, k]),
        "at_latitude_deg": float(np.degrees(phis[j])),
        "at_pressure_Pa": float(np.exp(xs[k])),
        "vertical_shear": bool(np.any(np.diff(u, axis=1))),
        "turnover_loci": len(found),
        "largest_locus_jump": max((f["jump_dudphi_ms_per_rad"] for f in found), default=0.0),
        "loci_sample": found[:5],
        "linear_largest_latitude_jump": linear_jump(phi, x, u),
    }
    results["files"][path.as_posix()] = entry
    column = int(np.argmin(np.abs(np.exp(x) - reference)))
    phic = np.radians(PHI_C_DEG)
    results["phi_c"][path.as_posix()] = float(
        PchipInterpolator(phi, u[:, column])(phic) - np.interp(phic, phi, u[:, column]))
    print(f"{path.as_posix()}: grid {entry['grid']}, max |u| {entry['largest_abs_u']:.1f} m/s, "
          f"largest |u_L2 - u_L| {entry['largest_diff_ms']:.3f} m/s at {entry['at_latitude_deg']:.2f} deg, "
          f"{entry['at_pressure_Pa']:.4g} Pa; vertical shear {entry['vertical_shear']}; turnover loci "
          f"{entry['turnover_loci']} (largest jump {entry['largest_locus_jump']:.3g}, decision L's latitude "
          f"node jump {entry['linear_largest_latitude_jump']:.3g}); at phi_c on the reference level "
          f"{results['phi_c'][path.as_posix()]:+.4f} m/s", flush=True)

# 3. The experiments' grids
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    run2 = xr.open_datatree("forward/lindal_transfer/output/lindal_transfer_profile.nc",
                            engine="netcdf4").to_dataset(inherit=False)
labels = np.sort(np.asarray(run2["pressure_label_Pa"].values, dtype="float64"))
for build in sorted(Path("forward").glob("shear_*/shear_*_build.toml")):
    text = build.read_text(encoding="utf-8")
    shear = text[text.index("[shear]"):]
    value = {key: (re.search(rf"^{key}\s*=\s*([^\n#]+)", shear, flags=re.M) or [None, None])[1]
             for key in ("case", "shear_reference_pressure_Pa", "stop_pressure_Pa", "shape", "scale")}
    grid = np.sort(np.asarray([float(v) for v in re.search(r"pressure_grid_Pa = \[(.*?)\]", text, flags=re.S)
                               .group(1).replace("\n", " ").split(",") if v.strip()]))
    rows = {}
    for name in ("shear_reference_pressure_Pa", "stop_pressure_Pa"):
        if value[name] is None:
            continue
        p = float(value[name])
        if not (grid[0] <= p <= grid[-1]):
            rows[name] = {"pressure_Pa": p, "note": "outside the wind grid"}
            continue
        on_node = bool(np.any(np.isclose(grid, p, rtol=1e-9)))
        i = int(np.clip(np.searchsorted(grid, p, side="right") - 1, 0, grid.size - 2))
        intervals = [(grid[i - 1], grid[i]), (grid[i], grid[i + 1])] if on_node and i > 0 else [(grid[i], grid[i + 1])]
        rows[name] = {"pressure_Pa": p, "on_wind_node": on_node,
                      "intervals": [{"from_Pa": float(a), "to_Pa": float(b), "d_ln_p": float(np.log(b / a)),
                                     "anchor_levels_inside": int(np.sum((labels > a) & (labels < b))),
                                     "anchor_layers_spanned": float(np.log(b / a) / np.median(
                                         np.diff(np.log(labels[(labels > a / 3) & (labels < b * 3)]))))}
                                    for a, b in intervals]}
    results["experiments"][build.parent.name] = {"case": value["case"], **rows}
    print(build.parent.name, value["case"], json.dumps(rows))

OUT.mkdir(parents=True, exist_ok=True)
(OUT / "preexecution.json").write_text(json.dumps(results, indent=1), encoding="utf-8")
