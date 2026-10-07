"""SPEC_08 pre-execution measurements (REPORT_08_preexecution). Reads the tree; writes only under
reports/step08_1/preexecution/. The SPEC_08 schema is simulated in memory: the four parts made
optional in the kind W table and the two wind checks removed from `validate`; `lib.control` keeps
its own references to them, so the loaders run as they are on disk."""
import dataclasses, os, sys, time
from pathlib import Path
from types import MappingProxyType
import numpy as np
import xarray as xr

REPO = Path(__file__).resolve().parents[2]
SCR = REPO / "reports/step08_1/preexecution"
SCR.mkdir(parents=True, exist_ok=True)
sys.path.insert(0, str(REPO / "src"))
os.chdir(REPO)
from casspian.lib import io as cio, schema as sch, control as ctl
from casspian.tools.wind import shear

WINDS = ["occul_data/lindal/lindal_wind.nc",
         "forward/lindal_closure/inputs/lindal_closure_wind.nc",
         "forward/lindal_closure/inputs/lindal_closure_wind_source.nc",
         "forward/lindal_transfer/inputs/lindal_transfer_wind.nc",
         "forward/lindal_transfer/inputs/lindal_transfer_wind_source.nc",
         "tests/step03_3/fixtures/before/lindal_wind.nc",
         "tests/step04_0/fixtures/swept/forward/lindal_closure/inputs/lindal_closure_wind.nc",
         "tests/step04_0/fixtures/swept/occul_data/lindal/lindal_wind.nc"]

def short(e): return f"{type(e).__name__}: {str(e)[:230]}"

print("== 1. registered and fixture wind files")
for w in WINDS:
    try:
        d = cio.read(w, "wind").load()
    except Exception as e:
        print(f"{w}: READ REFUSED {short(e)}"); continue
    lat = d["latitude_planetocentric_deg"].values; p = d["pressure_Pa"].values
    parts = [n for n in ("u_total_ms", "u_reference_ms", "u_shear_ms", "reference_level_pressure_Pa") if n in d]
    ref = float(d["reference_level_pressure_Pa"]); node = bool(np.any(p == ref))
    u = d["u_total_ms"].values
    poles = [float(np.max(np.abs(u[np.abs(lat - s) <= 1e-9]))) for s in (-90, 90)]
    s = float(np.nanmax(np.abs(u - (d["u_reference_ms"].values[:, None] + d["u_shear_ms"].values))))
    extra = sorted(set(d.data_vars) - {"u_total_ms", "u_reference_ms", "u_shear_ms", "u_total_uncertainty_ms", "value_provenance", "reference_level_pressure_Pa"})
    print(f"{w}: {lat.size}x{p.size}, parts {parts}, ref {ref:g} Pa node {node}, |u| at poles {poles}, "
          f"sum departure {s:.2e}, dims u_total {d['u_total_ms'].dims}, extra vars {extra}")
    d.close()

print("\n== 2. step1 check 3 message, now")
edited = REPO / "reports/step1/products/acceptance_raw_edited_to_wind.nc"
try: cio.read(edited, "wind"); print("not refused")
except Exception as e: print(short(e))

# --- simulate SPEC_08 deliverable 1 in memory -----------------------------------------
new_vars = tuple(dataclasses.replace(v, required=False) if v.name in
                 ("u_total_ms", "u_reference_ms", "u_shear_ms", "reference_level_pressure_Pa") else v
                 for v in sch._WIND.variables)
sch.KINDS["wind"] = dataclasses.replace(sch._WIND, variables=new_vars)
sch.check_wind_components = lambda ds, where: None   # validate looks these up at call time
sch.check_wind_poles = lambda ds, where: None          # lib.control keeps its own references

print("\n== 3. step1 check 3 message, simulated SPEC_08 schema")
try: cio.read(edited, "wind"); print("not refused")
except Exception as e: print(short(e))

print("\n== 4. acceptance 2 files under the simulated schema")
base = cio.read(WINDS[0], "wind").load()
def subgrid(ds, lat_keep, p_keep):
    return ds.isel(latitude_planetocentric=lat_keep, pressure=p_keep)
lat = base["latitude_planetocentric_deg"].values; p = base["pressure_Pa"].values
inner = np.flatnonzero(np.abs(lat) < 85)
k_ref = int(np.flatnonzero(p == float(base["reference_level_pressure_Pa"]))[0])
print(f"base pressure nodes {p.tolist()}")
cases = {
  "shear_only": (subgrid(base, inner, [max(k_ref-1,0), k_ref, k_ref+1]).drop_vars(["u_total_ms", "u_total_uncertainty_ms", "u_reference_ms"])),
  "total_only": (subgrid(base, inner, [k_ref]).drop_vars(["u_reference_ms", "u_shear_ms", "reference_level_pressure_Pa"])),
  "reference_only": (subgrid(base, inner, [k_ref]).drop_vars(["u_total_ms", "u_total_uncertainty_ms", "u_shear_ms"])),
}
for name, ds in cases.items():
    try:
        out = cio.write(SCR / f"{name}.nc", ds, "wind", created_by="tests/step08_1/preexecution_measure.py")
        back = cio.read(out, "wind").load()
        eq = all(np.array_equal(back[v].values, ds[v].values, equal_nan=True) for v in ds.variables)
        print(f"{name}: written and read back, vars {sorted(ds.data_vars)}, array-equal {eq}")
        back.close()
    except Exception as e:
        print(f"{name}: {short(e)}")
# an extra companion for the reference part (SPEC_11's data carry uncertainties)
ds = cases["reference_only"].copy()
ds["u_reference_uncertainty_ms"] = (("latitude_planetocentric",), np.full(inner.size, 5.0),
    {"units": "m s-1", "long_name": "test", "provenance": "measured", "uncertainty_kind": "1sigma"})
try:
    cio.write(SCR / "reference_with_companion.nc", ds, "wind", created_by="tests/step08_1/preexecution_measure.py"); print("reference_only + u_reference_uncertainty_ms: written")
except Exception as e: print("reference_only + companion:", short(e))

print("\n== 5. the loaders on stripped copies (simulated schema, loaders as on disk)")
def stripped(src, drop, dest):
    d = xr.open_dataset(src).load(); d.close()
    if drop == "pole":
        d["u_total_ms"].values[np.abs(d["latitude_planetocentric_deg"].values - 90) <= 1e-9, :] = 1e-6
    else:
        d = d.drop_vars(drop)
    enc = {n: {"_FillValue": np.nan} for n, v in d.data_vars.items() if np.issubdtype(v.dtype, np.floating)}
    d.to_netcdf(dest, engine="netcdf4", encoding=enc)  # attributes kept, so the commit stays clean
    return dest

manifest = ctl.read_reduction_manifest("occul_data/lindal/lindal_reduction.toml")
closure = ctl.read_run_namelist("forward/lindal_closure/lindal_closure.toml")
transfer = ctl.read_run_namelist("forward/lindal_transfer/lindal_transfer.toml")
LOADERS = [
  ("load_reduction_inputs", manifest, ctl.load_reduction_inputs),
  ("load_run_inputs (closure)", closure, ctl.load_run_inputs),
  ("_load_transfer_inputs", transfer, ctl.load_run_inputs),
]
for label, obj, fn in LOADERS:
    t = time.perf_counter()
    try: fn(obj); print(f"{label}: registered inputs admitted in {time.perf_counter()-t:.1f} s")
    except Exception as e: print(f"{label}: registered REFUSED {short(e)}")
    src = obj.inputs["wind"]
    for drop in (["u_total_ms", "u_total_uncertainty_ms"], ["reference_level_pressure_Pa"], "pole"):
        tag = drop if drop == "pole" else drop[0]
        dest = stripped(src, drop, SCR / f"{label.split()[0]}_{tag}.nc")
        replaced = dataclasses.replace(obj, inputs=MappingProxyType({**obj.inputs, "wind": dest}))
        try:
            fn(replaced); print(f"  without {tag}: ADMITTED")
        except Exception as e:
            print(f"  without {tag}: {short(e)}")

print("\n== 6. the shear tool on reduced inputs (simulated schema, tool as on disk)")
src = cio.read("forward/lindal_transfer/inputs/lindal_transfer_wind_source.nc", "wind").load()
P5 = dict(shear_reference_pressure_Pa=1e5, shape="linear_ln_p", stop_pressure_Pa=700.0, stop_fraction=0.0)
for label, drop in (("no u_reference_ms", ["u_reference_ms"]),
                    ("no u_reference_ms, no u_shear_ms", ["u_reference_ms", "u_shear_ms"]),
                    ("no u_shear_ms", ["u_shear_ms"]),
                    ("no reference_level_pressure_Pa", ["reference_level_pressure_Pa"]),
                    ("no u_total_ms", ["u_total_ms", "u_total_uncertainty_ms"])):
    try:
        out = shear.construct(src.drop_vars(drop), "decay_above", P5)
        print(f"{label}: constructed, vars {sorted(out.data_vars)}")
    except Exception as e:
        print(f"{label}: {short(e)}")
