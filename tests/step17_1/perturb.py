"""The linearization test of SPEC_17 v0.2 Step 1 item 5: a known shear added to the base wind, run in
a copy of the tree.

The base is `lindal_iris_v_sl_ii_110`'s wind (hypothesis 1, which the 74.2 S and 3.0 S runs of
hypothesis 1 share). The added shear is

    Delta s(phi, p) = c A(phi) B(p),

`A` one over segment 2 (2.37 S to 26.06 S planetocentric), falling to zero over the 2 deg inside each
end as a raised cosine, so that nothing reaches the other segments; `B` one from 10 to 30 mbar, falling
to zero as a raised cosine over half a scale height outside each end (6.07 and 49.5 mbar). The wind
change is the shear integrated in `ln p` from the reference level, as the tool integrates it:
`Delta u = c A(phi) int_{ln p_ref}^{ln p} B dln p'`, zero below 49.5 mbar and `-c A (ln 3 + 0.5)`
above 6.07 mbar.

Each copy is a run directory under `reports/step17_1/perturbed/`: the namelist and inputs of a
registered run, the wind replaced and the anchor's path made to reach the repository's data from
there. Nothing under `forward/` is touched. Run from the repository root.
"""

import shutil
from pathlib import Path

import numpy as np

from casspian.forward import transfer as forward_transfer
from casspian.lib import io as cio
from casspian.tools.lindal import lindal_wind as lw

BASE = "lindal_iris_v_sl_ii_110"
WORK = Path("reports/step17_1/perturbed")
SEGMENT_DEG = (-26.06095353762054, -2.3736176374048172)
TAPER_DEG = 2.0
BAND_Pa = (1000.0, 3000.0)
TAPER_LN_P = 0.5


def raised(x):
    """0 at x <= 0, 1 at x >= 1, a raised cosine between."""
    x = np.clip(x, 0.0, 1.0)
    return 0.5 - 0.5 * np.cos(np.pi * x)


def shape_latitude(latitude_deg):
    south, north = SEGMENT_DEG
    lat = np.asarray(latitude_deg, dtype="float64")
    return raised((lat - south) / TAPER_DEG) * raised((north - lat) / TAPER_DEG)


def shape_pressure(pressure_Pa):
    x = np.log(np.asarray(pressure_Pa, dtype="float64"))
    low, high = np.log(BAND_Pa[0]), np.log(BAND_Pa[1])
    return raised((x - (low - TAPER_LN_P)) / TAPER_LN_P) * raised(((high + TAPER_LN_P) - x) / TAPER_LN_P)


def integrated_pressure(pressure_Pa):
    """`int_{ln p_ref}^{ln p} B dln p'`, on a fine grid in ln p."""
    fine = np.linspace(np.log(0.1), np.log(2.0e6), 400001)
    b = shape_pressure(np.exp(fine))
    cumulative = np.concatenate([[0.0], np.cumsum(0.5 * (b[1:] + b[:-1]) * np.diff(fine))])
    at = lambda x: np.interp(x, fine, cumulative)
    return at(np.log(np.asarray(pressure_Pa, dtype="float64"))) - at(np.log(lw.REFERENCE_PRESSURE_Pa))


def delta_s(latitude_deg, pressure_Pa, c):
    """`(pressure, latitude)`, m/s per unit ln p."""
    return c * shape_pressure(pressure_Pa)[:, None] * shape_latitude(latitude_deg)[None, :]


def delta_u(latitude_deg, pressure_Pa, c):
    """`(latitude, pressure)`, m/s, as the wind file holds it."""
    return c * shape_latitude(latitude_deg)[:, None] * integrated_pressure(pressure_Pa)[None, :]


def make_copy(registered, label, c):
    """A run directory for `registered`'s target under the base wind plus `c`'s shear; its namelist."""
    source = Path("forward") / registered
    name = f"{registered}_{label}"
    target = WORK / name
    if target.exists():
        shutil.rmtree(target)
    (target / "inputs").mkdir(parents=True)
    for kind in ("composition", "gravity", "rotation"):
        shutil.copyfile(source / "inputs" / f"{registered}_{kind}.nc", target / "inputs" / f"{name}_{kind}.nc")
    wind = cio.read(Path("forward") / BASE / "inputs" / f"{BASE}_wind.nc", "wind").load()
    change = delta_u(wind["latitude_planetocentric_deg"].values, wind["pressure_Pa"].values, c)
    wind["u_total_ms"] = wind["u_total_ms"] + change
    wind["u_shear_ms"] = wind["u_shear_ms"] + change
    wind.attrs["title"] = f"{wind.attrs.get('title', '')}; SPEC_17 Step 1 item 5, Delta s amplitude {c:g} m/s"
    cio.write(target / "inputs" / f"{name}_wind.nc", wind, "wind", created_by="tests/step17_1/perturb.py")
    text = (source / f"{registered}.toml").read_text(encoding="utf-8")
    text = text.replace(f'"{registered}', f'"{name}').replace(f"inputs/{registered}_", f"inputs/{name}_")
    text = text.replace('path   = "../../occul_data/', 'path   = "../../../../occul_data/')
    text = text.replace("figures = true", "figures = false")
    namelist = target / f"{name}.toml"
    namelist.write_bytes(text.encode("utf-8"))
    return namelist


def run_copy(registered, label, c):
    return forward_transfer.run(make_copy(registered, label, c)).product


if __name__ == "__main__":
    import sys
    import time

    t0 = time.time()
    print(run_copy(sys.argv[1], sys.argv[2], float(sys.argv[3])), f"{time.time() - t0:.0f} s")
