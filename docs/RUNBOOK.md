# CASSPIAN runbook: from a clean clone to the transfer product

Version 0.2, 28 September 2026. Author of record: S. Rafkin. Status: sections 2 to 5 run by the
author on a clean clone of `e71d59e` on a CentOS 7 machine (glibc 2.17, bash 4.2, no C++
compiler, system Python 3.7.4, Anaconda on the path), every value reproduced; section 6, the
regression, waits on the acceptance scripts being committed under `tests/` (SPEC_00 v0.22). This is the document
a person follows to install the package and reproduce the test case without an agent.

The repository holds source, control files, static data, specifications and reports, and no
netCDF. Every product, from the raw bundle to the transfer profile and its figures, is rebuilt from
the control files by the console tools in the order below. That is the design (SPEC_00 §2), and it
is what makes a clone a test: if the rebuild reproduces the registered values, the clone is the same
code on the same data.

## 1. What the machine needs

Python 3.11 or later (`tomllib` is used throughout; a system Python 3.7 cannot run the code and is
not to be replaced, since the operating system depends on it), `git` on the path (every product
records the commit it was written at, and `refrac` and `forward` refuse an input written from a
dirty tree, so the rebuild is done on a clean checkout), and `sha256sum` for the regression script.
The declared dependencies are `numpy`, `scipy`, `xarray` and `netCDF4`; `matplotlib` is imported by
the plotting and reduction code and was not declared at `e71d59e` (fixed in the commit that adds
the acceptance scripts). No compiler is needed provided every package is installed from a
prebuilt wheel, which section 2 forces.

## 2. Clone and install

The interpreter comes from `uv`, which installs a prebuilt Python in the user's home directory
beside the system one. conda was tried on the author's machine and its dependency solve did not
finish in reasonable time; `uv` took under a minute.

```
curl -LsSf https://astral.sh/uv/install.sh | sh     # once per machine; installs to ~/.local/bin, then open a new shell
git clone git@github.com:ScotRafkin/SaturnProfile.git
cd SaturnProfile
git status --porcelain                              # prints nothing
uv venv --python 3.12 --seed .venv                  # --seed puts pip in the environment
source .venv/bin/activate
which python pip; python --version                  # both under .venv/bin/, Python 3.12
pip install --only-binary=:all: numpy scipy xarray netCDF4 matplotlib
pip install -e . --no-deps
casspian-forward --help                             # prints its usage
python -c "import numpy, scipy, netCDF4, xarray, matplotlib; print(numpy.__version__, scipy.__version__)"
```

Three things the author's machine taught. Without `--seed`, `uv venv` makes an environment with no
pip in it, and `pip` then resolves to another installation's (here Anaconda's or the system's),
whose old version refuses an editable install from `pyproject.toml` alone; `which pip` catches it.
The newest numpy and contourpy publish wheels only for glibc 2.28 and later, so on glibc 2.17 pip
falls back to compiling from source and fails for want of a C++ compiler; `--only-binary=:all:`
makes it choose the newest versions that still ship wheels for the platform, which on the author's
machine were numpy 2.2.6 and scipy 1.16.3 under Python 3.12.14. `--no-deps` on the package install
keeps pip from re-resolving numpy from source. A newer machine (glibc 2.28 or later, or a C++
compiler present) can use `pip install -e .` alone.

Every command below is run from the repository root with the environment active; a new shell
needs `source .venv/bin/activate` again.

## 3. Rebuild the Lindal reduction (the anchor)

The chain writes into `occul_data/lindal/`; the order is the order the files depend on each other
(SPEC_01, SPEC_03 §"Rebuild and sweep"). Each tool takes the build control file and reads its own
section.

```
casspian-lindal-raw                                                  # raw/lindal_raw.nc from Table I and the scalars
casspian-gravity-file      occul_data/lindal/lindal_build.toml       # lindal_gravity.nc, Null et al. 1981
casspian-rotation-file     occul_data/lindal/lindal_build.toml       # lindal_rotation.nc, System III
casspian-wind-from-curve   occul_data/lindal/lindal_build.toml       # lindal_wind.nc, Ingersoll and Pollard 1982
casspian-composition-lindal occul_data/lindal/lindal_build.toml      # lindal_composition.nc
casspian-lindal-inputs     occul_data/lindal/lindal_build.toml       # lindal_thermo.nc, lindal_geodesy.nc, lindal_reduction.toml
casspian-refrac            occul_data/lindal/lindal_reduction.toml   # lindal_refractivity.nc, kind N, and figures/
```

```
python -c "import netCDF4 as nc; d=nc.Dataset('occul_data/lindal/lindal_refractivity.nc'); N=d['refractivity'][:]; s=d['refractivity_uncertainty'][:]; print(len(N), (s/N).min(), (s/N).max(), float(d['anchor_isobar_radius_m'][:]), float(d['latitude_planetocentric_deg'][:]))"
```

Measured on the author's Linux clone: `66 0.023318454671386792 0.0233184546713868
58516188.288884744 30.80556842739218`.

What says it worked: `occul_data/lindal/lindal_refractivity.nc` validates as kind N with 66 levels;
`refractivity_uncertainty / refractivity` is 2.331845e-2 at every level; the anchor isobar is
1.0e4 Pa at planetocentric latitude 30.8056°; `anchor_isobar_radius_m` is 58,516,188 m; the
figures F1 to F4 are written under `occul_data/lindal/figures/`. The registered SHA-256 at
`e71d59e` is `920304440ee4554648c3aa6321b713741ee26eb0eb1ff223f5a81a8f9b20fd37`; a rebuild on
another machine is expected to differ in its bytes (the files record `created_at`, and the last
digits of floating point can differ between platforms), so the hash is recorded for the report and
the values above are the check.

## 4. Rebuild the closure run and compare with the registered values

```
casspian-run-inputs forward/lindal_closure/lindal_closure_build.toml   # inputs/: composition, gravity, rotation, wind
casspian-forward    forward/lindal_closure/lindal_closure.toml         # output/lindal_closure_profile.nc, F5 and F6
```

What says it worked: the closure product's produced pressure at the gauge level (index 29) is
9998.46545058 Pa, its temperature there 83.38720186 K, and the pressure at the top level
19.95262315 Pa; the closure residual is 9.74e-4 above 2 mbar and −3.28e-3 at the bottom row, as
SPEC_03 records. Registered SHA-256 at `e71d59e`:
`fcc2e2c4ae71554537d8aa53f07b247724b718a371abab73de31939bae9ae6af`, with the same caveat.

A one-line check of the three values:

```
python -c "import netCDF4 as nc; d=nc.Dataset('forward/lindal_closure/output/lindal_closure_profile.nc'); p=d['pressure_Pa'][:]; T=d['temperature_K'][:]; print(p[29], T[29], p[0])"
```

Reproduced on the author's Linux clone to every printed digit.

## 5. Run the transfer to 10° N (the test case)

```
casspian-run-inputs forward/lindal_transfer/lindal_transfer_build.toml
casspian-forward    forward/lindal_transfer/lindal_transfer.toml       # 91 s on the author's Linux machine, 70 to 100 s on Windows
casspian-plots      forward/lindal_transfer/output/lindal_transfer_profile.nc --out /tmp/f   # the same figures by hand
```

What says it worked (SPEC_04 §7, v0.22): the product `forward/lindal_transfer/output/
lindal_transfer_profile.nc` validates as kind `profile` in transfer mode with the groups `anchors/
lindal`, `anchors/lindal/transfer`, `reference_surface`, `isobars`, `estimate`, `inputs`,
`namelist` and `transfer_record`; `altitude_m` is 411,135 m at the top level, 98,187 m at the
gauge level and −15,290 m at the bottom (to 5 m); `reference_surface_radius_m` is 60,128,613 m
(to 1 m); the delivered temperature lies below the anchor's by 1.1022e-2 at the top, 1.0884e-2 at
the gauge and 1.0833e-2 at the bottom (to 1e-4); the largest `|pressure_identity_residual|` is
5.78e-7 (bounded at 1e-6); `transfer_record` carries `passes = 2`, `shear_kernel_largest_abs_s_
over_g_per_rad` 0.0976 and `isobar_shift_largest_lindal_m2s2` −31,241. Figures F5 to F9 are in
`forward/lindal_transfer/output/figures/`, and the copies `casspian-plots` writes are identical
apart from the footer.

```
python -c "import netCDF4 as nc; d=nc.Dataset('forward/lindal_transfer/output/lindal_transfer_profile.nc'); z=d['altitude_m'][:]; r=d['pressure_identity_residual'][:]; print(z[0], z[29], z[-1], float(d['reference_surface_radius_m'][:]), abs(r).max())"
```

Measured on the author's Linux clone: `411134.8829284598 98186.67066940351 -15290.278049858825
60128612.966417134 5.776665666923364e-07`, the Windows product's values to the last digit of the
double; wall time `real 1m31.1s`.

## 6. The regression

Once the acceptance scripts are in the repository under `tests/` (SPEC_00 v0.22; the incoming
agent's first commit), the regression is the proof that the clone is the same code:

```
git pull                                     # on the clone, after that commit is pushed
nohup bash tests/run_regression.sh > regression_console.txt 2>&1 &
                                             # about seven hours on Windows, four and a half of them step04_4;
                                             # nohup so that closing the terminal does not stop it
```

Every row must read `N of N checks pass` at its reference count; the rows are written to
`reports/regression/regression.txt` with each suite's log beside it, and the reference counts are
the ones in `docs/specs/STATE.md`. The script sets the
registered products aside and restores them after each suite, so it is run on the products of
sections 3 and 4. A suite that needs a fixture the clone does not have (a `before/` baseline, a
`-dirty` copy) is a finding for the runbook, not a failure of the code.

## 7. What to report

`reports/REPORT_05_step0.md` (the clone as SPEC_05's Step 0) records: the machine, the
interpreter and the package versions of section 2; the values of sections 3 to 5 as measured
here; the SHA-256 of the rebuilt kind N and closure product beside the registered ones, which are
not expected to agree since every file records `created_at`, with the largest numerical difference
across each file's variables in their place; the regression rows; and the timings beside the
Windows ones.

## Revision history

| Version | Date | Change |
|---|---|---|
| 0.2 | 2026-09-29 | Sections 2 to 5 run by the author on a CentOS 7 clone: the `uv` route with `--seed`, `--only-binary=:all:` and `--no-deps` for glibc 2.17 without a compiler, the measured versions, values and wall time; section 6 with `nohup` and the regression driver at `tests/run_regression.sh` |
| 0.1 | 2026-09-28 | Written from the repository at `e71d59e` |
