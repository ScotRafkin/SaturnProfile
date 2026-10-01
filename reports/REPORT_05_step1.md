# REPORT 05, Step 1. The shear tool, and the registered transfer product

Coding agent, 30 September 2026. Specification: `docs/specs/SPEC_05_Shear_Experiments.md` v0.4,
sections 1.3 to 1.7, 1a and 2, with the section 7 rulings. Working tree on `main` at `7276f68`;
nothing of this step is committed.

**Result.** `tests/step05_1/accept_step05_1.py` passes 11 of 11 (checks 1 to 9, 9a, 10) in 20 s.
Every closed form is met at 0 to 1.459e-16 relative against the 1e-12 bound. The input file is read
and never written: its SHA-256 is the same before and after the run. §1a is done and committed
(`7276f68`), and it required one change to the regression driver, which is part of this review at
the author's direction (section 3).

## 1. What is built

- `src/casspian/tools/wind/shear.py`, in the three layers of §1.5.
  - The parts: `u_at_pressure` reads `u_s(phi)` through `lib.windfield.WindField` at every latitude
    node (ruling 9), turning `WindField`'s refusal outside the grid into the named `p_s` refusal;
    `position_ln_p` and `position_p` are the two shapes of `x(p)`; `ramp` is `F(x)`.
  - The cases: `identity`, `uniform`, `decay_above`, `increase_below`, in the table `CASES`, each
    with the parameters it requires and accepts. The two ramp cases share one function and differ
    only in the side of `p_s` their stop pressure must lie on.
  - `assemble` is the one function every case except `identity` passes through: it sets `u_total`,
    forms `u_shear = u_total - u_reference`, sets `value_provenance` to 2 in every cell, replaces
    `vertical_structure` by the case name, and, only when `uncertainty_ms` is given, replaces the
    uncertainty's values and its `long_name` and `uncertainty_method` by one statement. Nothing is
    added (§1.7).
  - `construct(source, case, parameters)` touches no file. `build(control_path, section="shear")`
    reads the section and the source, calls `construct`, sets `input_hashes` (the source and the
    control file) and writes with `lib.io.write`.
- `pyproject.toml`: the console entry `casspian-wind-shear`, taking the control file and
  `--section` (default `shear`).
- The `[shear]` keys are `source`, `output`, `case`, the case parameters of §1.4 and the optional
  `uncertainty_ms`. Any other key, `role` included, is refused by `lib.control.load_section` as
  unknown, which is the existing rule for every section.
- `tests/step05_1/accept_step05_1.py`, the acceptance suite, and its row in the driver table,
  `step05_1/accept_step05_1 11`.
- `tests/run_regression.sh`: the registered transfer product set aside and restored (section 3).

## 2. The acceptance

The input is `forward/lindal_transfer/inputs/lindal_transfer_wind.nc`, 361 latitudes by 61
pressures, 1 Pa to 1 MPa, `altitude independent`, SHA-256 `266fdc38...` before and after. Each case
is written through `build` from a control file under `reports/step05_1/cases/`; check 5's file is
written by the console entry.

| Check | What | Measured | Bound | Result |
|---|---|---|---|---|
| 1 | `identity`: every variable array-equal; closure content check | 12 of 12 variables equal; content check identical | exact | pass |
| 2 | `uniform`, c = 1, p_s = 1e5 | `u_total` array-equal; source data array-equal | exact | pass |
| 3 | `uniform`, c = 0 | `u_total` zero at all 22,021 nodes; `u_shear = -u_reference` exactly; reads back | exact | pass |
| 4 | `uniform`, c = 1, p_s = 1.1e5 between nodes, input sheared in pressure | 0.000e+00 relative | 1e-12 | pass |
| 5 | `decay_above`, `linear_ln_p`, p_stop 700, f 0 (console entry) | 0.000e+00 relative; `u_s` exactly at the 11 levels at and below 1e5; zero exactly at the 29 at and above 700 | 1e-12 | pass |
| 6 | `decay_above`, `linear_p`, p_stop 0, f 0 | 1.459e-16 relative on the 50 levels above p_s | 1e-12 | pass |
| 7 | `increase_below`, `linear_ln_p`, p_stop 1e6, f 1.5 | 0.000e+00 relative; exactly `1.5 u_s` at 1e6 | 1e-12 | pass |
| 8 | outputs of 2 to 7: sum identity, poles, provenance, source data, attributes | all six: poles exactly zero, provenance 2, flags kept, source data equal, no unexpected attribute | exact | pass |
| 9 | `construct` in memory against the file `build` writes (check 5's case) | 12 of 12 variables equal | exact | pass |
| 9a | uncertainty carried; `uncertainty_ms = 25` | array-equal in all six; 25 at every node | exact | pass |
| 10 | the named refusals | all ten cases refused, nothing written | | pass |

Notes on what the numbers mean.

- **Check 2 answers ruling 9's question.** 1e5 Pa is node 50 of the wind grid, and `u_total` came
  back array-equal, so `WindField` returns a node value exactly. That follows from its weights: at a
  node the weight is exactly 0 or 1, and at the last node the clipped interval gives a weight of
  exactly 1.
- **Checks 4, 5 and 7 measure zero, not round-off.** The script's closed forms use the same
  arithmetic as the tool (the ln p interpolant in the form `a (1 - w) + b w`; `F` from the same
  `x`), so they agree bit for bit. They are still independent statements of the formula, written in
  the script from §1.3 and §1.4 and not imported from the tool. Check 6's closed form, `u_s p / p_s`,
  is algebraically but not arithmetically the tool's `1 - x`, and measures 1.459e-16.
- **Check 8's attribute rule** (v0.4 wording): the globals that differ are exactly the four
  `lib.io.write` stamps (`created_by`, `created_at`, `casspian_git_commit`, `history`) and the two
  the tool replaces (`vertical_structure`, `input_hashes`); no variable attribute differs. In check 4
  `casspian_git_commit` does not differ, because the synthetic input was written from the same
  working tree.
- **Check 10**: unknown case, unknown shape, missing parameter, a parameter the case does not use,
  `uncertainty_ms` with `identity`, `p_s` above the grid's top (0.5 Pa) and beyond its bottom (2e6
  Pa), `decay_above` with the stop below `p_s`, `increase_below` with the stop above it, and
  `linear_ln_p` with a stop of 0. Each message names the key and the rule, and no file is written.
  Nothing else is refused (ruling 10).

## 3. The registered transfer product, and the driver fix it required

**§1a.** After the rulings were committed (`a315171`, `7ef249b`) the tree was clean, and
`casspian-forward forward/lindal_transfer/lindal_transfer.toml` rebuilt the product in 81 s, stamped
`7ef249b`. Compared with the `-dirty` copy on disk (the last `step04_5` run, at `4fc1cc7-dirty`),
all 88 computed variables of the root, `reference_surface`, `isobars`, `estimate` and
`anchors/lindal/transfer` groups are array-equal, and no attribute of the root, `estimate` or
`transfer_record` differs other than the writer's stamps. `r0(10 N)` is 60,128,612.966 m (the
variable is `reference_surface/reference_surface_radius_m` at the target node, for Step 3 check 3)
and the pressure identity 5.777e-07, both the accepted Step 5 values. The product and its
`.gitignore` exception are committed at `7276f68`.

**The driver.** `step04_5` run 1 writes `forward/lindal_transfer/output/lindal_transfer_profile.nc`,
and the driver set aside only the transfer run's inputs, so once the product was committed every
regression would have left it modified and `-dirty`. The driver now copies the product aside with
the inputs, restores it after every suite, and includes it in the closing digest, whose line now
reads `transfer inputs and product restored`. Proven by `bash tests/run_regression.sh step04_5`:
11 of 11 at its reference count in 3689 s, the transfer inputs and product restored equal, and the
product never listed among the dirty paths after the suite.

## 4. The regression

Step 1 adds a module no existing code imports, so under SPEC_04 §0 it reruns its own suite only
(SPEC_05 Step 1). The driver change reaches the restore of one file that one suite writes, so it
reruns that suite, `step04_5`. Both were run through the driver:

| Suite | Checks | Reference | Time |
|---|---|---|---|
| `step05_1/accept_step05_1` | 11 of 11 | 11 | 20 s |
| `step04_5/accept_step04_5` | 11 of 11 | 11 | 3689 s |

The other 25 suites are not reached: no code they import changed, and the one committed file this
step adds (the transfer product) is one only `step04_5` writes.

## 5. Findings

1. **`pip install -e . --no-deps` after pulling.** The console entry exists only after a reinstall.
   On the Linux clone that is needed before `casspian-wind-shear` runs; it is Step 2 deliverable 6's
   runbook line.
2. **The regression driver restored every registered file but one.** Found while registering the
   transfer product (section 3) and fixed in this step at the author's direction.

## 6. What is committed on acceptance

`src/casspian/tools/wind/shear.py`, `pyproject.toml`, `tests/step05_1/accept_step05_1.py`,
`tests/run_regression.sh` (the transfer product's restore and the `step05_1` row), and this report.
