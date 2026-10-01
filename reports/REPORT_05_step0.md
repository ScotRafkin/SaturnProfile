# REPORT 05, Step 0. The scripts under tests/, the Linux clone, and the regression of record

Coding agent, 1 October 2026. Specification: `docs/specs/SPEC_05_Shear_Experiments.md` v0.4,
section 0, and `Claude outputs/INSTRUCTIONS_git_2026-09-29.md` (items 1 to 7). The Linux runs are
the author's, on the CentOS 7 clone of `docs/RUNBOOK.md`; their logs are in
`reports/regression/linux_2026-10-01/` on the Windows machine.

**Result.** The acceptance scripts are under `tests/` (`84638d6`), the regression driver asserts
reference counts, and the registered products are committed (`4d32efe`, decision D1). The Linux
regression of record ran 21 of 25 suites at their reference counts. The four that did not are two
causes, neither a wrong calculation: line endings in three build control files, a real provenance
break, now fixed by a rebuild whose 518 computed variables are identical to the products they
replace (section 4); and the closure rerun's p and T, which differ from the Windows-registered
product by 1 and 2 units in the last place on Linux (section 5), ruled on by
REVIEW_05_step0 (a bound of 1e-14 relative). At the author's direction Step 0 closes on Windows
and the clone is not rerun (section 10).

## 1. What Step 0 did, by commit

| Commit | What |
|---|---|
| `9cab59d` | The runbook v0.2 and the 28 September handoff (item 1). |
| `84638d6` | 48 files under `tests/`: every acceptance script, the helpers they import, the cited diagnose scripts, the most complete version of the sweep, rebuild and candidate tools, the `step03_3` fixtures, and `tests/run_regression.sh` with its table of reference counts; `.gitignore`, `.gitattributes` (`*.sh` LF), `pyproject.toml` (`matplotlib`); SPEC_00 v0.22 (items 2 to 5). |
| `348bcff` | The runbook v0.3; the author's two working documents ignored. |
| `4d32efe` | D1: the 17 registered products and the 13-file `step04_0` swept fixture committed; the `step03_3` anchor path; `lib.io.git_commit` without `git -C`; the driver lists dirty paths. |
| `7276f68` | The registered transfer product (SPEC_05 §1a). |
| `7a6107d`, `38cd50d` | The sweep tool for committed products (section 4). |

The note to the author of 29 September lists every file moved and every path changed in the move;
it is not repeated here.

## 2. The first Linux run, and what it found

The first full run on the clone (26 suites, `step04_4` included) gave 22 of 26. Its four
mismatches and the fixes, all in `4d32efe`:

- **`step02_1` and `step02_6` check 1** compared SHA-256 with hashes recorded from Windows
  builds, and the clone had rebuilt the products. Fixed by D1: the registered products are
  committed, so every platform reads the same bytes.
- **`step03_3` crashed**: an anchor path hard-coded relative to the old script location
  (`../../../candidate/`), which the Windows proof missed because the stale outputs under
  `reports/step*/` were still present. Fixed by computing the path to the committed fixture; the
  Windows proof was rerun with every `reports/step*/` directory moved out of the tree.
- **`step04_0` crashed**: a fixture (`SWEPT`) defaulted to an absolute path in an earlier session's
  scratch directory on Windows. Fixed by committing its 13 files under
  `tests/step04_0/fixtures/swept/`.
- **Every clone product recorded `casspian_git_commit = "unknown (not a git checkout)"`**: CentOS
  7's git is 1.8.3.1, and `git -C` arrived in 1.8.5. Fixed by C1, `cwd=` in place of `-C`; the
  run of record shows the clone's products now carry the real commit.

`step03_3` check 9 also failed in that run, comparing the clone's rebuilt wind with the
Windows-built candidate fixture: 2 of 361 latitudes, at most 7.6e-16 relative in `u_total_ms` and
6.6e-15 in its uncertainty, every other variable identical, the planetocentric to planetographic
conversion rounding differently in the last bit. Under D1 the clone no longer rebuilds what that
check reads, and the run of record passes it.

## 3. The Linux regression of record

At `4d32efe`, every suite but `step04_4`, which passed 11 of 11 in 27,247 s against the Linux
built products in the first run (the author's ruling to carry it). One dirty path throughout,
`regression_console.txt`, the run's own log.

| Suite | Checks | Time (s) | | Suite | Checks | Time (s) |
|---|---|---|---|---|---|---|
| `step1/accept_step1` | 6 of 6 | 1 | | `step02_4` | 9 of 9 | 476 |
| `step1/verify_review_changes` | 14 of 14 | 1 | | `step02_5` | 7 of 7 | 416 |
| `step2` | 8 of 8 | 1 | | `step02_6` | 7 of 7 | 988 |
| `step3` | 5 of 5 | 0 | | `step03_1` | 9 of 9 | 1 |
| `step4` | 7 of 7 | 1 | | `step03_2` | 8 of 8 | 0 |
| `step5` | 7 of 7 | 121 | | `step03_3` | **15 of 16** | 54 |
| `step6` | 6 of 6 | 84 | | `step03_4` | 13 of 13 | 8 |
| `step7` | 8 of 8 | 225 | | `step04_0` | **14 of 15** | 60 |
| `step8` | 9 of 9 | 1 | | `step04_1` | **12 of 13** | 30 |
| `step9` | 6 of 6 | 0 | | `step04_2` | 13 of 13 | 3767 |
| `step02_1` | 6 of 6 | 0 | | `step04_3` | 9 of 9 | 8189 |
| `step02_2` | 7 of 7 | 140 | | `step04_5` | **10 of 11** | 6082 |
| `step02_3` | 9 of 9 | 299 | | `step04_4` (first run) | 11 of 11 | 27247 |

The registered kind N, the closure inputs and product, and the transfer inputs were restored
exactly. The clone runs the long suites 1.7 to 1.9 times slower than Windows.

## 4. Cause 1: line endings in three build control files

**Diagnosis.** `git ls-files --eol` listed 13 tracked files whose Windows working copies were CRLF
while git holds LF; they predate `.gitattributes` pinning text to LF, and git never rewrote them.
Three are build control files the products hash: `occul_data/lindal/lindal_build.toml`,
`forward/lindal_closure/lindal_closure_build.toml` and
`forward/lindal_transfer/lindal_transfer_build.toml`. The other ten (eight reports and
specifications, `STATE.md`, `src/casspian/lib/control.py`) are hashed by nothing. Of the 18
registered files, 17 recorded the CRLF bytes of one of the three; only `raw/lindal_raw.nc`, which
hashes no build file, was clean. A clean checkout on either platform gives LF, so the recorded
hashes disagreed: warnings only, no value moved. `step03_3` check 7 asserts that an unedited kind T
reads silently, so it failed, rightly, and `step04_0` check 14 reruns `step03_3`.

**Fix.**

1. The 13 working copies rewritten to LF, each verified byte-equal to `git show HEAD:<file>`;
   `git ls-files --eol` reports no CRLF, and nothing was staged because git already held these
   bytes. No untracked source or control file carries CRLF.
2. The sweep tool made to work on committed products (`7a6107d`, `38cd50d`). Rewriting the first
   committed product dirtied the tree, and `lib.io.git_commit` would have stamped every later file
   `-dirty`, which `refrac` and `forward` refuse. At the author's choice (option A) the tool marks
   the 18 registered products `skip-worktree` for the sweep only and clears the mark in a
   `finally`, so the stamp checks every code and control file and ignores only the outputs being
   rebuilt. It also rebuilds the transfer product, registered since the tool was written, through
   `forward.transfer.run`; the first attempt called the closure driver, stopped after 17 files,
   and was discarded by returning the tree to its committed bytes.
3. `python tests/step04_0/sweep.py` on the clean tree at `38cd50d`, 292 s: all 18 files carry
   `38cd50df9620...` and none `-dirty`; the marks were cleared.

**Evidence that nothing moved.** Against copies set aside before the sweep, every variable of
every group of all 18 files, the embedded input copies included, 518 in all, is array-equal. The
attributes that differ are the writer's stamps and those that carry a hash or a resolved path
(`input_hashes`, `input_sha256_*`, `latitude_conversion_inputs`, `control_file`, `raw_bundle`,
`resolved`). All 117 recorded hashes in the 18 files now match the files on disk, and the build
control files among them match git's LF bytes. The sweep's four warnings are
`casspian-lindal-inputs` reading the previous thermo and geodesy files, which carried the old
hashes, before rewriting them.

**What changes with it.** Every registered file has new bytes, `lindal_refractivity.nc` included
(`920304440ee4...` becomes `b79c60cb1cba...`). No suite pins a hash in code; `step02_1` check 1 and
`step02_6` check 1 read the registered hashes from the reports' hash tables, so section 7's table
is the record and both suites now read this report as well (section 6).

## 5. Cause 2: the closure rerun's p and T on Linux

`step04_1` check 10 and `step04_5` check 11 are one comparison: the closure namelist rerun and its
computed columns compared with the registered closure product by `array_equal`. On the clone N,
geopotential, refractivity and molar mass came back bit-identical; p and T did not. Measured by
the author on the clone from the rerun `step04_5` left on disk:

| Variable | Levels differing | Largest relative | Largest in units in the last place |
|---|---|---|---|
| `pressure_Pa` | 3 of 66 | 1.789e-16 | 1 |
| `temperature_K` | 2 of 66 | 2.705e-16 | 2 |

p is integrated from N through `exp` and `log`, which round differently in the last bit in the two
platforms' math libraries, and T follows from p. Both checks now print, beside the bit-identity,
the largest relative difference and the largest in units in the last place of every compared
variable, so the size is in every future run (section 6).

**For the review.** The proposed ruling (the reviewing agent, 1 October) is to keep the bit-identity
on the platform that registered the products and to bound the difference elsewhere, for example at
1e-14, 37 times the largest measured. REVIEW_05_step0 ruled: one bound on every platform, 1e-14 relative,
with the bit-identity and the last-place difference still printed (section 8).

## 6. Suites changed in this step

| Suite | Change | Why |
|---|---|---|
| `step03_3` | the run cases' anchor path computed to the committed fixture | section 2 |
| `step04_0` | `SWEPT` defaults to `tests/step04_0/fixtures/swept/` | section 2 |
| `step02_1`, `step02_6` | `reports/REPORT_05_step0.md` added to the reports whose hash tables check 1 reads | section 4 |
| `step04_1` check 10, `step04_5` check 11 | the largest relative and last-place difference printed beside each bit-identity; pass condition unchanged | section 5 |
| `tests/step04_0/sweep.py` | `skip-worktree` on the registered products for the sweep; the transfer product rebuilt | section 4 |
| `tests/run_regression.sh` | dirty paths listed; the transfer product set aside and restored (Step 1) | sections 1 and 2 |

No reference count changes.

## 7. Hashes, the sweep at `38cd50d`

Every file below was rebuilt by `tests/step04_0/sweep.py` on the clean tree at `38cd50d` and
carries `38cd50df9620190224fa018b1ed99e94dce5126c`, none `-dirty`. The reduction chain's rows are
named by file, the form the `step02_1` and `step02_6` suites read, and supersede REPORT_04_step0
section 6. The manifest's hash is unchanged.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_reduction.toml` | manifest | `97e542d6665f525002d2a144a52da5633892f007beb28f778434af245a88dfee` |
| `lindal_raw.nc` | raw | `04d0c34e829223af5d83cc7ccefb1307deb4555c003c7c5d3862f79925667468` |
| `lindal_thermo.nc` | thermo | `8a5a8eca517eeface81ab14934a4c22b911bbad1beee08112753a6d85e2ae227` |
| `lindal_geodesy.nc` | geodesy | `54da2d1c9669d662daa254d632018599ccc2021399d9ab7777b6366730f94cb5` |
| `lindal_gravity.nc` | gravity | `3d0a8f247ab454b08ce3d2e2c779ebcca8f51d121bea2406bb0b9a5348c2a8cc` |
| `lindal_rotation.nc` | rotation | `c007a790747939b267a72c788d2778684691760c4e318e4b813f649c73df8eae` |
| `lindal_wind.nc` | wind | `1c2397d49dba7ff32aa9b8369ee35e972edb5bad2492d1be144223f79ef54d91` |
| `lindal_composition.nc` | composition | `370c350b34e6c8d5ba65d058dfe47b8963dc8373cbbd20575255e5eec72049d3` |
| `lindal_refractivity.nc` | refractivity | `b79c60cb1cba26181b86a59184846dd388380f7d7b4671cb8535cd7042d84120` |

The forward runs' registered files, by path:

| File | Kind | SHA-256 |
|---|---|---|
| `forward/lindal_closure/inputs/lindal_closure_gravity.nc` | gravity | `9660980aac2ceb87ffe7e9a86a1d477d15f2e5e7b633d0702c3ce6829cbff8a3` |
| `forward/lindal_closure/inputs/lindal_closure_rotation.nc` | rotation | `93977739e61a2aced301fcb6f2cbbf143912ac6706937d0e056edea660108f75` |
| `forward/lindal_closure/inputs/lindal_closure_wind.nc` | wind | `d4733d45bdbbb6baae3c8ac2c2aeb6dfd1e4c327f12ea384066fa77b92c24327` |
| `forward/lindal_closure/inputs/lindal_closure_composition.nc` | composition | `b2e56da52084f7b23d4b6b91b3c90a41d95a0d8e3522a738779b0478f7d49d2e` |
| `forward/lindal_closure/output/lindal_closure_profile.nc` | profile | `f902919795f368fa6aa4c1e5ed4067f8ced2a469a8d165164fcbb3d49d14c78b` |
| `forward/lindal_transfer/inputs/lindal_transfer_gravity.nc` | gravity | `225512c0ed0ce774a6deb31cd3965ba88d8b27bf0e48081006412966656401cf` |
| `forward/lindal_transfer/inputs/lindal_transfer_rotation.nc` | rotation | `1932ccdb733fea224686f1f341c3db06dba14f0f44d00d38525460b2a1cddb2b` |
| `forward/lindal_transfer/inputs/lindal_transfer_wind.nc` | wind | `d08e3531ab8e91e7727cd7755cadf32d5e3de54a4916af9e5388b51fc77f3e54` |
| `forward/lindal_transfer/inputs/lindal_transfer_composition.nc` | composition | `b5ff248a7fa4546cc84f90fb3db4181043a02d481ef2e2f9f9200a9b9a0add44` |
| `forward/lindal_transfer/output/lindal_transfer_profile.nc` | profile | `dd5c940364683310410800c8645f16788f9a0d2e9b4aae9f7e1d2f98cc93b29d` |

## 8. The regression for this change

**On Windows**, the suites the rebuild and the suite changes reach, through the driver on the
rebuilt products: `step02_1` and `step02_6` (they read the new hash table), `step03_3` and
`step04_0` (the provenance of check 7 and its rerun), `step04_1` and `step04_5` (the departure
output, and the bit-identity against the rebuilt closure product). The products' values are
unchanged (section 4), so no other suite is reached.

| Suite | Checks | Reference | Time (s) |
|---|---|---|---|
| `step02_1/accept_step02_1` | 6 of 6 | 6 | 11 |
| `step02_6/accept_step02_6` | 7 of 7 | 7 | 666 |
| `step03_3/accept_step03_3` | 16 of 16 | 16 | 173 |
| `step04_0/accept_step04_0` | 15 of 15 | 15 | 296 |
| `step04_1/accept_step04_1` | 13 of 13 | 13 | 33 |
| `step04_5/accept_step04_5` | 11 of 11 | 11 | 3885 |

6 of 6 suites at their reference counts; the registered kind N (`b79c60cb...`), the closure inputs
and product, and the transfer inputs and product restored exactly. `step03_3` check 7 now reads the
unedited kind T silently. The new departure lines print `0.000e+00 relative (0 ulp)` for every
compared variable of `step04_1` check 10 and `step04_5` check 11, as they must on the platform that
registered the products.

**The ruling of REVIEW_05_step0, on Windows.** `step04_1` check 10 and `step04_5` check 11 pass
when every compared variable agrees with the registered product to 1e-14 relative, and exactly
where the registered value is zero; each still prints its bit-identity and its last-place
difference. Rerun on Windows through the driver: `step04_1` 13 of 13 (36 s) and `step04_5` 11 of 11
(3735 s), both at their reference counts, every compared variable within the bound at
`0.000e+00 relative (0 ulp)`, and every registered file restored exactly.

**On the clone: not rerun, at the author's direction** (section 10). Its outcome is known from the
evidence in hand: the causes of the four mismatches are fixed and verified on Windows, and the two
re-pointed checks will print the 1 and 2 last-place departures of section 5, under the 1e-14 bound.

## 9. Findings

1. **Line endings were a real provenance break** (section 4). Fixed; the CRLF scan covers every
   tracked text file.
2. **Bit-identity across platforms cannot hold for values computed through `exp` and `log`**
   (section 5). For the review.
3. **Committed products and the dirty stamp.** Under D1 a sweep rewrites tracked files, and the
   stamp would mark them `-dirty`. The sweep tool handles it (`skip-worktree`); Step 2, which
   rebuilds and recommits both registered products, will use the same tool.
4. **Two fixtures and one path were missed in the move** and were caught only on a clean clone
   (section 2). The Windows proof now runs with `reports/step*/` moved out of the tree.
5. **git 1.8.3.1 on the clone** stays; C1 removed the only dependency on a newer git.
6. **Where the regression's time goes, for the decision R specification** (the author's direction,
   1 October: the vectorized columns go into that specification). One production transfer run at
   10 N, profiled with figures off on Windows, took 94 s. `lib.mesh.build_columns` took 78 s of it:
   843 columns (421 latitudes, two outer-loop passes), each integrated on its own, with 273,068 calls
   to `mesh.slopes`, each evaluating `lib.gravity.g_eff_vector` (the harmonic series and its Legendre
   polynomials, 58 s) and `lib.windfield.wind_at` (15 s) at one point. The reference surface took
   8 s and reading the inputs 7 s. The cost is Python's per-call overhead, about 0.27 ms per slope,
   not the arithmetic; at the 5,000 m2/s2 spacing `step04_4` is pinned to, each column takes about
   ten times the steps, and its halved-spacing checks four times that again. Integrating every
   column together, one array evaluation of gravity and wind per step for all latitudes, removes
   most of it. It changes `lib`, and has to reproduce every product to the bit or go through a
   ruling.

## 10. The author's direction, 1 October 2026

1. **Step 0 closes on Windows.** REVIEW_05_step0 made closure conditional on a rerun of the six
   suites on the clone. The author waives it: Step 0 closes on the Windows result of section 8 and
   the Linux measurements of sections 3 and 5.
2. **No routine Linux runs.** The clone has shown what it was for: two missed fixtures, a
   hard-coded path, `git -C`, the line-ending provenance break, and a platform difference of 1 and
   2 units in the last place, now bounded. A Linux run is made only when the author chooses (a
   release, a port), and the protocol does not require one.
3. **What a small change reruns.** The rule stays SPEC_04 section 0 at v0.22, rerun what a change
   can reach, made explicit for the long suites: `step04_2` to `step04_5` are rerun only when a
   change touches `lib`, `refrac`, `forward` or the registered products. A change to a tool, a test,
   a figure or a document reruns its own suite.
