# REPORT 01, Step 0. Repository reset and skeleton

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.2, Step 0, against
`SPEC_00_Architecture_and_Data_Files.md` v0.3. Status: reported, awaiting review.

---

## 1. What was built

**Commit 1, `Phase 1 prototype as of September 2026, superseded by SPEC_00`, tagged
`phase1-prototype` at `d7decfd1da5098d65544eeb4b1c65b038d9081eb`.** The working tree was
committed exactly as it stood, including the seven uncommitted Phase 1 source and test files
and the modification to `composition_utils.py`, plus `data/lindal1985_voyager2.nc` added with
`git add -f` because the then current `.gitignore` excluded it under the `*.nc` rule.
`git status --ignored` was checked before tagging: the only remaining ignored entries were
`.pytest_cache/`, `occultation_recon/__pycache__/`, and `tests/__pycache__/`, none of which is
part of the record. The commit and the tag were pushed. The netCDF carries the superseded
36.5 degree latitude attribute and is preserved as history, not as data.

**Commit 2, the reset.** Removed `occultation_recon/`, `tests/`, `docs/this_spec.md`, and
`docs/CASSPIAN_OccultationReconstruction_BuildSpec.md`. Moved `data/lindal1985_table1.csv` to
`occul_data/lindal/raw/lindal_table1.csv` with `git mv`, so the history of that file follows
it, and deleted `data/`. Left `docs/Lindal_et_al_1985_AJ90_1136.pdf`, `docs/specs/`,
`data_static/`, and `occul_data/lindal/raw/lindal_scalars.toml` untouched. No `attic` directory
was created.

Created the skeleton: `src/casspian/{lib,tools,refrac,forward}/__init__.py`,
`src/casspian/tools/{lindal,wind,composition,gravity}/__init__.py`, `forward/README.md` at the
repository root stating what goes there, and `occul_data/lindal/raw/` as above. Each
`__init__.py` carries a one line docstring naming the SPEC_00 section that governs it. No
functional code was written.

Replaced `pyproject.toml` (name `casspian`, `src` layout, `requires-python >= 3.11`,
dependencies `numpy`, `scipy`, `xarray`, `netCDF4`, version `0.1.0`, no console entry points,
with a comment saying they are added as the tools appear), `.gitignore` (caches, build
artifacts, `forward/*/output/`, and `occul_data/*/*.nc` and `occul_data/*/raw/*.nc` with the
comment that these two are removed once a reduction is accepted; `*.nc` is no longer ignored
globally), and `README.md` (ten lines: what the repository is, SPEC_00 as the architecture,
STATE.md as the current state, proprietary to Southwest Research Institute).

Created `docs/specs/STATE.md` with one row per step of SPEC_01, columns step, deliverable,
status, report, date. Step 0 is its first `reported` row.

## 2. Deviations from the specification, and why

1. **`src/casspian/__init__.py` was created, which item 3 does not list.** Item 3 names the four
   subpackage `__init__.py` files and the four tool `__init__.py` files but not the top level
   one. Without it `import casspian`, which the acceptance check requires, cannot succeed. Taken
   as an omission in the enumeration rather than as an instruction.

2. **That file declares `__version__ = "0.1.0"`,** matching the `pyproject.toml` version.
   Nothing in Step 0 needs it. It is written now because SPEC_00 section 5 requires `created_by`
   to carry a tool name and version in every file, and Step 1 is the first writer. If the author
   prefers the version to have a single source of truth read from package metadata, that is a
   one line change in Step 1.

3. **Commit messages carry no AI attribution trailer,** per SPEC_01 section 0. This overrides
   the default commit format of the coding agent. Recorded because the difference is otherwise
   invisible.

4. **`data_static/` and `occul_data/` entered git in commit 1.** Both were untracked when Step 0
   began. Item 1 says to commit the working tree exactly as it is, so they are in the Phase 1
   commit and in the tag, even though they are current data rather than Phase 1 artifacts.
   Item 2 leaves them in place, so they carry forward unchanged. This appears to be the intended
   reading; it is flagged because it puts current data inside a tag named for superseded work.

5. **`.pytest_cache/` was removed from the working tree** along with `tests/`. It was ignored,
   so it was never tracked and its removal appears in neither commit.

## 3. Acceptance results

| Check | Result |
|---|---|
| `git log --oneline` shows two commits above the Phase 1 history | Pass. `1e0bb85` (Phase 1 history), then `d7decfd` (Phase 1 prototype), then the reset commit that carries this report. |
| The tag resolves | Pass. `phase1-prototype` resolves to `d7decfd1da5098d65544eeb4b1c65b038d9081eb`. `git ls-tree -r phase1-prototype` lists `data/lindal1985_voyager2.nc`, so the tag is complete. |
| `pip install -e .` succeeds | Pass. Built `casspian-0.1.0-0.editable-py3-none-any.whl`. All four declared dependencies were already satisfied in the environment. |
| `python -c "import casspian"` runs | Pass. Reports `casspian 0.1.0` from `src/casspian/__init__.py`. All eight subpackages import as well. |
| The tree matches SPEC_00 Figure 1 for every directory that exists at this stage | Pass. The root holds `pyproject.toml`, `README.md`, `.gitignore`, `src/`, `data_static/`, `occul_data/`, `forward/`, `docs/`. `src/casspian/` holds `lib/`, `tools/`, `refrac/`, `forward/`; `tools/` holds `lindal/`, `wind/`, `composition/`, `gravity/`. `occul_data/lindal/raw/` holds `lindal_table1.csv` and `lindal_scalars.toml`. Not yet present, and not due at this stage: `raw/notes.md`, `raw/lindal_raw.nc`, the six standard input files, the manifest, the product, and any `forward/<run>/`. |
| `STATE.md` exists | Pass. `docs/specs/STATE.md`, ten rows, Step 0 the first `reported` row. |

Environment: Python 3.14.3 on Windows 11. The interpreter is newer than the declared floor of
3.11, which is satisfied.

## 4. Questions and observations for review

1. **Line endings.** The repository has no `.gitattributes`, and git reports that LF will be
   replaced by CRLF for every text file on this Windows checkout. Nothing in Step 0 depends on
   it, but the CSV and TOML transcriptions under `data_static/` and `raw/` are about to become
   hashed inputs, and SPEC_00 section 6 makes provenance by the SHA-256 of the file on disk. Two
   checkouts on different platforms would hash the same transcription differently. A
   `.gitattributes` pinning the text transcriptions to LF would remove that risk before Step 2
   records the first hash. Raised, not acted on.

2. **`raw/notes.md`.** SPEC_00 section 2.2 lists `notes.md` under `occul_data/lindal/raw/` as
   transcription notes. Step 0 does not create it and Step 2 names it neither as an input nor as
   a deliverable. Confirming whether it is due, and from whom, before Step 2 is reported.

3. **`data_static/species_master.toml` and the harmonic GM values** are present on disk. SPEC_00
   section 9 lists both as literature tasks in progress that block SPEC_01. Steps 1 and 2 read
   neither, so Step 1 can proceed; Step 3 reads the harmonic sets and Step 8 the species table.
   Flagged so that the open items are closed before those steps rather than at them.

4. Three defects in the specification texts were raised before Step 0 was executed, and were
   fixed by the author in SPEC_01 v0.2 and SPEC_00 v0.3: Step 0 as first written would have
   destroyed `data/lindal1985_voyager2.nc`, since an ignored file does not reach a tag; SPEC_01
   declared a dependency on SPEC_00 v0.2 while v0.3 was the current text; and SPEC_00 cited the
   companion document under a filename it does not have. Step 0 was executed from the revised
   text only.

## 5. Next step

Step 1, `lib.constants` and the netCDF conventions layer, is not started and is not to be begun
until this report has been reviewed.
