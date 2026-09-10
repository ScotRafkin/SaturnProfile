# REPORT 01, Step 1. `lib.constants` and the netCDF conventions layer

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.3, Step 1, against
`SPEC_00_Architecture_and_Data_Files.md` v0.3. Status: reported, awaiting review.

**The working tree is not committed.** Per the instruction of 10 September, a step is committed
only after its report is accepted. Item 0 below is the exception, because the specification
makes it its own commit before any code.

---

## 1. What was built

**Item 0, line endings. Committed and pushed as `27d73f5`.** `.gitattributes` at the repository
root pins `*.csv`, `*.toml`, `*.md`, `*.py`, `*.txt`, `.gitignore` and `.gitattributes` to
`text eol=lf`, and marks `*.nc`, `*.png`, `*.svg` and `*.pdf` as `binary`. `git add
--renormalize .` staged nothing beyond `.gitattributes` itself, which shows the index blobs were
already LF; the exposure was on checkout, since `core.autocrlf` is `true` on this machine and
would have written CRLF on a fresh clone. Every tracked text file was then checked byte by byte
for CR: none carries one. `git check-attr` confirms the rules resolve as intended, including
`binary` on the figure PNG.

A second commit, `251bca3`, carries the author's own uncommitted documents (SPEC_00 v0.3,
SPEC_01 v0.3, `REVIEW_01_step0.md`, `STATE.md`, `notes.md`). It is separate so that item 0
stayed the single clean commit the specification asks for.

**Deliverable 1a, `src/casspian/lib/constants.py`.** The eight CODATA 2018 values named in the
specification, each with a comment giving its CODATA name and its units, plus
`CODATA_RELEASE = "2018"`. Nothing else, and no Saturn quantity. The molar gas constant is
spelled `MOLAR_GAS_CONSTANT` because SPEC_00 section 3.1 forbids anything named `R`.

Self consistency, checked and reported rather than assumed:
`STANDARD_PRESSURE / (BOLTZMANN_CONSTANT * STANDARD_TEMPERATURE)` gives 2.686780111798444e25
against the tabulated Loschmidt constant 2.686780111e25, a relative difference of 3.0e-10, which
is the published truncation. `BOLTZMANN_CONSTANT * AVOGADRO_CONSTANT` gives 8.31446261815324
against the tabulated molar gas constant 8.314462618, a relative difference of 1.8e-11, likewise
the truncation. Both tabulated values are exact in CODATA 2018; the differences are in the digits
the specification asked to be written down, not in the constants.

**Deliverable 1b, `src/casspian/lib/schema.py`.** The registry and the validation layer. It
holds the fixed strings (`CONVENTIONS`, `HARMONIC_CONVENTION`), the vocabularies (provenance,
uncertainty kind, absent meaning, role, vertical coordinate, positive, direction), the forbidden
bare names, the split of the section 5 globals into those `write` fills and those the caller must
supply, and a `KindSpec` for each of the eight kinds giving its schema version, dimensions,
variables, required globals, required groups, and whether per variable attributes are enforced.
All eight schema versions are 1.

`validate` runs, in order: the naming rule; the globals; the dimensions, including the
`<dimension>_absent_meaning` rule and the demand that a `point` file state its single coordinate;
the required variables and their uncertainty companions; the per variable attributes; the
coordinate rules (strict monotonicity, `positive`, and `direction` on pressure); then the two
kind specific checks, mole fractions summing to one for kind C and the component sum identity
for kind W. Every refusal names the item at fault and cites the SPEC_00 section it comes from.

**Deliverable 1b, `src/casspian/lib/io.py`.** `sha256`, `input_hash_entry`, `input_hashes`,
`history_append`, `package_version`, `git_commit`, `group_paths`, `write` and `read`.

`write(path, dataset, kind, groups=None, created_by="casspian")` fills the seven writer globals,
appends a `history` line, checks attribute types and variable types, validates, and writes
netCDF-4. `groups` is a mapping of group path to Dataset, which is how kinds C, N and `raw` carry
their structure; nested paths such as `scalars/latitude` work. `read(path, kind)` checks the
declared kind, checks the schema version is not newer than the reader, checks the writer globals
are present, validates, and checks that every required group exists.

The v0.3 version rule is implemented in `package_version`: the version comes from
`importlib.metadata`, and if it disagrees with `casspian.__version__` the call refuses. `write`
obtains its version only through that function, so no file can be written while the two disagree.

## 2. Deviations, decisions, and defects found

**Defects in the specification, raised not patched.**

1. **The `provenance` vocabulary conflicts with kind T.** SPEC_00 section 5 lists seven allowed
   values: `measured`, `derived`, `assumed`, `inferred`, `interpolated`, `parameterized`,
   `extrapolated`. The kind T table in section 6.1 gives the provenance of `pressure_Pa` and
   `temperature_K` as "`derived` or `modeled`". `modeled` is not in the section 5 list. The
   registry implements the section 5 list, because section 8 calls section 5 the one validation
   layer, so a kind T file whose provenance is `modeled` is currently refused. This must be
   settled before Step 9 writes the first kind T file. Either `modeled` is added to section 5 or
   section 6.1 should say `derived`.

2. **The `uncertainty_kind` vocabulary is incomplete.** Section 5 allows `1sigma`, `range` and
   `stated`. SPEC_01 Step 7 requires the Lindal wind file to carry
   `uncertainty_kind = "bin_std"`. I added `bin_std` to the allowed set so that Step 7 can be
   built as specified, which is a change to a section 5 vocabulary made in code. Confirm it, or
   name the value differently.

3. **`HARMONIC_CONVENTION` had to be authored.** Step 3 says the string is "fixed as a constant
   in `lib.schema`", and section 6.4 says what it must state (the potential form of handoff
   section 9A.5, that the Legendre argument is the sine of planetocentric latitude, and that
   `g_N` is positive inward), but the string itself is written nowhere. I wrote one to that
   description. **This needs the author's word for word approval before Step 3 writes the first
   kind G file**, because the reader refuses on a single altered character, so changing the
   string later invalidates every kind G file already written. The current text is in
   `schema.py` and is quoted here in full:

   > V = -(GM/r) [1 - sum_l J_l (R_norm/r)^l P_l(sin(phi_c))], even degrees l; the Legendre
   > argument is the sine of planetocentric latitude; J values are unscaled; g_N = -dV/dr is
   > positive inward. Handoff section 9A.5.

**Decisions taken, and why.**

4. **`read` returns a `Dataset` for a file with no groups and an `xarray.DataTree` for a file
   with groups.** The specification says it returns a Dataset with groups accessible. A Dataset
   has no group concept, so a single return type cannot do both; a DataTree always would make
   the common case awkward. Flagged because it is an interface the later steps depend on.

5. **The kind T source profile globals are declared but not enforced.** Section 6.1 lists eleven
   globals that apply "for a source profile" (spacecraft, event, observation date, frequency
   bands, latitude definition, height datum, source top boundary, and the three source
   citations, plus `raw_bundle`), and no attribute distinguishes a source profile from a
   retrieval or from the model's delivered field. They are recorded in the registry as a named
   conditional set and are not checked, since guessing the discriminator would refuse valid
   files. A discriminator is needed before Step 9. The obvious candidate is that a source
   profile is the instance that is `point` in latitude, but the model's delivered field for one
   target latitude may also be a point, so I did not assume it.

6. **Kind N requires only its six `inputs/*` groups.** Section 6.7 also names `manifest` and
   `reduction_record`, but describes both as carrying attributes rather than variables, and does
   not say whether they are groups or attribute prefixes. Not required until that is settled.

7. **Tolerances that the specification states in words.** The kind W component identity is
   required to hold "to round-off"; I implemented 1e-12 relative to the largest absolute value
   in `u_total`, so the check scales with the field rather than with an absolute wind speed. The
   kind C mole fraction check uses the 1e-9 of section 6.2. Note that SPEC_01 Step 8 asks the
   composition tool to produce fractions summing to one to 1e-12, which is a tighter demand on
   the tool than the reader's refusal threshold; both stand as written.

8. **Type rules are enforced on write.** Section 5 says float64 for physical quantities and int8
   for status flags with `flag_values` and `flag_meanings`. `write` refuses a float32 variable
   and refuses a flag variable that is not int8 or that lacks either attribute. This is stricter
   than the letter of section 8, which does not list type faults among the refusals, and it
   catches the mistake at the moment it is made rather than at the reader.

9. **`raw` is exempt from the per variable attribute and coordinate rules**, per section 2.2.1,
   which gives it no standard schema beyond the section 5 globals. The naming rule of section 5
   still applies to it, since Step 1 acceptance requires a bare `latitude` to be refused.

10. **`codata_release` and `input_hashes` are not demanded of every file.** Section 5 qualifies
    them ("for any file that used a constant of nature", "for any derived file"). In practice
    `write` always fills `codata_release`, so it is always present.

11. **`git_commit` reports a dirty tree as `<sha>-dirty`.** This matters more under the new
    commit on acceptance rule than it would have before: every file a tool writes during a step
    is written from an uncommitted tree, so `casspian_git_commit` will carry the `-dirty` suffix
    on every product built before its step is accepted. The alternative, recording a bare
    commit that does not describe the code that ran, would make the provenance attribute say
    something false. See question 2 below.

## 3. Acceptance results

The five checks of Step 1, plus one beyond it, run by
`accept_step1.py`. Every check prints its actual result.

| Check | Result |
|---|---|
| Write a two level `raw` file with one variable and read it back with every section 5 global present and correctly typed | **Pass.** `level` size 2, one variable `example_value_Pa` reading back `[1.0, 2.0]`. All twelve globals present: `Conventions = 'CF-1.10, CASSPIAN-0.1'`, `casspian_kind = 'raw'`, `casspian_schema_version = 1` (integer), `title`, `profile_or_run`, `role = 'reduction'`, `source`, `created_by = 'accept_step1 0.1.0'`, `created_at = '2026-09-10T18:12:37Z'`, `casspian_git_commit = '27d73f5b84ce6c553d65210abece6af8e2ec46c8-dirty'`, `codata_release = '2018'`, `history = '2026-09-10T18:12:37Z written by accept_step1 0.1.0 as kind raw'`. |
| `read` refuses it when asked for kind `thermo` | **Pass.** Message: `acceptance_raw.nc: declared kind is 'raw' but 'thermo' was asked for (SPEC_00 section 8).` |
| `read` refuses a file whose `casspian_kind` was edited to `wind` but which lacks the wind variables | **Pass.** The attribute was edited in place with netCDF4 on a copy. Message: `acceptance_raw_edited_to_wind.nc: required global attribute 'rotation_system_name' is missing (kind wind, SPEC_00 section 5 and 6).` The refusal comes on the first missing item, which is a global rather than a variable; the variables are equally absent and would be named next. |
| `sha256` matches `sha256sum` on the command line | **Pass.** Both give `536b8fd8a705f17e463bc3584ff27c8b94c74ccb45b2e16e69f01382eb15ee31`. One artifact worth recording: GNU coreutils prefixes its output line with a backslash whenever the filename contains one, which every Windows path does, so the comparison strips a leading backslash. The digest itself is untouched. |
| A dataset with a variable named `latitude` is refused by `write` | **Pass.** Message: `must_not_exist.nc: the name 'latitude' carries no convention and no unit. SPEC_00 section 5 requires latitude_planetocentric_deg or latitude_planetographic_deg, radius_m or height_m, pressure_Pa.` No file was created. |
| Beyond the specification: groups survive a write and read round trip | **Pass.** A `raw` file written with `table1`, `scalars`, and the nested `scalars/latitude` reads back with all three groups present and `scalars/latitude` carrying `planetographic_deg = 36.3`. Run because Step 2 writes a group bearing raw bundle with exactly this shape. |

Also checked, outside the listed acceptance: the version rule refuses as intended. With
`casspian.__version__` forced to `9.9.9` against an installed `0.1.0`, `package_version` raises
`version disagreement: importlib.metadata reports '0.1.0' but casspian.__version__ is '9.9.9'`.

No em dash or en dash appears in any file written in this step.

## 4. Questions for review

1. **Where do acceptance scripts live?** SPEC_00 Figure 1 has no directory for them, and Step 0
   deleted `tests/`. The Step 1 script therefore ran from a scratch directory outside the
   repository and is not preserved. That is tolerable while the checks are structural, but
   Steps 4, 5 and 6 have physical acceptance numbers (8.951 m/s², 60,244 km, 30.8185 degrees)
   that will want to be reproducible by a second reader, and Step 4 asks for a deliberately
   wrong coefficient to be run as a discriminating check. Proposing a `checks/` directory at
   the repository root, outside `src/`, holding one script per step. Not created, since it
   would be a change to the layout of SPEC_00 section 2.

2. **Does a product written from a dirty tree need to be rebuilt after acceptance?** Under the
   commit on acceptance rule, every netCDF a tool writes during a step records
   `casspian_git_commit` as `<sha>-dirty`, where `<sha>` is the previous accepted step. From
   Step 2 the products are real data files, not throwaways. Either the tools are rerun after
   the acceptance commit so the products carry a clean commit, or `-dirty` is accepted in the
   record for products built during development. Worth settling before Step 2 writes the raw
   bundle.

3. **The three specification defects in section 2 above**, in order of urgency:
   `HARMONIC_CONVENTION` wording before Step 3; the `provenance` conflict before Step 9; the
   `uncertainty_kind` addition before Step 7.

## 5. Next step

Step 2, stage one of the Lindal tool and the raw bundle, is not started and is not to be begun
until this report has been reviewed.

---

## 6. Changes applied after review

Added 10 September 2026, after `REVIEW_01_step1.md`, so that this report and the committed code
agree. The four changes the review asked for were made, and verifying them exposed two further
defects in SPEC_00 that are raised here rather than patched around silently.

**The four changes.**

1. `HARMONIC_CONVENTION` is now the code `"CASSPIAN-J1"`, with the defining sentence in the
   comment beside it, and `KNOWN_HARMONIC_CONVENTIONS` holds the codes the reader accepts. A
   kind G file carrying any other code is refused by code, not by string comparison.
2. `UNCERTAINTY_KINDS` is back to `{"1sigma", "range", "stated"}`. `uncertainty_method` is
   carried as free text on any uncertainty variable and is never validated.
3. `PROVENANCE_VALUES` gains `modeled`, with a comment saying it belongs only to files
   `casspian.forward` writes.
4. Kind T requires `thermo_instance` in `{source_profile, retrieval, model_output}`. The
   thirteen source profile globals are collected in `SOURCE_PROFILE_GLOBALS`, required when the
   instance is `source_profile`, and forbidden otherwise. A source profile must also state its
   latitude as a global in the source's convention.

Also applied from the review's other decisions: `read` now returns a `DataTree` for C, N and
`raw` and a `Dataset` for T, D, G, R and W, fixed by the kind through a `grouped` field on
`KindSpec` rather than by whether a particular file happens to have groups; kind N requires the
`manifest` and `reduction_record` groups alongside the six `inputs/*`; and `reports/` is in
`.gitignore`, with the Step 1 scripts and their printed output under `reports/step1/`.

**Two further defects in SPEC_00, raised not patched.**

5. **The `positive` attribute cannot be required of every coordinate.** Section 5 says every
   coordinate variable "carries a `positive` attribute (`up` or `down`)". In CF that attribute
   is defined for vertical coordinates. Kind G's `degree`, kind C's `species` and kind D's
   `surface` are coordinates with no vertical sense, and demanding `positive` of them would
   make the Step 3 acceptance impossible to satisfy without writing a meaningless attribute.
   The rule is now applied only to coordinate names containing `pressure`, `geopotential`,
   `height` or `radius`. Strict monotonicity is still checked on every coordinate. **This
   narrows a section 5 rule and needs a ruling.**

6. **Section 5 and the section 6 tables disagree on how an uncertainty companion is named.**
   Section 5 says `<variable>_uncertainty`. Section 6.4 names it `GM_uncertainty_m3s2` and
   section 6.6 names it `u_total_uncertainty_ms`, both keeping the unit suffix last, which is
   the better reading since `GM_m3s2_uncertainty` puts the unit in the middle of the name.
   `VarSpec` now carries an optional `uncertainty_name` so a kind table can give the companion
   its own name, and the two above are set that way; everything else keeps the section 5
   default. The attribute check matches on the `_uncertainty` token rather than on the ending,
   so both spellings are recognized. **This needs a ruling too**, since it decides the variable
   names in kinds G, W and N before Steps 3, 7 and 9 write them.

**Verification.** `reports/step1/verify_review_changes.py` exercises all six items and prints
its results to `reports/step1/verify_output.txt`. Thirteen checks, all passing: the known code
validates and `CASSPIAN-J2` is refused; `UNCERTAINTY_KINDS` holds exactly the three section 5
values and `bin_std` is now refused while `1sigma` with a free text `uncertainty_method` is
accepted; `modeled` is in the vocabulary; a kind T file without `thermo_instance` is refused, a
complete `source_profile` validates, one missing `raw_bundle` is refused, a `model_output`
carrying source profile globals is refused, and a `model_output` without them validates; `read`
returns a `Dataset` for kind R and a `DataTree` for kind `raw` when neither file has any group,
which is the point of fixing the type by kind; and a harmonic `degree` coordinate validates
without a `positive` attribute. The six Step 1 acceptance checks were rerun after every change
and still pass.
