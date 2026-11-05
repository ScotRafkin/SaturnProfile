# SPEC 18. Kliore et al. (1980) as a second registered profile

Version 0.4, 5 November 2026. Author of record: S. Rafkin. **Status: closed, Step 1 accepted
(REVIEW_18_step1).**
Depends on SPEC_16.

## 0. Purpose

Reduce the Pioneer 11 exit profile of Kliore et al. (1980, JGR 85, 5857) Table 1 to a kind N, and
show it works as an anchor with a closure run. Then Kliore is ready for the day time propagation
exists, when it becomes a second anchor beside Lindal and a cross-validation of the two. It is
not used against the Voyager-epoch wind (author, 30 October 2026).

**Three layers (author, 30 October 2026).** A tool is specific to its data when the data are
unique, and general when it operates on the same form of data.
1. **Ingest, specific to the source.** Lindal's Table I and Kliore's Table 1 are different data, so
   each has its own preprocessing into the standard kinds. Kliore's undoes the equatorial
   mapping, writes the 100 mbar level (§2) and states his composition.
2. **From the standard kinds on, general.** Kind T is the standard profile. The tools that write the
   composition, gravity, rotation and wind files and the manifest operate on the same form for
   every profile, and are one tool each. Where a Lindal tool mixes the two layers, its general
   part is separated out and shared.
3. **The model, `refrac` onward,** is unchanged. It changes only to fix an error.

**Lindal's answer is unchanged.** Lindal's registered files are never rewritten. Any tool in
Lindal's chain that changes is shown, by rebuilding Lindal's files in a copy, to give them equal
in value to the registered ones. If that holds, nothing is broken and no full regression is run.
A change that cannot meet it is not made; the need is met in Kliore's own code.

Protocol as SPEC_04 §0.

## 1. Settled

1. **Directory and prefix:** `occul_data/kliore1980/`, files prefixed `kliore1980_`. The
   transcription `pioneer11_exit.csv` stays as published, with the 17.5 mbar correction
   (REVIEW_16_step1).
2. **Composition:** Kliore's own, 90 percent hydrogen and 10 percent helium by number, no
   ammonia. The reduction then reproduces his measured refractivity (SPEC_16's self-check). A run
   that uses this anchor chooses its own forward composition, as now.
3. **Latitude:** 9.5 S planetocentric, his table's value.
4. **Radius:** Kliore's equatorial-mapped radius back to the local radius with f = 0.096 and the
   planetocentric reading (SPEC_16).
5. **Gravity:** Null et al. (1981), the Pioneer 11 solution, as Lindal's reduction. **Rotation:**
   System III.
6. **Reduction wind:** Lindal's, the Ingersoll and Pollard (1982) cloud wind, altitude
   independent. A comment says it is a Voyager wind used two years later.
7. **Geodesy (kind D):** Lindal's fitted 100 mbar and 1 bar surfaces. Kliore fitted none, and
   Lindal's fit used Pioneer 11 tracking among its points (Lindal Fig. 9). SPEC_16 found Kliore's
   radii within 7 km of it. The anchor rule is Lindal's, `equatorial_radius`.
8. **Coverage:** 1.8 to 182 mbar. Nothing is extrapolated. `p_b_rule = "anchor_profile_top"`.
9. **No new keys.** Epochs, sources and each choice above go in comments. A key is added only if
   code reads its value.

## 2. The 100 mbar level (ruled, author, 30 October 2026)

The reduction requires its anchor isobar to be a level of kind T and a surface of kind D (SPEC_00
§8, kept by SPEC_R1). Lindal's anchor isobar is 100 mbar. Kliore tabulates at steps in radius,
with 99.6 and 101.9 mbar on either side of it.

Kliore's preprocessing writes one level at exactly 100 mbar into his kind T, its temperature and
radius linear in ln p between the two neighbors (about 89.8 K; the interpolation error is about
0.01 K). The level is named in the build file's comment and the tool's code comment. The
transcription keeps Kliore's 128 levels untouched. `refrac` is the same code for both profiles.

## 3. Step 0. Pre-execution report

Before code, `reports/REPORT_18_preexecution.md` proposes:
1. **The line, tool by tool:** for each tool in Lindal's chain (`build_raw`, `build_inputs`, the
   gravity, rotation, wind and composition builders), whether it is source-specific, general, or
   mixed, and for a mixed one what is separated out. For example, the manifest writing in
   `build_inputs` is general and its reading of Lindal's raw bundle is not; the composition
   builder's ammonia fill from Lindal's raw bundle is Lindal's.
2. **Kliore's preprocessing:** what it reads and writes, what his raw bundle holds (the table, the
   scalars of §1, and the reference geoid's polar radius the wind tool reads), and where it lives.
3. **The Lindal rebuild** in a copy that shows each changed tool gives Lindal's files equal in
   value.

The author and reviewing agent rule on it before Step 1.

## 4. Step 1. The reduction and the closure

**Deliverables.**
1. Kliore's preprocessing, and the general tools as ruled at Step 0.
2. `occul_data/kliore1980/`: the transcription, the raw bundle, the build file, kinds T, C, D, G,
   R and W, the manifest, and the kind N from `refrac`. Committed as Lindal's are.
3. **A closure run,** `forward/kliore1980_closure/`, made from `lindal_closure` with Kliore's
   inputs: the forward production at Kliore's own latitude, as SPEC_03 did for Lindal.

**Acceptance (`tests/step18_1/accept_step18_1.py`, under `reports/step18_1/`).**
1. **Lindal's answer unchanged.** No file under `occul_data/lindal/` or `forward/lindal_*` changes.
   Lindal's files rebuilt in a copy by every changed tool are equal in value to the registered
   ones. With that, no other suite is rerun.
2. **Kliore's refractivity.** Kind N's refractivity against Kliore's tabulated column: RMS and
   largest, in n-units, against the table's 0.05 n-unit rounding. The inserted 100 mbar level is
   excluded, and reported.
3. **Kliore's radius.** Kind N's radius against Kliore's local radius (§1 item 4): the mean
   offset and its spread, against SPEC_16's +9.2 km.
4. **The closure.** The closure product's residuals in temperature and pressure, with the same
   measures and bins SPEC_03 Step 4 used for Lindal, against the envelope of Kliore's printed
   precision (computed in the acceptance script; F6 is not changed). Lindal's bounds do not apply:
   they come from Lindal's table. Pressure identity reported.
5. **Kliore's table, hydrostatically.** His temperatures integrated down in radius with our gravity
   at 9.5 S and his 90/10 composition, against his printed pressures. The drift the dry run found
   below 10 mbar (up to 0.9 percent, mean +5.5e-3 below 100 mbar) is explained before Step 1
   closes: if this integration reproduces it, it is Kliore's gravity against ours, a property of
   the source; if not, the report says what was tried next.

## 5. Ruled on the pre-execution report (author, 5 November 2026)

The tool-by-tool line of its §1 stands: `build_composition` and `build_inputs` are split, the
general parts of `build_inputs` going to `tools/profile_inputs.py`; every other tool is unchanged.
Kliore's preprocessing is `tools/kliore/`, as its §2 proposes. Then:

1. **The label latitude** is computed by stage two as `phi_c + psi(phi_c)` on the geoid `refrac`
   marches, from 9.5 S planetocentric: about 11.89 S planetographic, and `refrac` returns 9.5 S.
2. **Height** is `(r_local - r_local at 100 mbar) / cos psi`, so kind N's radii keep Kliore's spacing.
3. **Uncertainties** are NaN in kind T, with Kliore's ranges kept in the raw bundle. How they become
   an anchor error is ruled when Kliore becomes an anchor.
4. **Kind D's source:** Lindal's two surfaces and reference polar radius transcribed into
   `kliore1980_scalars.toml`, cited, and checked equal to Lindal's kind D. Kliore's build reads
   none of Lindal's files.
5. **No `inserted_level` variable.** The inserted level is identified by its pressure, 100 mbar
   exactly, and documented in comments.
6. **The closure** is measured against Kliore's printed precision (Step 1 acceptance 4), and its
   drift is explained by the hydrostatic check of Step 1 acceptance 5 before Step 1 closes.
7. **Season:** Ls 353.7 for 1 September 1979, by the method that reproduces Lindal's, in comments.

## 6. Revision history

- v0.1, 30 October 2026: first draft.
- v0.2, 30 October 2026: accepted; the three layers; Lindal's answer unchanged as the test, with
  no full regression; the inserted 100 mbar level ruled.
- v0.3, 5 November 2026: §5 rules on the pre-execution report; acceptance 4 restated against
  Kliore's precision; acceptance 5, the hydrostatic check, added.
- v0.4, 5 November 2026: closed; the drift is Kliore's molar mass times gravity, 0.6 percent below
  ours, likely 2.2 amu.
