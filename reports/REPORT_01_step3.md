# REPORT 01, Step 3. The gravity and rotation tools, kinds G and R

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.5, Step 3, against
`SPEC_00_Architecture_and_Data_Files.md` v0.6. Status: reported, awaiting review.

**The working tree is not committed**, per the commit on acceptance rule. Both products carry
`48969a72fd03f2224d5bdff7ce2a2baa857fe608-dirty` and are to be rebuilt after the acceptance
commit.

**Step 2 closed first.** Step 2 was committed as `48969a7` and `casspian-lindal-raw` was rerun,
so `lindal_raw.nc` now carries the clean commit `48969a7` with no suffix, and its `raw_sources`
picks up the amended `notes.md`. The eight Step 2 acceptance checks were rerun against the
rebuilt bundle and all still pass.

---

## 1. What was built

**`src/casspian/tools/gravity/control.py`.** Control file reading: loads one section, checks
its keys against the set the tool knows, refuses an unknown key as an error rather than a
warning (SPEC_00 section 7), refuses a bare number in a control file (principle 2), and
resolves relative paths against the directory holding the control file, as section 7 requires.

**`src/casspian/tools/gravity/build_gravity.py`**, entry point `casspian-gravity-file`. Writes
kind G from a `data_static` harmonic transcription. **`build_rotation.py`**, entry point
`casspian-rotation-file`. Writes kind R from a rotation entry. Both entry points are declared
in `pyproject.toml`.

**`occul_data/lindal/lindal_build.toml`**, the control file, with `[gravity]` and `[rotation]`
only. It names the static source, the output, and the prefix, and nothing else; the later steps
add their own sections. No physical value appears in it, and the tool refuses one if it does.

Products: `occul_data/lindal/lindal_gravity.nc` from `null1981` and
`occul_data/lindal/lindal_rotation.nc` from `system_iii`.

## 2. Decisions and observations

1. **`fitted` and `status` are kept apart.** The static files carry both. `fitted` says how the
   number was obtained and becomes the `J_status` flag variable (0 fitted, 1 assumed, 2
   constrained, CF style with `flag_values` and `flag_meanings`). `status` says whether a second
   reader has checked the transcription against the paper, reads `verified` throughout, and
   becomes a per-variable attribute, which is what Step 3 asks for. They answer different
   questions and are not merged.

2. **The netCDF asserts the convention code; the source's prose is kept beside it.** The file
   carries `harmonic_convention = "CASSPIAN-J1"` and, in
   `harmonic_convention_source_statement`, the sentence the transcription states:
   `V = -(GM/r)[1 - sum J_2n (R/r)^2n P_2n(sin phi_c)], g_N positive inward, phi_c
   planetocentric`. I compared that against the definition of `CASSPIAN-J1` in `lib.schema` term
   by term and they agree. **No machine check is possible here**: the tool cannot verify prose
   against a code, so the assertion that a given transcription obeys `CASSPIAN-J1` rests on a
   human having read both. Recording the prose in the file is what makes that check repeatable
   by the next reader.

3. **`angular_rate_rad_s` is recomputed, not transcribed**, as the static file itself instructs.
   The transcribed value 1.63785e-4 differs from `2 pi / 38362.4` by 9.823e-12 rad/s, being
   rounded for a human reader. The transcribed number is kept as the attribute
   `angular_rate_transcribed_rad_s` with a note saying which one to use, so the difference is
   visible rather than lost.

4. **`[check_values]` in `null1981.toml` is not written into kind G.** Those numbers are
   Lindal's Table II gravities and the trap about the `(degree + 1)` coefficient. They are the
   acceptance aid for Step 4, not part of the harmonic set, and kind G has no place for them.
   SPEC_00 section 2.1 forbids hard-coding any value that appears in a static file, so the
   Step 4 acceptance script must read them from `data_static/harmonics/null1981.toml` rather
   than repeating them.

5. **`uncertainty_scaling` is written only when the source states one.** SPEC_00 section 6.4
   lists it among the kind G globals, but it is a property of a particular solution: `iess2019`
   declares three times formal, `null1981` declares nothing. It is written when present and
   omitted otherwise, and the registry does not require it.

6. **Control file parsing lives in the tool, not in `lib.control`.** SPEC_00 section 3.1 gives
   `lib.control` the reduction manifest and the run namelist, which are Step 9 and later.
   Building the general parser now would be implementing ahead, so `tools/gravity/control.py`
   does the narrow job and will be absorbed when `lib.control` arrives.

## 3. Defects and questions raised

1. **`J6` in `null1981.toml` carries `uncertainty = 0.0`, which claims perfect knowledge of an
   assumed value.** The entry is `fitted = false`, and the source's own words are that
   `J6 = 84 x 10^-6` comes from interior models and "is used for all calculations in this
   paper." Its uncertainty is unstated, not zero. SPEC_00 section 5 makes NaN the way a file
   says no information at that location, and 0.0 says something quite different: a later error
   propagation reading this file would treat J6 as exact and would understate the uncertainty
   on `g`. The tool transcribes what the file says rather than patching it, so
   `J_uncertainty` currently reads `[1.8e-05, 3.8e-05, 0.0]`. **Proposing that the 0.0 becomes
   NaN in the static file**, which is a one line change there and needs no code change here.
   Raised rather than acted on, since `data_static` is the author's by hand.

   *Resolved by `REVIEW_01_step3.md`: `null1981.toml` now carries `nan`, and SPEC_00 v0.7
   section 5 makes it a rule that an uncertainty the source does not state is NaN and never 0.0,
   in transcriptions as much as in netCDF. After the rerun `J_uncertainty` reads
   `[1.8e-05, 3.8e-05, nan]`.*

2. **The provenance vocabulary has no term for a coordinate that is a label.** `degree` holds
   2, 4, 6: an index, not a measurement. SPEC_00 section 5 requires a `provenance` on every
   variable and offers eight values, none of which fits an index. It is written as `assumed`,
   which is the nearest, and is wrong in spirit. The same will arise for `species_name` in kind
   C and `surface` in kind D. A ninth value such as `index`, or an exemption for coordinate
   variables that are labels, would settle it. Not urgent, and no number depends on it.

3. **`J` is written with `provenance = "inferred"`** on the reading that a harmonic coefficient
   fitted to tracking data is deduced from other measurements rather than measured directly.
   `derived` would also be defensible. Confirm which the author intends, since every later kind
   G file will follow whichever is chosen.

## 4. Acceptance results

Run by `reports/step3/accept_step3.py`; full output in `reports/step3/output.txt`.

| Check | Result |
|---|---|
| `lindal_gravity.nc` reads back with `J` = [16479e-6, -937e-6, 84e-6] on `degree` = [2, 4, 6] | **Pass**, by exact equality: `[0.016479, -0.000937, 8.4e-05]`. |
| `J_status` = [fitted, fitted, assumed] | **Pass.** Stored `[0, 0, 1]`, which resolves through `flag_meanings` to `['fitted', 'fitted', 'assumed']`. |
| `normalization_radius_m` = 6.0e7 | **Pass.** `60000000.0`. |
| `GM_m3s2` = 3.7929085e16 | **Pass.** `3.7929085e+16`, the planet GM. The system GM `3.7938544e+16` is an attribute, not a variable, as Step 3 requires. |
| A second run pointed at `iess2019` gives six degrees to 12 and `normalization_radius_m` = 6.033e7 | **Pass.** `degree = [2, 4, 6, 8, 10, 12]`, `normalization_radius_m = 60330000.0`, all six `fitted`, `uncertainty_scaling` carried through, `GM_m3s2 = 3.7931206234e+16`. The run used a generated control file under `reports/step3/products/`, so nothing outside that directory was written. |
| `lindal_rotation.nc` gives `period_s` = 38362.4 and `angular_rate_rad_s` = 2 pi / 38362.4 to double precision | **Pass.** `0.00016378499017735037`, bit identical to `2 * math.pi / 38362.4`. |
| The reader refuses a G file whose `harmonic_convention` is a code it does not know | **Pass.** With the attribute edited in place to `CASSPIAN-J2` on a copy: `harmonic_convention is 'CASSPIAN-J2', which this reader does not know (SPEC_00 section 6.4). The codes it knows are: CASSPIAN-J1.` |
| Beyond the specification: an unknown control file key is an error | **Pass.** A `[gravity]` section carrying `normalization_radius_m` is refused by name, citing section 7, and the message lists the keys the section accepts. This also demonstrates the principle 2 check, since the offending key is a physical value. |

Five of five pass. No em dash or en dash appears in any file written in this step.

## 5. Next step

Step 4, `lib.gravity`, is not started and is not to be begun until this report has been
reviewed. It is the step the `(degree + 1)` trap in `null1981.toml` is waiting for, and its
acceptance numbers are to be read from that file rather than repeated in the script.
