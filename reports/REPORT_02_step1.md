# REPORT 02, Step 1. `lib.control`: the reduction manifest and input loading

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.2, Step 1, against
`SPEC_00_Architecture_and_Data_Files.md` v0.11 and the closed `SPEC_01_Lindal_Tool_Chain.md` v0.18
(commit `ece58d2`). Status: reported, awaiting review. **Six of six acceptance checks pass.**

**On the status of SPEC_02.** Its header still reads "draft for author markup; no step to begin
before the draft is accepted", and `STATE.md` had no SPEC_02 rows. Step 1 was begun on the
author's instruction that the specification was ready. Section 3 below carries findings from a
read-only check of the draft's numbers against the real files, for the markup; none of them
bears on Step 1 or Step 2, and two of them would fail Step 3's acceptance as written.

**The working tree is not committed**, per the commit on acceptance rule. This step writes no
data product, so the directory sweep after the acceptance commit will find nothing to rebuild.

---

## 1. What was built

**`src/casspian/lib/control.py`**, the general control file reader SPEC_00 section 3.1 names.

- `ControlFileError`, `load_section` and `reject_physical_values`, absorbed from
  `tools/gravity/control.py` with their behavior unchanged, including the named-key path
  resolution. The old module is **deleted, not shimmed**: `casspian.tools.gravity.control` is no
  longer importable, so there is one control reader and not two. The five tools that imported it
  and the one SPEC_01 acceptance script that did (`reports/step3/accept_step3.py`) now import
  `lib.control`.
- `read_reduction_manifest(path)`, returning a frozen `ReductionManifest`: every section and key
  checked against the SPEC_00 section 7.1 vocabulary, types and ranges checked, paths resolved
  against the manifest's directory, the manifest text kept verbatim with its SHA-256.
- `load_reduction_inputs(manifest)`, returning a frozen `ReductionInputs`: the six inputs read
  under their kinds with `lib.io.read`, pulled into memory with their files released, each with
  the SHA-256 and `casspian_git_commit` it carried on disk, and checked against each other before
  any arithmetic.

## 2. Decisions

1. **An unknown key is refused either way; what the message calls it depends on the key.** The
   manifest legitimately carries numbers, such as `pressure_Pa`, `convergence_m` and
   `max_iterations`, because those keys are declared choices the vocabulary knows. So "holds a
   number" cannot be the test for a physical value. The test is: a key the vocabulary does not
   know, whose name ends in a unit suffix and whose value is a number, is refused as **a physical
   value in a control file**; any other unknown key is refused as **unknown**. The acceptance's
   `polar_radius_m = 54438000` needed the first message and a stray text key gets the second.
   The unit suffixes are `lib.schema.UNIT_SUFFIXES`, the same list the uncertainty naming rule
   uses, so there is one definition of what a unit suffix is.

2. **The manifest record has no radian form of the tolerance.** It returns
   `fixed_point_tolerance_deg` in degrees as written. SPEC_02 section 0 has the conversion happen
   once, at the `refrac` boundary; a convenience property here would be a second place that
   converts, which is the thing that rule exists to prevent.

3. **The optional sections get only the defaults SPEC_00 section 7.1 declares.** `figures` and
   `enabled` default to false. When either is switched on, the keys that feature needs
   (`format` and `dpi`; `wind_draws` and `wind_draw_seed`) are required, because section 7.1
   shows them in its example but declares no default for them, and inventing one would be a
   choice made silently in the parser.

4. **A manifest may name three of the four anchor rules.** `lib.geoid.wind_geoid` also knows
   `"latitude"`, but that rule needs an anchor latitude and section 7.1 has no key to carry one,
   so the parser refuses it with that reason rather than failing later inside the march.

5. **Manifest paths must carry the profile prefix.** SPEC_00 section 8 lists a prefix refusal
   only for namelists; section 2.2 says every file under a profile directory carries the prefix.
   Enforcing it on manifests is stricter than the section 8 list, in the direction the Step 1
   review of SPEC_01 called the floor rather than the ceiling. Flagged so it is deliberate.

6. **Two load-time checks beyond Step 1's list.** The geoid anchor surface must be a surface of
   the geodesy file, and the named anchor quantity must be a variable of it. Step 2 reads both,
   and a manifest that passed Step 1 and then failed inside the march would be a worse failure.

## 3. Findings for the SPEC_02 markup

A read-only check before building, of Step 1's preconditions on the real inputs and of the
draft's stated numbers by direct arithmetic on the same files.

**Step 1's preconditions all hold.** Every product still matches the SHA-256 recorded in
REPORT_01_step9 section 6 and carries the clean commit `ece58d2`. The anchor isobar, 10000 Pa, is
both a geodesy surface and tabulated level 29 of the thermo file, with `h_ref` = 90,000 m. The
composition and thermo levels are bit-identical. The wind's reference level is a node of its
pressure grid, and its rotation system and rate match kind R exactly.

1. **The step working directories collide with SPEC_01's.** SPEC_02 section 0 puts acceptance
   scripts under `reports/step<N>/`, but `reports/step1/` to `reports/step9/` already exist and
   hold SPEC_01's scripts, which the regression still runs. This step uses `reports/step02_1/`,
   which the existing `reports/step*/` ignore rule already covers. **Proposing SPEC_02 section 0
   name `reports/step02_<N>/`.**

2. **Step 3's composition uncertainty term is 0.02332, not 0.0227, at every level.** The draft
   states the fractional uncertainty of N "at every level equals the composition term, 0.0227 +-
   0.0002", with the term `0.03 x (R_H2 - R_He) / R_bar`. By direct arithmetic on the species
   group:

   | quantity | value |
   |---|---|
   | `R_H2` | 5.061821e-30 m3 |
   | `R_He` | 1.302675e-30 m3 |
   | `R_bar`, dry | 4.836272e-30 m3 |
   | `0.03 x (R_H2 - R_He) / R_bar` | **0.02332** |
   | the same from the STP values, 0.03 x (136 - 35) / 129.94 | 0.02332 |

   It is also **constant across all 66 levels, to 3.5e-18**, because ammonia dilutes the H2 and
   He shares by the same factor and so cancels between numerator and denominator. The acceptance
   as written therefore fails at every level. The denominator that would give 0.0227 is
   4.968e-30 m3, which is not any `R_bar` the file contains, so the 0.0227 does not appear to come
   from a different but defensible convention. **Proposing 0.0233 +- 0.0002.**

3. **Step 3's 1 bar refractivity assumes zero ammonia where Table I tabulates 10.9 ppm.** The
   draft states N = 2.59859e-4 at 1 bar, "(R_bar = 4.836272e-30 m3 at zero NH3)". The arithmetic
   is right for the dry `R_bar`: 5.3731235e25 times 4.8362722e-30 is 2.598589e-4. But Table I
   tabulates NH3 at 1000.00 mbar, 10.9 ppm, so the product will write **2.598560e-4**, 1.1e-5
   relative below the stated value. **Proposing** either quoting 2.598560e-4 at 1 bar, or
   moving the zero-ammonia check to a level at or above 794.33 mbar, where ammonia is exactly zero
   by the v0.17 rule and the dry `R_bar` is the right one.

4. **The rest of Step 3's numbers check out.** Number density at 1 bar 5.373124e25, at the top
   level 1.04441e22 and at the bottom level 6.43287e25 m-3, all matching to the figures given. The
   bottom level's `R_bar` is `1 - 79.2528e-6` times the dry value, of which the draft's
   `1 - 79.3e-6` is the rounding. I first reported that ratio as a discrepancy, from a check that
   used the 1 bar level as the dry reference; it is not one, and finding 3 is why.

## 4. Acceptance results

Run by `reports/step02_1/accept_step02_1.py`; full output in `reports/step02_1/output.txt`.

| Check | Measured |
|---|---|
| `lindal_reduction.toml` loads with all six inputs and the hashes of REPORT_01_step9 section 6 | **Pass.** All six SHA-256 and the manifest's match the recorded table; all six carry commit `ece58d2`. Anchor isobar 10000 Pa, anchor rule `mean_polar_radius` on the 10000 Pa surface by `radius_polar_m`; tolerance 1e-06 returned in degrees; composition loaded as a DataTree and the other five as Datasets; both optional sections at their declared defaults. |
| A manifest naming a `-dirty` product is refused with the file named | **Pass.** `lindal_wind.nc carries casspian_git_commit = 'ece58d2...-dirty'. SPEC_00 section 8: a file offered as an input to refrac may not carry -dirty.` |
| `pressure_Pa = 9000` refused as neither a geodesy surface nor a tabulated level | **Pass**, with both failures named in one message. The manifest itself parses, 9000 being a positive pressure; the refusal is a load-time rule because only the inputs can say it is wrong. |
| An extra `polar_radius_m = 54438000` under `[geoid]` refused as a physical value | **Pass**, at parse time: `[geoid] polar_radius_m = 54438000 is a physical value in a control file.` |
| Beyond the specification: other parse-time refusals | **Pass.** An unknown text key is refused and named as unknown, not physical; `wind = "../lindal_wind.nc"` is refused as outside the manifest's directory; `anchor_rule = "latitude"` is refused with the reason. |
| Beyond the specification: the old module absorbed, not shimmed | **Pass.** `casspian.tools.gravity.control` is not importable; `load_section` on `lindal_build.toml [gravity]` still resolves `source` to an absolute path that exists. |

Six of six pass. The ten SPEC_01 suites were rerun after the control reader moved, since five
tools and one of those suites import it: all pass, 6, 14, 8, 5, 7, 7, 6, 8, 7 and 6, including
Step 3's suite, which exercises the unknown-key refusal through a tool. No product under
`occul_data/` was rebuilt, so all still carry `ece58d2`. No em dash or en dash appears in any file written in this step.

## 5. Next step

Step 2, the frozen `phi_c` and `r0`, is not started and is not to be begun until this report has
been reviewed. It is the first step that reads the loaded inputs, and the first whose numbers
depend on the `anchor_rule` the manifest now carries.
