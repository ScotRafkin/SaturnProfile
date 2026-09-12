# REPORT 01, Step 8. The composition tool for Lindal, kind C

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 11 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.16, Step 8, against
`SPEC_00_Architecture_and_Data_Files.md` v0.8. Status: reported, awaiting review.
**Six of seven acceptance checks pass.** The one that fails is an acceptance number that does
not follow from the rule the same paragraph declares; finding 1.

**The working tree is not committed**, per the commit on acceptance rule.
`lindal_composition.nc` carries `-dirty` and is to be rebuilt after the acceptance commit.

**Step 7 closed first.** Committed as `c618e80`, the product rebuilt clean, and the review's
`anchor_latitude` fix committed as `8fe73ad`: under `anchor_rule = "latitude"` the anchor
latitude is now unioned into the march nodes, so the rule that the anchor radius is never
interpolated holds for every anchor rule. Verified on a deliberately coarse 1 degree caller
grid: the anchor resolves to 3.0e-08 m, matching a grid that already contains the target.

---

## 1. What was built

**`src/casspian/tools/composition/build_composition.py`**, entry point
`casspian-composition-lindal`, and the `[composition]` section of `lindal_build.toml`.
Product: `occul_data/lindal/lindal_composition.nc`, kind C, 66 levels with a `species` group.

Ammonia is taken from the Table I column in the raw bundle, filled linearly in pressure through
the nine tabulated values, extrapolated beyond them, and clamped at zero. It is assigned first
and the remainder is split 0.94 H2 to 0.06 He, so the three mole fractions are `x`, `s(1-x)` and
`(1-s)(1-x)` and sum to one identically rather than by normalization.

## 2. Decisions

1. **The split ratio is declared in the control file and then checked against the source.**
   SPEC_01 has the control file declare 0.94 and 0.06 as the ratio to hold. That is a choice
   about dividing the remainder and is not the same statement as the source's "94 percent H2";
   here they coincide. The tool refuses if the declared H2 share disagrees with the
   `h2_fraction` the raw bundle transcribed from p. 1137, so the coincidence cannot drift apart
   unnoticed.

2. **`x_H2_uncertainty` is read from the raw bundle, not declared in the control file.** SPEC_01
   lists it among the things Step 8 writes, and 0.03 is a value the source states (p. 1137,
   "94 +/- 3%"). A physical value the source states is data, so it comes from the transcription
   rather than from a control file, which principle 2 keeps free of them.

3. **A per level `nh3_provenance` flag accompanies the variable level `provenance`.** SPEC_00
   section 5 puts one `provenance` on a variable, but the ammonia column is `measured` at nine
   levels, `interpolated` at one interior gap and `extrapolated` elsewhere. The variable carries
   the dominant rule and the companion flag carries the per level truth, which is the same
   pattern kind W uses for `value_provenance`.

4. **The absent dimension attribute is `latitude_planetocentric_absent_meaning`.** SPEC_01
   Step 8 writes `latitude_absent_meaning`; SPEC_00 section 5 defines the attribute as
   `<dimension>_absent_meaning`, and the dimension of kind C is `latitude_planetocentric`. The
   schema-correct name is written. Step 9 uses the same shorthand and will want the same
   reading.

5. **Recorded paths are relative to the repository root.** The first build wrote
   `master_table` and every `input_hashes` entry as an absolute Windows path, which puts one
   machine's directory layout into a data product. SPEC_00 section 6.2 says "path under [D0]".
   Corrected before the acceptance was run.

## 3. Findings

### 1. The stated clamp level does not follow from the declared rule

SPEC_01 Step 8 declares the fill rule (linear in pressure, `interp1d` with
`fill_value="extrapolate"`, clamped at zero) and then states as acceptance that "NH3 is zero at
and above 794.33 mbar (the clamp)". Those two cannot both hold.

The first two tabulated points are (831.76 mbar, 2.60 ppm) and (870.96 mbar, 5.10 ppm). Linear
extrapolation below the first crosses zero at **790.992 mbar**, which is **3.34 mbar below**
794.33. So:

| level | value |
|---|---|
| 831.76 mbar (first tabulated) | 2.600 ppm |
| **794.33 mbar** | **0.213 ppm, not zero** |
| 758.58 mbar | 0.0, clamped |
| every level above 758.58 | 0.0, clamped |

Interpolating linearly in log pressure instead, which is the other natural reading, moves the
crossing only to 792.86 mbar, still below 794.33. No variant of the declared rule makes 794.33
the clamp level.

**Proposing the acceptance be restated** as "NH3 is zero at and above 758.58 mbar, the first
grid level at or below the zero crossing at 790.99 mbar". Nothing in the product changes; the
number in the specification does. It is worth asking where 794.33 came from, since it is a grid
level and not a round number, so it was computed by something.

The rest of the ammonia column reproduces exactly: **15.931 ppm** at the 1047.13 mbar interior
gap against the specified 15.9 +- 0.1, and **79.253 ppm** at the 1298.48 mbar end gap against
the 79.3 Phase 1 found.

### 2. `is_polar` is a property of the molecule but the master table carries it per set

`data_static/species_master.toml` carries `is_polar` inside each refractivity set, and the
`lindal1985` entry for NH3 omits it, having only the value, status and source. Reading it from
that set alone and defaulting to false would have written **`is_polar = false` for ammonia**,
which is a statement the table never made and which is wrong: ammonia is polar, and the `modern`
set says so.

The tool now looks the property up across every set in the table, uses the value if the sets
agree, and refuses if they disagree, since a molecular property cannot differ by which
refractivity set is chosen. `is_polar` now reads `[0, 0, 1]` for H2, He, NH3.

**Proposing that `is_polar` move out of the per set tables and sit beside `molar_mass`**, which
is where the table already keeps molecular properties. `temperature_dependence` is the harder
case and should stay per set, since it describes the refractivity value rather than the
molecule.

### 3. What ammonia actually does in this reduction

Worth stating because it bounds finding 1. The `lindal1985` set carries NH3 refractivity
**zero**, on the Fig. 3 caption reading that Lindal computed density from H2 and He only. So
ammonia enters the reduction through the mean molar mass and nowhere else. Between a level with
no ammonia and the deepest level, where it reaches 79.3 ppm, the mean molar mass moves from
2.135083 to 2.136264 amu, **0.055 percent**. The clamp level dispute of finding 1 concerns
levels where the fill is at most 0.2 ppm, which is a part in 10^6 of the mean molar mass. It
should be settled for correctness, not because a number depends on it.

## 4. Acceptance results

Run by `reports/step8/accept_step8.py`; full output in `reports/step8/output.txt`.

| Check | Measured |
|---|---|
| Mole fractions sum to one at every level to 1e-12 | **Pass. 1.110e-16** over 66 levels. Exact by construction, not by normalization. |
| NH3 zero at and above 794.33 mbar | **Fail as stated.** Zero at and above **758.58 mbar**; **0.213 ppm** at 794.33; crossing at **790.992 mbar**. Finding 1. |
| NH3 at the 1047.13 mbar interior gap, 15.9 +- 0.1 ppm | **Pass. 15.931 ppm.** |
| NH3 extrapolated at 1298.48 mbar, reported | **79.253 ppm** (Phase 1 found 79.3). |
| Mean molar mass with NH3 zero, 2.1351 amu | **Pass. 2.135083 amu.** |
| Mean molar mass at the deepest level, reported | **2.136264 amu** at 1298.48 mbar with 79.3 ppm ammonia, 0.055 percent above the dry value. |
| Mean molecular refractivity at a zero-ammonia level, to six figures | **Pass. 4.836272e-30 m3**, against `(0.94 x 136 + 0.06 x 35) x 1e-6 / n_Loschmidt` = 4.836272e-30, relative difference **0.0e+00**. Per species: H2 5.061821e-30, He 1.302675e-30, NH3 0. |
| `read` as kind C succeeds | **Pass.** `composition_role = "reduction"`, `vertical_coordinate = "pressure_Pa"`, `latitude_planetocentric_absent_meaning = "point"` with `latitude_planetographic_deg = 36.3`. |
| `read` refuses a copy with one mole fraction perturbed by 1e-6 | **Pass.** `mole fractions depart from one by 1.000e-06, which exceeds the tolerance 1.0e-09 of SPEC_00 section 6.2`. |
| Beyond the specification: the species group | **Pass.** All seven section 6.2 fields present, `master_table = "data_static/species_master.toml"` with its SHA-256, `is_polar = [0, 0, 1]`. |

Six of seven pass. No em dash or en dash appears in any file written in this step.

## 5. Next step

Step 9, stage two of the Lindal tool and the manifest, is not started and is not to be begun
until this report has been reviewed. It is the step that completes `occul_data/lindal/` and the
first that reads this file, kind G, kind R and kind W together.

---

## 6. Changes applied after review

Added 11 September 2026, after `REVIEW_01_step8.md`. Seven of seven now pass.

**1. The ammonia rule is restated in three parts, and finding 1 dissolves.** The review answers
the question the report should have asked: is zero aloft Lindal's statement or our assumption?
It is our assumption. Table I is blank above 831.76 mbar, and blank means not measured. What
justifies zero there is physics rather than the table, and the file now says so in an
`nh3_aloft_assumption` attribute: the troposphere at those levels is saturated with ammonia,
and the saturation mixing ratio over NH3 ice on Lindal's own temperatures falls from a few ppm
near 0.8 bar, consistent with the tabulated 2.6 ppm at 831.76, to below 1e-7 by 0.5 bar.

So there is no upward extrapolation and no clamp. Ammonia is **exactly zero** above the highest
tabulated level, flagged `assumed`, and the acceptance number holds exactly: 794.33 mbar is the
first grid level above 831.76, and it reads 0.0. The 0.213 ppm of finding 1 was an artifact of
extrapolating a measurement into a region where nothing was measured, then clamping the result.
The number 794.33 was never a computed crossing, which is what made it puzzling.

I should record what I got wrong there. I treated the mismatch as an arithmetic discrepancy
between a rule and a number, checked both readings of the rule, and proposed changing the
number. The rule was the thing at fault, and the question that would have found it, what does a
blank cell in Table I mean, is one I did not ask.

**2. `is_polar` moved to its own table.** The author's edit puts it in an `[is_polar]` table
beside `molar_mass` and removes the per set copies; the tool reads it from there and refuses a
species the table is silent about, rather than looking across sets as the interim version did.
`temperature_dependence` stays per set, since it describes the refractivity value and not the
molecule. The same edit promotes the `lindal1985` NH3 entry to `verified_reading` with the
caption and footnote quoted.

**3. The per level flag gains `assumed`**, so `nh3_provenance` now carries
`measured interpolated assumed extrapolated` on values 0, 1, 2, 3, and all four codes appear in
the product.

**Acceptance after the changes.**

| Check | Measured |
|---|---|
| Mole fractions sum to one | **1.110e-16** over 66 levels |
| Ammonia zero above the highest tabulated level, flagged `assumed` | **Pass.** 794.33 mbar reads 0.0 exactly and is flagged `assumed`; 55 levels aloft, all zero, all assumed; 831.76 mbar reads 2.600 ppm, flagged `measured` |
| Interior gap at 1047.13 mbar | **15.931 ppm** |
| Extrapolated at 1298.48 mbar | **79.253 ppm** |
| Mean molar mass, ammonia zero | **2.135083 amu** |
| Mean molar mass, deepest level | **2.136264 amu** |
| Mean molecular refractivity | **4.836272e-30 m3**, relative difference 0.0e+00 |
| `read` as kind C, and the refusal | **Pass**, refused at 1.000e-06 against the 1e-9 tolerance |
| Species group | **Pass**, `is_polar = [0, 0, 1]`, `master_table = "data_static/species_master.toml"` with its hash |

Seven of seven.


---

## 7. Addendum after closure: the v0.19 amendment

Added 12 September 2026, per `REVIEW_02_step3.md` findings 4 and 5 and SPEC_01 v0.19 Step 8.
SPEC_02 Step 3 found two things the reduction needs from this file that it did not carry.

**What changed.**

1. **The closure is declared.** The composition file now carries `closure_rule =
   "share_of_remainder"` and `closure_species = "H2 He"` (SPEC_00 section 6.2 v0.13), with the
   share species first. `refrac` reads them and checks them against the values; it no longer
   infers the structure.
2. **Unstated per molecule uncertainties are NaN.** `refractivity_uncertainty_m3` is now NaN for
   all three species, with `uncertainty_method` saying NaN is written where the master table
   states no uncertainty. It used to be 0.0.

**Two decisions.**

1. **The master table says `nan`, not the tool.** `data_static/species_master.toml` wrote an
   unstated uncertainty as `uncertainty_e6 = 0.0`, with a comment saying it was not stated. A
   tool reading that table cannot tell a stated zero from a placeholder, so the fix belongs in
   the table. Every `uncertainty_e6 = 0.0` placeholder is now `nan`, a TOML float, with its
   comment kept. That is two entries in `lindal1985` (H2, He) and four in `modern`. The `modern`
   set feeds no product yet; it was changed so the table follows one rule. A species with no
   `uncertainty_e6` key at all (the `lindal1985` NH3 entry) also becomes NaN in the tool, so
   leaving a value out cannot turn it into a zero. A finite value, zero included, still passes
   through as stated.
2. **The closure attributes are constants of the tool, not control file keys.** The tool always
   assigns ammonia first and splits the remainder between H2 and He. The declaration describes
   what the code does, so it is written from the code, like `remainder_split` beside it. The
   split ratio stays a control file choice, as before.

**Acceptance.** `reports/step8/accept_step8.py`, nine of nine.

| Check | Measured |
|---|---|
| 1 to 7, the original Step 8 checks | **Pass, every number unchanged.** The only line of output that differs from the pre-amendment run is `master_table_hash`, which moved from `dc667eaa...` to `a1a16e33...` because the table itself changed. |
| 8. The closure declaration reads back and agrees with the values | **Pass.** `closure_rule = 'share_of_remainder'`, `closure_species = 'H2 He'`. `x_H2 / (x_H2 + x_He)` runs 0.94 to 0.9400000000000001, spread 1.1e-16. |
| 9. `refractivity_uncertainty_m3` is NaN for all three species | **Pass.** `[nan, nan, nan]` for H2, He, NH3. |

**Nothing else was rebuilt.** Of the files under `occul_data/lindal/`, only the composition
hashes the master table or depends on the composition. The others still carry `ece58d2`.

**The new hash.** The composition was rebuilt on a clean tree after the amendment was committed.
This row supersedes the `lindal_composition.nc` row of REPORT_01_step9 section 6, and the SPEC_02
Step 1 acceptance reads it from here.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_composition.nc` | composition | `a1250b60c4b68c42109bb1fbe22244b218af9866aa1e3007f442a9d2ba26654f` |
