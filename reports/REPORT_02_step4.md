# REPORT 02, Step 4. The kind N product and `casspian-refrac`

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 12 September 2026. Specification: `SPEC_02_Refrac_and_Diagnostics.md` v0.6, Step 4,
against `SPEC_00_Architecture_and_Data_Files.md` v0.13 and `SPEC_01_Lindal_Tool_Chain.md` v0.19.
Status: accepted with one change (REVIEW_02_step4), applied; see section 7. **Nine of nine acceptance checks pass** after the change.

**The working tree is not committed**, per the commit on acceptance rule.
`occul_data/lindal/lindal_refractivity.nc` was written by this step from an uncommitted tree, so
it carries `fc40d22...-dirty`. The directory sweep after the acceptance commit will rebuild it,
and nothing else.

## 0. Work done before this step, per REVIEW_02_step3

The review's order of work was followed and each stage committed:

| Commit | What |
|---|---|
| `fc7cef3` | SPEC_00 v0.13, SPEC_01 v0.19, SPEC_02 v0.6, and the Step 3 review |
| `6ee8690` | SPEC_01 v0.19 Step 8 amendment. The composition tool declares its closure and writes NaN per molecule uncertainties. The not stated placeholders in `species_master.toml` became `nan`. |
| `31c30ef` | The rebuild of `lindal_composition.nc` on a clean tree, hash `a1250b60...`, recorded in REPORT_01_step8 section 7. Step 8 passes 9 of 9 and its other numbers are unchanged. |
| `fc40d22` | SPEC_02 Step 3 with the v0.6 changes. 9 of 9; addendum in REPORT_02_step3 section 7. |

---

## 1. What was built

**`src/casspian/refrac/product.py`**: `build_product(manifest_path)` and `main`, registered as
**`casspian-refrac`** in `pyproject.toml`. One call reads the manifest, loads the six inputs,
freezes the anchor, reduces, computes the anchor rule spread, and writes kind N to the manifest's
`[output] product`.

The file, as written for Lindal:

- **On `level`.**
  - `radius_m` is the coordinate, with `positive = "up"`.
  - Variables: `height_above_anchor_isobar_m`, `number_density_m3`, `refractivity`,
    `mean_refractivity_m3`, `mean_molar_mass_kg_mol`.
  - Six companions, one for `radius_m` and one for each variable.
- **Six scalars, each with a companion:** `latitude_planetocentric_deg`, `psi_deg`,
  `latitude_planetographic_deg` (with its swath, swath source, value source and definition from
  kind T), `anchor_isobar_pressure_Pa`, `anchor_isobar_radius_m`, `anchor_isobar_height_m`.
- **Groups.**
  - `inputs/thermo`, `inputs/composition` with `inputs/composition/species`, `inputs/geodesy`,
    `inputs/gravity`, `inputs/rotation`, `inputs/wind`.
  - `manifest`, with `text` and `sha256`.
  - `reduction_record`, 58 attributes (section 3).
- **Global `input_hashes`** over the six inputs and the manifest, paths relative to the
  repository root.

**`src/casspian/lib/schema.py`**, the kind N registry and the naming rule:

- `mean_refractivity_m3` and `mean_molar_mass_kg_mol` are added, each requiring its companion
  (SPEC_00 section 6.7 v0.13).
- `height_above_anchor_isobar_m` and all six scalars now require companions. Section 6.7 lists
  "scalars with uncertainties" and the registry had marked only two.
- `_kg_mol` is added to `UNIT_SUFFIXES`. Without it the naming rule gave
  `mean_molar_mass_kg_mol_uncertainty` instead of section 6.7's
  `mean_molar_mass_uncertainty_kg_mol`.

**`src/casspian/lib/io.py`**: `read` of kind N checks the hashes (decision 1).

**`src/casspian/refrac/reduce.py`** gains the `height_above_anchor_isobar_uncertainty_m`
companion. It carries the two height terms and uses the Step 3 rule that `h - h_ref` has no
height term at the anchor level. It is NaN for Lindal, both terms unstated.

## 2. Decisions

1. **How `read` detects an edited `input_hashes`.** A hash recorded once cannot catch its own
   edit, so `refrac` records each input's SHA-256 twice: in the global `input_hashes` and as
   `input_sha256_<input>` in `reduction_record`. `read` of kind N refuses when a recorded hash is
   no longer listed in `input_hashes`, naming the input, or when the entry count is not seven.
   It **warns**, and does not refuse, when an entry no longer matches the file now on disk at
   that path, as SPEC_00 section 8 says. The embedded copy is authoritative. This catches an
   accidental or careless edit, not a deliberate forgery of both attributes; the embedded
   groups are the audit copy for that.
2. **Copied companions keep their declared kind.** `latitude_planetographic_uncertainty_deg` is
   copied from kind T as declared: 0.2, `range`, `uncertainty_method` saying it is copied and
   not converted. The same holds for `anchor_isobar_pressure_uncertainty_Pa` and
   `anchor_isobar_height_uncertainty_m` (NaN, `stated`). The SPEC_00 section 5 v0.13 conversion
   applies to derived companions. A copied value is not derived, and converting it here would
   lose what the source declared. **The propagated companions are all `1sigma`** and carry
   `uncertainty_terms_included`, `uncertainty_terms_unstated` and
   `uncertainty_kind_conversions`. The copied ones carry none of the three, since no terms were
   combined.
3. **List attributes.** Term labels are space separated; they contain no spaces. Conversion
   texts are newline separated, as `input_hashes` is. An empty list is an empty string, so every
   propagated companion carries the same three attributes.
4. **Provenance of the scalars.**
   - `latitude_planetocentric_deg`, `psi_deg`, `anchor_isobar_radius_m`: `derived`.
   - `latitude_planetographic_deg`: `measured`, with its `value_source` copied (the Fig. 4
     label).
   - `anchor_isobar_pressure_Pa`: `index`, since it is the manifest's choice of level, a label.
   - `anchor_isobar_height_m`: `measured`, the tabulated `h_ref`.
5. **The anchor rule spread is taken at the frozen `phi_c`**, the latitude held, as "two extra
   marches". The fixed point is not rerun under each rule. It is recorded as
   `anchor_rule_spread_r0_north_pole_m` and `..._south_pole_m` with a note that it is a choice
   and enters no companion.
6. **Radians become degrees once, at the write**, for the stored scalars and the iterates. The
   partials and terms in `reduction_record` keep the units their names carry (`_rad`,
   `_m_per_rad`), exactly as Step 3 computed them.
7. **The composition's nested group is embedded as a nested group**,
   `inputs/composition/species`, so the embedded copy has the same shape as the file.
8. **`[diagnostics]` is not acted on.** Step 5 builds the figures. The Lindal manifest carries no
   `[diagnostics]` section, so nothing is skipped silently for this profile.

## 3. Acceptance results

Run by `reports/step02_4/accept_step02_4.py`; full output in `reports/step02_4/output.txt`.

| Check | Measured |
|---|---|
| 1. `casspian-refrac` writes the product, which reads back as kind N | **Pass.** The entry point exits 0 and writes `occul_data/lindal/lindal_refractivity.nc`. `read` returns a DataTree; kind `refractivity`, role `reduction`. Ten groups as listed in section 1; seven `input_hashes` entries. |
| 2. `read` refuses a copy missing `inputs/wind` | **Pass.** `kind refractivity requires the group 'inputs/wind', which is not in the file`. |
| 3. `read` refuses a copy with the thermo `input_hashes` entry edited | **Pass.** The first hex digit was changed: `input_hashes does not list the thermo hash the reduction recorded (ff9c12991399f94c...)`. |
| 4. Every embedded input group identical to its file (`xarray.Dataset.identical`) | **Pass**, all seven: thermo, composition, composition/species, geodesy, gravity, rotation, wind. |
| 5. The six scalars match Step 2 and Step 3 | **Pass, bit for bit** against a fresh reduction. `phi_c` 30.804949289059618 deg (+- 0.109208, 1sigma). `psi` 5.495050679420528 (+- 0.00661488). Label 36.3 (0.2, range, copied). Anchor pressure 10000 Pa (NaN). `r0` 58,519,883.28568849 m (+- 16,376.8, 1sigma). `h_ref` 90,000 m (NaN). Per level `radius_m`, `number_density_m3` and `refractivity` are identical too. |
| 6. The SPEC_00 section 2.2 listing is complete | **Pass**, all thirteen entries present. |
| 7. Beyond the specification: companions and `reduction_record` | **Pass.** Twelve companions with their kinds, terms and conversions (output lists each). The record has every attribute Step 4 names plus `fixed_point_closure_deg` (-3.15e-8 deg) per REVIEW_02_step2. |
| 8. Beyond the specification: no false warning | **Pass.** Reading the product with every input unchanged raises no `input_hashes` warning. |

**The anchor rule spread at the anchor.** At `phi_c` = 30.8049 deg, `north_pole` anchoring gives
`r0` = 58,537.389 km (+17.506) and `south_pole` gives 58,502.408 km (-17.475), against
58,519.883 km under `mean_polar_radius`. The SPEC_01 Step 7 figure of about +-19 km was at the
equator. At the anchor latitude it is +-17.5 km, larger than the propagated 16.4 km radius
uncertainty and recorded beside it, not in it.

## 4. Findings

1. **The anchor rule spread exceeds the propagated uncertainty at the anchor.** The spread is
   +-17.5 km; the declared uncertainties propagate to 16.4 km. Keeping them apart is the ruled
   design (SPEC_02 decision 5), and both numbers are in the file. It is worth knowing now that
   the choice of anchor rule is the largest single contribution to where the profile sits in
   radius. The Monte Carlo wrapper is where that choice gets exercised.

2. **The naming rule needed `_kg_mol`.** `UNIT_SUFFIXES` had no molar mass suffix. The only
   variable it affects today is the new companion; kind C's `molar_mass_kg_mol` has no
   companion. SPEC_00 section 5 names the rule but not the suffix list, so no specification text
   changes. Noted because it touches `lib.control`'s physical value detection too, which now
   also recognizes `_kg_mol` keys.

3. **SPEC_00 section 8 lists the `input_hashes` warning for every reader.** It is implemented for
   kind N only, the one kind that records its inputs' hashes twice and embeds them. For the
   other derived kinds, whose inputs (the raw bundle, the master table) live outside their
   directory, the same warning could be added to `read` generally. Not done here, since no step
   asks for it.

## 5. Regression

The ten SPEC_01 suites and the SPEC_02 Step 1 to 3 suites were rerun after the library changes of this step (`reports/step02_4/regression.txt`). All pass: 6, 14, 8, 5, 7, 7, 6, 8, 9, 6 for SPEC_01, and 6, 7, 9 for SPEC_02 Steps 1 to 3.

The first run showed SPEC_01 Step 9 at 5 of 6. Its check 6 required `lindal_refractivity.nc` to be absent from the profile directory, which was true until this step wrote it. The check now expects the product once it exists, the last entry of the SPEC_00 section 2.2 listing, and still flags any other unexpected file. Rerun, it passes 6 of 6.

No product under `occul_data/` other than `lindal_refractivity.nc` was written. No em dash or en dash appears in any file written in this step.

## 6. Next step

Step 5, the standard diagnostics, is not started. It waits for the review of this report.

---

## 7. Addendum after review: the anchor rule spread by a full rerun (SPEC_02 v0.7)

Added 12 September 2026, per `REVIEW_02_step4.md` ruling 1.

**The change.** The anchor rule spread is now the whole Step 2 fixed point, rerun under
`north_pole` and under `south_pole`. It is no longer two marches at the frozen latitude.
`_anchor_rule_spread` calls `freeze_anchor` on the same inputs with a copy of the manifest whose
`anchor_rule` is replaced, so each rule's result is built exactly as the product's own pair was.
`reduction_record` now carries `r0` and `phi_c` under each rule:
`anchor_rule_spread_r0_north_pole_m`, `anchor_rule_spread_r0_south_pole_m`,
`anchor_rule_spread_latitude_planetocentric_north_pole_deg`, `..._south_pole_deg`. The note says
the fixed point was rerun and that the spread enters no companion. The two latitude held
marches are gone from the code.

**Acceptance, rerun.** Nine of nine: checks 1 to 8 as before, plus check 9 for the v0.7 expected
values.

| rule | `r0` | against `mean_polar_radius` | `phi_c` | shift |
|---|---|---|---|---|
| `north_pole` | **58,537.683 km** (58,537.7 +- 0.1) | +17.800 km | **30.801965 deg** (about 30.8020) | -0.002985 deg |
| `mean_polar_radius` | 58,519.883 km | | 30.804949 deg | |
| `south_pole` | **58,502.116 km** (58,502.1 +- 0.1) | -17.767 km | **30.807925 deg** (about 30.8079) | +0.002976 deg |

The half range is 17.784 km, against the latitude held +17.506 and -17.475 km of the first
build. The difference is the coupling through `dphi_c/dr_anchor`, as the review predicted. The
spread is still larger than the 16.4 km the declared uncertainties propagate to.

**Regression.** Section 5 stands. The change touches only `refrac/product.py`, which no earlier
suite imports.

**The sweep.** After the acceptance commit, `lindal_refractivity.nc` was rebuilt on a clean tree
by `casspian-refrac occul_data/lindal/lindal_reduction.toml`. This is the first recorded hash of
the file.

| File | Kind | SHA-256 |
|---|---|---|
| `lindal_refractivity.nc` | refractivity | `327f814103bbfdefe86ca4bee262c6ff272e5f4cdc556a31a2e83fb1a725f584` |
