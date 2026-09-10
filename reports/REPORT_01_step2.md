# REPORT 01, Step 2. Stage one of the Lindal tool, the raw bundle

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 10 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.5, Step 2, against
`SPEC_00_Architecture_and_Data_Files.md` v0.6. Status: reported, awaiting review.

**The working tree is not committed**, per the commit on acceptance rule. The bundle on disk
therefore carries `-dirty` in `casspian_git_commit` and is to be rebuilt after the acceptance
commit.

---

## 1. What was built

**Housekeeping commit `3625059`, before any Step 2 code**, as SPEC_01 section 0 requires. The
four report and review files moved from `docs/specs/` to `reports/` with `git mv`, so the
history follows them, and the ignore rule narrowed from `reports/` to `reports/step*/` so that
the build record is committed while each step's working directory stays out. `docs/specs/` now
holds only the two specifications, `STATE.md` and the figure, matching SPEC_00 v0.6 section 2.

Two commits preceded it: `92482c0` carrying the author's v0.6 and v0.5 documents and the review
addendum, and `2573207` correcting Step 1 to the SPEC_00 v0.5 uncertainty naming rule. On that
correction: `lib.schema.uncertainty_companion` now derives every companion name by inserting
`_uncertainty` before the unit suffix, both per-variable overrides were removed, and the
registry now carries none. The names it produces match the section 6 tables exactly and are
printed by check 7 of `reports/step1/verify_review_changes.py`, which runs fourteen checks, all
passing.

**Deliverable: `src/casspian/tools/lindal/build_raw.py`**, with the console entry point
`casspian-lindal-raw` declared in `pyproject.toml` (Step 0 item 4 said entry points are added as
tools appear; this is the first). It takes `--raw-dir`, defaulting to `occul_data/lindal/raw`,
and `--output`, defaulting to `lindal_raw.nc` in that directory.

It writes twelve groups. `table1` holds `pressure_Pa`, `temperature_K`, `nh3_mole_fraction` and
`height_m` on a `level` dimension ordered by increasing pressure, each with `units`,
`long_name`, `provenance` and `value_source = "Table I"`, with `provenance = "measured"` on
`height_m` and `"derived"` on the other three. `scalars` holds one sub-group per TOML table:
`source`, `latitude`, `composition`, `geodesy`, `gravity_check_values`, `rotation`,
`wind_used_by_source`, `top_boundary`, and beneath `geodesy` the two nested tables
`surface_100mbar` and `surface_1bar`.

## 2. Decisions and observations

1. **Unit conversion is done in exact decimal arithmetic on the text as written in the CSV,
   then rounded once to float64.** This is not cosmetic. Multiplying the parsed double instead
   disagrees in **15 of the 207 transcribed cells**: 2.51 mbar becomes 250.99999999999997 Pa
   rather than 251.0, 131.83 mbar becomes 13183.000000000002 rather than 13183.0, and 5.1 ppm
   becomes 5.0999999999999995e-06 rather than 5.1e-06. **None of the fifteen is a cell the Step 2
   acceptance names**, so the acceptance checks alone would not have caught the difference. The
   rule the tool follows is that the bundle holds the double nearest the decimal value the
   source printed, and nothing else.

2. **Nested TOML tables become nested groups.** `[geodesy.surface_100mbar]` is a table, and the
   specification says every table of the TOML is reproduced as a sub-group of the same name, so
   it becomes `scalars/geodesy/surface_100mbar` while the non-table keys of `[geodesy]` stay as
   attributes on `scalars/geodesy`. Flagged as a reading rather than a certainty.

3. **`input_hashes` and `raw_sources` are not the same list.** `input_hashes` carries the two
   files the tool actually read, per SPEC_00 section 5. `raw_sources` carries those two plus
   `notes.md`, per Step 2, which asks for the notes to be hashed so the record is complete even
   though the tool does not read them. Both are in the root globals below.

4. **Attribute order inside a group is not the order of the TOML.** netCDF does not preserve
   insertion order through a write and read cycle, so a second reader checking section 4 below
   against the paper will find the keys of each group shuffled relative to
   `lindal_scalars.toml`. No value is affected. Said here so it is not mistaken for a defect.

5. **`table1` carries three extra group attributes** that the specification does not ask for:
   `source_columns`, `level_order`, and `rows_reordered_from_file`. The last records whether the
   tool had to sort the rows to satisfy the increasing pressure requirement. It reads `no`: the
   CSV was already in that order, and the check is in the file rather than only in this report.

6. **`pressure_Pa` in `table1` carries `positive = "down"` and `direction = "increasing"`**
   although kind `raw` is exempt from the coordinate rules of SPEC_00 section 5. They cost
   nothing and mean exactly what section 5 says they mean.

7. **The bundle carries `-dirty`**, as the rule expects: `casspian_git_commit` reads
   `36250593cd522c8a78be0aa41d0ad5086fadfeb1-dirty`, the housekeeping commit plus the
   uncommitted Step 2 code. On acceptance the tool is rerun so the committed bundle carries a
   clean commit.

## 3. Acceptance results

Run by `reports/step2/accept_step2.py`; full output in `reports/step2/output.txt`.

| Check | Result |
|---|---|
| 66 levels, ordered by increasing pressure | **Pass.** `level` size 66, strictly increasing in pressure. The CSV was already in that order, so no reordering was needed. |
| The 1 bar row reads exactly 100000.0 Pa, 134.8 K, 0.0 m | **Pass**, by exact equality, at index 59. |
| The top row reads 20.0 Pa, 138.7 K, 376700.0 m | **Pass**, by exact equality. |
| The bottom row reads 129848.0 Pa, 146.2 K, -14100.0 m | **Pass**, by exact equality. |
| Exactly nine finite NH3 values, from 831.76 mbar (2.6 ppm) to 1258.93 mbar (66.9 ppm), NaN at 1047.13 and 1298.48 mbar | **Pass.** Nine finite values; first at 83176.0 Pa = 2.6e-06; last at 125893.0 Pa = 6.69e-05; both stated gaps are NaN. No filled value appears in the bundle. |
| `scalars/latitude` carries `planetographic_deg = 36.3` and the swath | **Pass.** `planetographic_deg = 36.3`, `swath_planetographic_deg = [36.3, 36.7]`, `value_source = "Fig. 4 label on the Voyager 2 ingress profile: 36.3 N"`. |
| `read(path, "raw")` succeeds | **Pass.** `casspian_kind = "raw"`, schema version 1, `created_by = "casspian-lindal-raw 0.1.0"`. |
| `read(path, "thermo")` refuses | **Pass.** `lindal_raw.nc: declared kind is 'raw' but 'thermo' was asked for (SPEC_00 section 8).` |

Eight of eight pass. No em dash or en dash appears in any file written in this step.

## 4. The scalars transcription, for a second reader

Every attribute of every `scalars` group as it reads back from the bundle, so that the
transcription can be checked against the paper without opening the netCDF. The root globals
follow. Ordering is as netCDF returns it, per decision 4 above.

```
[scalars]
  content = 'every scalar the source states, transcribed; one sub-group per TOML table'
  source_file = 'lindal_scalars.toml'

[scalars/source]
  table1_footnote = "This table gives the temperature (T), and the ammonia mixing ratio in parts per million as a function of the pressure level (p) in Saturn's atmosphere. Also provided is the local altitude (h) of the measurements relative to the 1-bar pressure level. The gas is assumed to consist of 94% H2 with the remainder being He and NH3. For other gas mixtures, one may scale the numerical data given here. The measurements were made during the ingress of Voyager 2."
  pressure_range_note = 'Profiles extend from about 0.2 mbar to 1.3 bar (p. 1138); X band extinguished below 1 bar, S band below 1.3 bar (p. 1138)'
  frequency_bands = 'S (2.3 GHz) and X (8.4 GHz)'
  frequency_bands_source = "Fig. 3 caption: 'Doppler frequency measurements at S and X bands'; Fig. 3 caption: upper NH3 segment from X band, lower from S band"
  observation_date = '1981-08-26'
  title = 'The atmosphere of Saturn: an analysis of the Voyager radio occultation measurements'
  observation_date_source = "Fig. 2 caption: data 'were acquired at DSS 43 on 26 August 1981'"
  spacecraft = 'Voyager 2'
  citation = 'Lindal, G. F., Sweetnam, D. N., Eshleman, V. R. 1985, AJ 90, 1136'
  event = 'ingress'

[scalars/latitude]
  value_source = 'Fig. 4 label on the Voyager 2 ingress profile: 36.3 N'
  uncertainty_note = 'Half the swath; the label is the defensible single value, not the true one (handoff 9A.8).'
  planetographic_deg = np.float64(36.3)
  swath_source = "p. 1138: 'The ingress gas data were acquired near the evening terminator in the region from 36.3 N planetographic latitude and 186.8 E system III longitude to 36.7 N latitude and 185.7 E longitude.'"
  manuscript_and_handoff_history = 'Carried as 36.5 through the project record until the Fig. 4 label was read on 9 Sept 2026 (handoff 9A.8).'
  uncertainty_deg = np.float64(0.2)
  uncertainty_kind = 'range'
  convention_source = "Fig. 4 caption: 'Each profile has been labeled with the planetographic latitude at which the lowest point on the radio ray touched the 100-mbar isobaric surface. (The lowest point on the ray is defined as the position where the ray is perpendicular to the local vertical.)'"
  longitude_system_iii_deg = [186.8, 185.7]
  swath_planetographic_deg = [36.3, 36.7]
  longitude_source = 'p. 1138, same sentence as the swath'

[scalars/composition]
  h2_fraction = np.float64(0.94)
  h2_uncertainty = np.float64(0.03)
  h2_source = "p. 1137: 'For the purpose of interpreting the refractivity data, we have adopted a hydrogen abundance of 94 +/- 3%; the remainder is assumed to be primarily helium.'"
  remainder = "He and NH3 (Table I footnote); 'primarily helium' (p. 1137); mean molecular mass 2.135 amu with 0.94/0.06 (Phase 1 check: 2.1351)"
  refractivity_conversion = "H2 and He only. Fig. 3 caption: 'In computing the total gas density from the refractivity, it was assumed that the atmosphere consists of 94% hydrogen with the remainder being mostly helium. For other nonpolar gas mixtures, one may scale the results presented here by noting that the number density is inversely proportional to the mean molecular refractivity.'"
  refractivity_sources = 'Essen 1953; Orcutt and Cole 1967; Bose et al. 1972 (p. 1138); no numerical values printed. Values used: see data_static/species_master.toml, set lindal1985.'
  alternative_noted = 'Conrath et al. 1984: 96.3% H2, 3.3% He, 0.4% CH4, m = 2.136 amu (p. 1137), mentioned and not adopted.'
  ortho_para_note = 'p. 1138: normal hydrogen (75% ortho) used in the refractivity measurements; ortho-para conversion does not produce a significant change in the refractive property of H2.'

[scalars/geodesy]
  fit_method = 'Least-squares fit of a computed geoid to the five observed 100 mbar radii (four Voyager occultations, one Pioneer 11); p. 1141 and Fig. 9 caption'
  fit_residual_km = [3.0, 4.0]
  fit_residual_source = "p. 1141: 'the standard deviations of these fits range from 3 to 4 km'; 'The isobaric and geodetic surfaces agree within about 4 km in the troposphere.'"
  fit_inputs = 'External gravity field (Null et al. 1981; Campbell 1984 pers. comm.), zonal winds at the cloud tops (Smith et al. 1982; Ingersoll and Pollard 1982; Garneau 1984 is also cited in the Fig. 9 caption), assumed independent of altitude (p. 1141, Appendix p. 1144), pole vector (Kozai 1957; Null et al. 1981; Simpson et al. 1983), System III rotation'
  fit_latitude_convention = "planetocentric (Fig. 9 abscissa; Appendix: 'r and phi denote the distance to the center of mass and the planetocentric latitude')"
  reference_geoid_polar_radius_km = np.float64(54438.0)
  reference_geoid_source = "Appendix Eq. 12: 'U_ref = U(r_p, pi/2), where r_p was set equal to 54 438 km'; Fig. 9 caption"

[scalars/geodesy/surface_100mbar]
  value_source = "p. 1141 (also abstract and Fig. 9 caption): 'the geodetic surface which in a least-squares sense best fits the 100-mbar data has a mean polar radius of 54 438 +/- 10 km, and an equatorial radius of 60 367 +/- 4 km. The oblateness of this surface (Req - Rp)/Req is 0.098 22 +/- 0.000 18.'"
  pressure_Pa = np.float64(10000.0)
  oblateness = np.float64(0.09822)
  radius_equatorial_km = np.float64(60367.0)
  radius_polar_mean_unc_km = np.float64(10.0)
  polar_asymmetry_note = 'Fig. 9 caption: south-polar radius may be of order 10 km greater than the north-polar radius, possibly seasonal.'
  radius_equatorial_unc_km = np.float64(4.0)
  oblateness_unc = np.float64(0.00018)
  radius_polar_mean_km = np.float64(54438.0)

[scalars/geodesy/surface_1bar]
  pressure_Pa = np.float64(100000.0)
  radius_equatorial_km = np.float64(60268.0)
  radius_equatorial_unc_km = np.float64(4.0)
  radius_polar_mean_km = np.float64(54364.0)
  radius_polar_mean_unc_km = np.float64(10.0)
  oblateness = np.float64(0.09796)
  oblateness_unc = np.float64(0.00018)
  value_source = "p. 1141 ('The corresponding values at the 1-bar level are 54 364 +/- 10 km, ...') and Table II p. 1143"

[scalars/gravity_check_values]
  g_equatorial_1bar_ms2 = np.float64(8.96)
  g_equatorial_unc_ms2 = np.float64(0.01)
  g_polar_1bar_ms2 = np.float64(12.14)
  g_polar_unc_ms2 = np.float64(0.01)
  value_source = "Table II, 'Acceleration of gravity' rows, Saturn column"
  atmospheric_temperature_1bar_K = np.float64(134.0)
  atmospheric_temperature_unc_K = np.float64(4.0)
  temperature_source = 'Table II, Saturn column; footnote: 94% hydrogen with the remainder being helium and ammonia'

[scalars/rotation]
  system = 'System III'
  period_hms = '10h 39m 22.4s'
  period_s = np.float64(38362.4)
  value_source = "Table II ('System III rotation period (Davies et al. 1983)') and Fig. 9 caption ('10h 39m 22.4s magnetic field or system III rotation period (Desch and Kaiser 1981; Davies et al. 1983)')"

[scalars/wind_used_by_source]
  citations = 'Smith et al. 1982; Ingersoll and Pollard 1982 (Appendix p. 1144); Garneau 1984 added in the Fig. 9 caption'
  latitude_convention = "planetographic (Appendix p. 1144: 'Published zonal-wind data ... give V_w as a function of the planetographic latitude')"
  vertical_structure = "assumed independent of altitude (Appendix p. 1144: 'the altitude variation was assumed to be negligible'; p. 1141)"
  frame = 'System III (Appendix Eq. 4: omega = omega_III + V_w / (r cos phi))'
  how_combined = 'not stated; project reading: Smith et al. 1982 Fig. 4 as the baseline table, others as cross-checks'

[scalars/top_boundary]
  statement = "not stated numerically. p. 1138: pressure profiles obtained 'by utilizing the equation for hydrostatic equilibrium and integrating the density from the top of the detectable atmosphere and downward'; the top of Table I is 0.20 mbar at 376.7 km, 138.7 K."
  uv_constraint_noted = 'p. 1141: ultraviolet data indicated 120 +/- 30 K at 0.05 microbar (Smith et al. 1983); above 0.2 mbar the temperature may decrease again.'

==============================================================================
Root globals
==============================================================================
  title = 'Lindal et al. (1985) Voyager 2 ingress, raw transcription bundle'
  profile_or_run = 'lindal'
  role = 'reduction'
  source = 'Lindal, G. F., Sweetnam, D. N., Eshleman, V. R. 1985, AJ 90, 1136'
  input_hashes = 'occul_data/lindal/raw/lindal_table1.csv sha256:1fc3a2c8b92805ddb00e295ce9aef19a7b851ac769a3ab8933c8193b5c2e30a1\noccul_data/lindal/raw/lindal_scalars.toml sha256:c879621c7663ea1cd2c9c3cc3202ddf2e2535bfe59d1bc0f1638d2d4f778c7ac'
  raw_sources = 'occul_data/lindal/raw/lindal_table1.csv sha256:1fc3a2c8b92805ddb00e295ce9aef19a7b851ac769a3ab8933c8193b5c2e30a1\noccul_data/lindal/raw/lindal_scalars.toml sha256:c879621c7663ea1cd2c9c3cc3202ddf2e2535bfe59d1bc0f1638d2d4f778c7ac\noccul_data/lindal/raw/notes.md sha256:34e353b9c52fe38f2711af15941b04a024d007c9ee74435660a4cd179941fcda'
  history = '2026-09-10T18:51:18Z casspian-lindal-raw: pressure converted from mbar to Pa by a factor of 100\n2026-09-10T18:51:18Z casspian-lindal-raw: ammonia converted from ppm to mole fraction by a factor of 1e-6\n2026-09-10T18:51:18Z casspian-lindal-raw: altitude converted from km to m by a factor of 1000\n2026-09-10T18:51:18Z casspian-lindal-raw: notes.md hashed into raw_sources as documentation, not read\n2026-09-10T18:51:18Z written by casspian-lindal-raw 0.1.0 as kind raw'
  Conventions = 'CF-1.10, CASSPIAN-0.1'
  casspian_kind = 'raw'
  casspian_schema_version = np.int64(1)
  created_by = 'casspian-lindal-raw 0.1.0'
  created_at = '2026-09-10T18:51:18Z'
  casspian_git_commit = '36250593cd522c8a78be0aa41d0ad5086fadfeb1-dirty'
  codata_release = '2018'
```

## 5. Questions for review

1. **Nested groups for nested TOML tables** (decision 2). If the intent was a flat
   `scalars/surface_100mbar` instead, that is a one line change, but nesting keeps the bundle
   the same shape as the transcription it came from.

2. **The scalars still want a second reader.** `notes.md` says the values were read from the PDF
   on 9 and 10 September and are marked for a second reader before the raw bundle is frozen.
   Section 4 above is that list in one place. Until it is checked, the bundle is correct as a
   transcription of the TOML but not yet confirmed as a transcription of the paper.

## 6. Next step

Step 3, the gravity and rotation tools and kinds G and R, is not started and is not to be begun
until this report has been reviewed. It is the first step to write a kind G file, so it is the
step that fixes `harmonic_convention = "CASSPIAN-J1"` into data.
