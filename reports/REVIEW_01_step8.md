# REVIEW 01, Step 8. The composition tool for Lindal, kind C

Review of `reports/REPORT_01_step8.md` against SPEC_01 v0.16 Step 8. Reviewer: S. Rafkin,
11 September 2026. **Disposition: accepted with changes.** SPEC_01 is revised to v0.17 and
`data_static/species_master.toml` is revised by the author. Make the changes, rerun, commit,
rebuild `lindal_composition.nc` clean, proceed to Step 9.

## The ammonia question

The author asked whether NH3 at zero aloft is Lindal's statement or our assumption. It is our
assumption, and the v0.5 wording of the rule made it look like a computation. Lindal's Table I
lists NH3 at nine levels from 831.76 to 1258.93 mbar, inferred from the S and X band
absorption where absorption was measurable (Fig. 3 caption: "ammonia vapor density profiles
inferred from the microwave absorption measurements"); above 831.76 mbar the column is blank,
which means not measured, not zero. What justifies zero there is physics, not the table:
Lindal states the troposphere at these levels is saturated with ammonia, and the saturation
mixing ratio over NH3 ice on his own temperatures runs from a few ppm near 0.8 bar (the
tabulated 2.6 ppm at 831.76 mbar is consistent with it) to below 1e-7 by 0.5 bar and below
1e-12 at the tropopause. The assumption costs under 1e-6 of the mean molar mass and nothing
in the refractivity, since the `lindal1985` set carries NH3 at zero on the Fig. 3 caption
reading (the conversion used H2 and He only; the Table I footnote puts NH3 in the remainder
for the temperature, not for the refractivity). v0.17 states the rule in three parts:
interior gap interpolated, below the deepest value extrapolated, above the highest value zero
by assumption with the reason written into an attribute, provenance `assumed`. No upward
extrapolation and no clamp; that construction is what put 0.2 ppm at 794.33 mbar and created
finding 1. The acceptance number 794.33 was never a computed crossing; it is the first grid
level above the highest tabulated value, and with the restated rule it holds exactly.

## Rulings on the findings and decisions

1. **Finding 1.** Resolved as above. The interior gap value (15.931 ppm) and the deep
   extrapolation (79.253 ppm) stand.

2. **Finding 2, `is_polar`.** Agreed, and a good catch: writing ammonia as nonpolar would have
   been wrong on a molecular fact. The property now lives in an `[is_polar]` table beside
   `molar_mass` in `species_master.toml` (author's edit, committed with this review) and the
   per-set copies are removed; the tool reads it from there. `temperature_dependence` stays
   per set. The same edit promotes the `lindal1985` NH3 entry from `unverified` to
   `verified_reading` with the caption and footnote quoted, since the report's reading is the
   right one and had been left unrecorded in the table.

3. **Finding 3.** Correct and worth keeping in the report: ammonia enters this reduction only
   through the mean molar mass, at 0.055 percent at the deepest level.

4. **Decisions 1 to 5.** All accepted: the declared split checked against the transcribed
   `h2_fraction`; `x_H2_uncertainty` from the raw bundle (a stated value is data); the
   per-level `nh3_provenance` flag; the schema-correct `latitude_planetocentric_absent_meaning`
   (v0.17 spells it out, and Step 9 is to use the same form); repository-relative paths in the
   product. The last one is the kind of thing that would have been found much later at much
   more cost.

## Actions

- Implement the v0.17 NH3 rule (zero above 831.76 mbar with `nh3_provenance = assumed` and the
  `nh3_aloft_assumption` attribute; interior interpolation and downward extrapolation as now).
- Read `is_polar` from the new `[is_polar]` table.
- Rerun the acceptance (seven of seven expected), commit, rebuild the product clean, Step 9.
