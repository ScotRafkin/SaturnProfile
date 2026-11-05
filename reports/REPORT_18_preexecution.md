# REPORT 18, pre-execution: Kliore et al. (1980) as a second registered profile

SPEC_18 v0.2, Step 0, before any code changes. 5 November 2026. For the author and the reviewing
agent. `main` is at `6802bed`, pushed. Nothing in the tree has changed except this report.

Pioneer 11's exit is 11.9 S planetographic (9.5 S planetocentric).

**In short.** The model takes Kliore unchanged. Two of Lindal's tools need splitting
(`build_composition`, `build_inputs`), and neither moves a Lindal value in the proposal. I tested
this with a throwaway dry run outside the repository (§4), using scratch scripts that are not
deliverables. It found four things the specification does not settle, which §5 asks to have ruled.
The one that matters most: **Kliore's closure does not meet Lindal's bounds**, by up to a factor of
ten, and part of that is not print rounding.

## 1. The line, tool by tool

| Tool | Layer | Change | Why |
|---|---|---|---|
| `tools/lindal/build_raw.py` | source-specific | none | Table I's columns and its pressure grid are Lindal's. `read_scalars`, `_attribute_value` and `_scaled` are general and are imported as they stand. |
| `tools/gravity/build_gravity.py` | general | none | Takes source and epoch from the build file. |
| `tools/gravity/build_rotation.py` | general | none | The same. |
| `tools/wind/build_wind.py` | general | none | Reads only `scalars/geodesy` `reference_geoid_polar_radius_km` from a raw bundle. Kliore's bundle carries that group. |
| `tools/wind/shear.py`, `tools/run/run_inputs.py`, `tools/run/new_run.py` | general | none | The closure uses `case = "identity"`. |
| `tools/composition/build_composition.py` | mixed | **split** | The writer is general. The ammonia fill and its text are Lindal's. |
| `tools/lindal/build_inputs.py` | mixed | **split** | Kind D's writer, the companion check and the manifest are general. The reading of Table I and the pressure grid are Lindal's. |
| `refrac`, `forward`, `lib`, `tools/plots` | model | none | §4 shows they run on Kliore as they are. |

**`build_composition`.** One tool, with this separation:
- **General:** the levels from `table1/pressure_Pa`, the remainder split and its check against the
  raw bundle's `h2_fraction`, the species group, the latitude grid and the writing.
- **Lindal's:** `fill_ammonia`, `x_NH3`, `nh3_provenance`, `nh3_aloft_assumption`,
  `nh3_tabulated_range_Pa`, the species group's `nh3_note`, and the "94 +/- 3 percent" wording.
- **The switch** is whether the raw table has an ammonia column. `nh3_fill_rule` becomes optional.
  No new key.
- `latitude_planetographic_deg` is written only if the raw bundle states it. Only the kind T figure
  reads that attribute under that name, and it reads kind T's, so nothing reads kind C's.

**`build_inputs`.** `_write_unless_unchanged`, kind D's writer, the companion check and the manifest
writer move to a new module, `tools/profile_inputs.py`, and `build_inputs` calls them. Kind T's
writer stays in each source's tool, because what it reads is the source's own.

**The name `casspian-composition-lindal` stays.** The docstring will say it is general. Renaming
would touch the runbook, the guides and an accepted suite, and would move no number.

## 2. Kliore's preprocessing

**Where:** `src/casspian/tools/kliore/`, as `build_raw.py` (`casspian-kliore-raw`) and
`build_inputs.py` (`casspian-kliore-inputs`), parallel to Lindal's.

**What goes where:**

| File | Holds |
|---|---|
| `occul_data/kliore1980/pioneer11_exit.csv` | The transcription, 128 levels, untouched. |
| `raw/kliore1980_scalars.toml` | §1's scalars, with sources in comments: citation, date, season, mixture, f = 0.096, 9.5 S; and Lindal's two surfaces and reference polar radius (ruling 4). |
| `raw/kliore1980_raw.nc` | The table on 129 levels, including `radius_local_m` and an `inserted_level` flag; the scalars as groups. |
| `kliore1980_build.toml` | The sections Lindal's has, plus the comments §1 asks for. |

**Stage one (raw bundle).** No other file needed.
- Undoes the equatorial mapping with f = 0.096 on the planetocentric reading (SPEC_16).
- Converts units in exact decimal, as Lindal's does.
- Inserts the 100 mbar row, linear in ln p between 99.6 and 101.9 mbar. In the dry run: 89.75 K,
  fraction 0.1756.

The inserted row goes in the bundle, not only in kind T, because kind C reads its levels from the
bundle and must match kind T's bit for bit.

**Build order:** raw, G, R, W, C, then stage two. Stage two needs G, R and W because of ruling 1.

**Stage two (kinds T and D, the manifest).**
- **Kind T:** pressure, temperature, height (ruling 2), the label latitude (ruling 1), and
  `pressure_printed_Pa`. The closure's F6 figure reads that last one, and the dry run failed without
  it. On the inserted level it holds 100.0 mbar, and the flag marks the row.
- **Kind D:** Lindal's surfaces, through the shared writer.
- **Manifest:** anchor isobar 1.0e4 Pa, rule `equatorial_radius`, through the shared writer.

**Season.** Lindal's date has a computed season (JPL approximate elements, IAU pole). The same method
reproduces Lindal's values within rounding (18.17 and 8.06 against 18.2 and 8.06), and for
1 September 1979 it gives **Ls 353.7, subsolar latitude -2.85**. The wind file keeps its own date,
1981-08-26, and the build file's comment says it is a Voyager wind used two years later.

**Closure run.** `forward/kliore1980_closure/`, made by `casspian-new-run` from `lindal_closure`, with
the paths, date, season and titles edited. `[composition]` points at Kliore's raw bundle with
species H2 and He.

## 3. The Lindal rebuild

**Scope:** the two changed tools only. Kind C, kind T, kind D and the manifest of `occul_data/lindal/`,
and each of the 23 runs' kind C files, since `run_inputs` calls the same tool.

**How:** in a copy, with the old files **deleted first**. `_write_unless_unchanged` otherwise keeps
what is there, and the comparison would prove nothing.

**The comparison:** every variable array-equal, NaN included. Every attribute the rebuild writes is
equal, except those the writer fills and the lineage ones (`history`, `control_file`,
`input_hashes`, `raw_bundle`). Global attributes the tools stopped writing at SPEC_R1 cannot be
compared, so they are left out. The manifest is compared byte for byte.

**Baseline, with today's tools, for `occul_data/lindal/` only:**

| File | Variables | Differ in value | Attribute differences |
|---|---|---|---|
| `lindal_composition.nc` | 15 | 0 | 3, all globals SPEC_R1 retired |
| `lindal_thermo.nc` | 7 | 0 | 14, all globals SPEC_R1 retired |
| `lindal_geodesy.nc` | 8 | 0 | 2, all globals SPEC_R1 retired |
| `lindal_reduction.toml` | | byte-equal | |

So the test can pass today, and a failure after the split will mean the split. The 23 run files were
not part of this baseline; Step 1 will run them.

## 4. The dry run

Kinds T, C and D were written by hand in scratch, G, R and W were copied from Lindal, and `refrac` and
`forward` ran unmodified with figures on.

**The reduction** runs, and F1 to F4 draw.

| Measure | Result |
|---|---|
| `phi_c` returned from the label | -9.50000008 deg (tolerance 1e-6) |
| Refractivity against Table 1, 128 levels | RMS 0.031 n-units, largest 0.086, mean +0.006 |
| Inserted 100 mbar level | N = 37.82 n-units, between the table's 37.6 and 38.4 |
| Kind N radius minus Kliore's local radius | +1.27 km, no spread |
| `refractivity_uncertainty` finite | 0 of 129 (ruling 3) |

- The RMS is lower than SPEC_16's 0.034 because that figure predates the 17.5 mbar correction.
- **The +1.27 km is not comparable with SPEC_16's +9.2 km.** That offset is the delivered radius under
  the Sanchez-Lavega and IRIS wind. This one is the reduction's own geoid under Lindal's
  altitude-independent wind. The 8 km between them is about what SPEC_17 found at the southern points
  (+7 to +8 km against Lindal's surface). I have not tested that as the explanation.

**The closure** runs, and F5 and F6 draw. T and p residuals agree to 4e-16.

| Bin (tabulated pressure) | Kliore, max abs residual | Lindal's measured | Lindal's bound |
|---|---|---|---|
| above 2 mbar | 0 (one level, the boundary) | 9.8e-4 | 1.5e-3 |
| 2 to 10 mbar | 3.26e-2 | 2.2e-3 | 3e-3, from 2 mbar down |
| 10 to 100 mbar | 5.69e-3 | 2.6e-3 | 3e-3, from 2 mbar down |
| below 100 mbar, without the bottom row | 9.07e-3 | 1.7e-3 | 3e-3, from 2 mbar down |
| bottom row, 182.2 mbar | +7.40e-3 | -3.3e-3 | 4e-3 |
| mean below 10 mbar | +3.56e-3 | -4.9e-4 | 1.5e-3 |

In temperature, that is 3.9 K at most above 5 mbar, 0.5 K from 10 to 100 mbar and 0.85 K below.

**There are two causes, and I have only established one:**
- **Print rounding, at the top.** Table 1 prints pressure to 0.1 mbar, which is 2.8 percent at 1.8 mbar
  and 0.5 percent at 10 mbar. Lindal's table is on a declared grid and has no such term.
- **A drift, deeper down, that rounding does not explain.** The mean residual grows from near zero at
  10 to 30 mbar to +3.4e-3 at 30 to 100 mbar and +5.5e-3 below. The rounding estimate there is about
  1.9e-3. Using that estimate (printed pressure plus 0.1 km in radius over a scale height, with g
  taken as 9.5 m/s2), 41 of 129 levels fall inside it. The drift could come from a difference between
  Kliore's own hydrostatic integration and ours, such as his gravity, but I have not tested that.

**F6's envelope will be too narrow for Kliore.** It infers the height column's precision from its
digits, and Kliore's heights are derived values, so the height term comes out near zero. I propose
leaving the figure alone and computing the envelope in the acceptance script from the printed 0.1 km.

## 5. Rulings asked for

**1. The label latitude.** §1 item 3 settles 9.5 S planetocentric, but `refrac` reads only a
planetographic label and solves for `phi_c`.
- **Proposed:** stage two computes the label as `phi_c + psi(phi_c)` on the geoid `refrac` will march.
  That gives -11.8934 deg, and `refrac` returns 9.5 S. It also lands next to Lindal's 11.9 S, which
  is a useful independent check.
- **Rejected:** Lindal's 11.9 S returns 9.51 S rather than §1's value. Kliore's own ellipse gives
  11.57 S, which is the wrong geoid for `refrac`.

**2. Height.** `refrac` computes `r = r0 + (h - h_ref) cos psi`, treating `h` as a distance along the
local vertical, but Kliore's column is radial.
- **Proposed:** `height_m = (r_local - r_local at 100 mbar) / cos psi`, so kind N's radii keep Kliore's
  spacing exactly. The datum is his own 100 mbar level.
- **Otherwise:** using the radial difference directly puts the spacing off by 142 m at the top and
  20 m at the bottom.

**3. Uncertainties.** Kliore gives asymmetric extreme ranges on every fifth row, and §1 does not rule on
them.
- **Proposed:** NaN in kind T, as for Lindal, with the ranges kept in the raw bundle. `x_H2_uncertainty`
  is also NaN, because Kliore's +/-3 K for composition is not a mole fraction uncertainty.
- **Consequence:** kind N's `refractivity_uncertainty` is NaN at every level. That does not affect this
  specification, but the day Kliore becomes a second anchor, `sigma_ln_N_measurement` is NaN and the
  combination will have to supply it.
- **Alternative:** filling the rows between means inventing values, and using only every fifth row
  leaves a column that is mostly NaN. I don't recommend either.

**4. Kind D's source.** §3 item 2 puts §1's scalars in Kliore's raw bundle.
- **Proposed:** transcribe Lindal's two surfaces and reference polar radius into
  `kliore1980_scalars.toml`, with Lindal's page cited. Acceptance checks them equal to Lindal's kind D.
  The transcription repeats a few numbers, but the shared writer and `build_wind` then read one form.
- **Alternative:** read them from Lindal's raw bundle, which avoids the repetition but ties Kliore's
  build to Lindal's file.

**5. The closure's acceptance (the main one).** §4 item 4 says "the same measures SPEC_03 Step 4 used".
- **Proposed:** the same measures and bins, reported against the printed envelope, **without Lindal's
  bounds**, which come from Lindal's table. Passing means that the production completes and that the T
  and p residuals agree to 1e-12. The drift is reported as measured.
- **Your call:** whether the drift should be explained before Step 1 closes. It could be split with one
  extra closure in a copy, for example without the wind. I would do that only if you ask, since this
  anchor is not used until time propagation exists, but an unexplained 0.5 percent in a registered
  closure is a fair thing to want resolved first.

## 6. What Step 1 will not do

It will not rerun other suites (SPEC_18 §0), rewrite anything under `occul_data/lindal/` or
`forward/lindal_*`, change `refrac` or `forward`, add a key that no code reads, or use Kliore against
the Voyager-epoch wind.
