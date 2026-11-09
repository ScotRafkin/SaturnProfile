# REPORT 20, pre-execution: sections, and the Voyager profiles as anchors

SPEC_20 v0.2, Step 0, before any code changes. 9 November 2026. For the author and the reviewing
agent. `main` is at `e4292c7`, pushed. No code has changed. Beside the author's uncommitted SPEC_20,
STATE.md and status note, only this report and STATE.md's Step 0 row are new. The checks of §2 were
read-only, on SPEC_19's products.

Latitudes are planetographic with planetocentric in parentheses: Voyager 2 egress 31.2 S (26.1 S),
Voyager 1 ingress 74.2 S (71.1 S), Voyager 1 egress 3.0 S (2.4 S).

**In short.**
- **The anchors take SPEC_18's line with no change to the model.**
- **Voyager 1's radius needs no outside point.** `refrac` uses a profile's heights only relative to
  100 mbar. The absolute radius comes from kind D and the geoid march, as for Lindal and Kliore. So
  for Voyager 1 I propose heights by hydrostatic integration from 100 mbar, and Fig. 9 as a check, not
  an input (§1.2).
- **The section is a byproduct of one march per hemisphere,** if the production is redone at each
  latitude. I checked this read-only on SPEC_19's runs: the march's path agrees with separate runs to
  1e-3 K (§2).
- **Nothing blocks.** Six choices are listed in §3 for a ruling.

## 1. The anchors

### 1.1 The chain, by SPEC_18's three layers

Each new anchor gets its own files under `occul_data/lindal/`, with a prefix per profile:
`voyager2_egress`, `voyager1_ingress`, `voyager1_egress`. Lindal's files are read and never written.

| Step | Layer | Tool | For the new anchors |
|---|---|---|---|
| Raw bundle | source-specific | **new** `tools/lindal/build_profile_raw.py` | Reads the traced CSV (SPEC_15) and writes the bundle the general tools need: `table1/pressure_Pa`, temperature, `height_m`, and the groups `scalars/composition`, `latitude`, `source` (epoch, season, the reliability notes) and `geodesy` (Lindal's two surfaces, as Kliore's transcribes them) |
| G, R, W | general | none | Lindal's reduction set, used as is (the manifest names files in its own directory, so these must sit beside it, as SPEC_20 places them) |
| C | general | `build_composition`, unchanged | Each anchor gets **its own C on its own levels.** `refrac` and the closure require the levels bit for bit (`lib/control.py:410, 1081-1087`), so Lindal's C cannot serve |
| T | source-specific | **new**, in the same module | Writes kind T from the bundle, as `tools/kliore/build_inputs.py` does |
| D, manifest | general | `tools/profile_inputs.py`, unchanged | **Each anchor's D under its own name.** `write_geodesy` rewrites whenever content differs (`profile_inputs.py:84-104, 253`), and a different epoch would overwrite `lindal_geodesy.nc` |
| N | model | `casspian-refrac`, unchanged | |
| Closure | general | a copy of `forward/kliore1980_closure` | Per anchor, under `forward/<prefix>_closure/`, editing only what the Kliore copy edited |

One new module, two commands (`casspian-voyager-raw` and `casspian-voyager-inputs`). It serves all three
profiles, since their CSVs share a format. No shared tool changes, so acceptance 1 needs no rebuild of
Lindal's files in a copy.

**What the profiles need on ingest:**
- **Voyager 2 egress** carries altitude, 391 rows from 1337.7 to 0.22 mbar. Its rows are in decreasing
  pressure; the closure needs them top down (`production.py:399`), so they are reversed. **It has no
  exact 100 mbar row,** and `refrac` requires one (`control.py:442-458`), so one is inserted, linear in
  `ln p`, as Kliore's was (`tools/kliore/build_raw.py:70-87`). Its pressures are ours: SPEC_15 integrated
  them from 1 bar with Lindal's C, G, R, W and D.
- **The two Voyager 1 profiles** are on 100 levels per decade, so 100.0 mbar is a row already. 385 rows
  from 1000 to 0.14 mbar (74.2 S) and 345 rows from 1148 to 0.42 mbar (3.0 S). They carry no height (§1.2).

### 1.2 The radius of the Voyager 1 profiles

The spec suggests Lindal's Fig. 9 100 mbar point with hydrostatic integration from it. I checked how
the radius enters the reduction:
- `refrac` reads `height_m` and uses only `h - h_ref`, the height above the exact 100 mbar level
  (`refrac/reduce.py:256-258, 335`).
- The absolute radius at 100 mbar, `r0`, is kind D's surface carried by the geoid march to the anchor
  latitude (`refrac/anchor.py:137-148, 210-213`). It is not the profile's own.

So the Fig. 9 point, whatever it says, does not set the anchor's radius; kind D does, as for Lindal and
Kliore. **The radius a Voyager 1 profile needs is only the one gravity is evaluated at** while
integrating its heights.

**Proposal:**
- **Heights:** `h = 0` at 100 mbar, integrated hydrostatically up and down in `ln p` from the profile's
  own T, with `m_bar` from its own kind C and `g` from Lindal's G, R and W. That is the setup of SPEC_15's
  `egress_profile.to_pressure`, run the other way: pressure to height instead of height to pressure.
- **Gravity's radius:** `r0` from kind D at the anchor latitude, the same `r0` that `refrac` will use,
  so the heights and the reduction agree.
- **Fig. 9 as a check:** its altitude at each profile against the delivered 100 mbar surface. REPORT_17
  §4 already found ours above Lindal's by +7.1 km at 74.2 S and +3.3 km at 3.0 S.

**What it costs:**
- **Gravity:** taking the Fig. 9 point instead would move gravity's radius by about 3 to 8 km. That is
  2.4e-4 in `g` per 7 km, about 70 m over the 300 km column. These are estimates, not measured; I
  expect the effect on delivered temperature to be under 0.01 K.
- **Heights are hydrostatic by construction,** so the pressure identity is about integration error, not
  the tracing. The closure (acceptance 2) still compares delivered T with the traced T at each pressure,
  and that does test the tracing.
- **Integration within the tool versus within Lindal's figure:** the heights are ours, as the Voyager 2
  egress pressures already are.

## 2. The section

### 2.1 One march per hemisphere is enough

The forward transfer already traces every anchor isobar across its whole mesh, and stores the
isobars' geopotential and `ln N` at every mesh latitude (`forward/transfer.py:287, 361, 603-615`;
product group `/isobars`). Only the hydrostatic production, the column, the altitude and the datum are
done at the target alone (`produce_at_target`, `transfer.py:663-734`). The mesh covers from one spacing
beyond the anchor to one spacing beyond the target (`lib/mesh.py:138-177`). So **one march toward each
pole covers the globe,** and the section is read off the two runs.

**Checked on SPEC_19's products, read-only.** The 74.2 S run's path, read at 31.2 S and 3.0 S, against
the separate runs ending there:

| | 31.2 S | 3.0 S |
|---|---|---|
| Hydrostatic production redone on the path, minus that run's delivered T | 3.5e-4 K | 9.4e-4 K |
| The isobar label instead, `T = p_label r_bar / (k N)`, minus delivered T | 0.22 K max, RMS 0.08 K | 0.28 K max, RMS 0.08 K |

The label's gap is the pressure identity residual, up to 2.75e-3. **So the section tool redoes the
production at each section latitude** (`fp.produce_on_geopotential`, as `produce_at_target` does). That
is cheap and meets acceptance 3's 0.1 K with room to spare; the label alone would not.

**What the tool adds per latitude:**
- **Production:** pressure, temperature, density (`n = N / R_bar`, `rho = n m_bar`) and molar mass. The
  composition is read at each latitude; in these inputs it does not vary, but it should be read per
  latitude in general.
- **The column** (`lm.column`), for `g` and height. It is needed because the lapse rate in K/km needs
  `g` at each latitude, and the product stores neither `g` nor the column off the target.
- **The wind,** from `WindField.wind_at`.

### 2.2 The grid, and its limits

- **Grid:** every 2.5 deg planetocentric. The mesh puts a node one spacing past the target, and the wind
  file stops at 90 deg (`lib/windfield.py:95-108`). So each march ends at 87.5 N or 87.5 S, and the
  section runs from 87.5 S to 87.5 N, with 71 latitudes.
- **Alternative:** `chain` and `outer_loop` accept a given mesh (`transfer.py:430, 1155`), so one
  pole-to-pole march would also work. I propose two marches, because that is the tested path; the
  given mesh has not been run.
- **Cost:** about a minute per march, so 10 marches for the five sections. The anchor matrix is 12
  ordinary runs, of which Lindal's three are SPEC_19's. These figures are estimated from earlier runs,
  not measured for this mesh.

### 2.3 Where it lives

- **The tool:** `src/casspian/tools/section/`, command `casspian-section`. It reads one control file
  naming the anchor's closure or transfer run to copy, the wind and the latitude step. It runs the two
  marches in a work directory, as `casspian-wind-fit` copies runs, and writes one product.
- **The sections:** a new top-level `sections/`, beside `fits/`, one directory per section:
  `sections/<anchor>_<wind>/<name>.toml` and `<name>_section.nc`, with `sections/*/work*/` ignored.
- **The product:** temperature, pressure, density, molar mass, wind, height and lapse rate on
  (latitude, level), on the anchor's isobar labels from 1.3 bar to 0.2 mbar. It also carries the
  march runs' identities and the pressure identity per latitude.
- **The anchor matrix runs** go in `forward/<anchor>_to_<target>/`, as every run does.

## 3. For a ruling

None of these blocks the step. Each one is mine to propose and yours to rule.

1. **Voyager 1's heights:** hydrostatic from 100 mbar, with gravity at kind D's `r0`, and Fig. 9 as a
   check (§1.2), in place of the spec's suggestion. Either costs under 0.01 K, by estimate.
2. **The composition:** "Lindal's 94/6" with or without his ammonia. The figures were computed for
   94/6 without NH3. Voyager 2 egress's pressures were integrated with Lindal's C, which includes it.
   I propose his C as for SPEC_19's "one composition for every input": 94/6 with his ammonia, on each
   anchor's own levels.
3. **The Voyager 1 season:** 13 November 1980 needs a solar longitude, computed by Kliore's method
   (`kliore1980_scalars.toml:14`). Lindal's G and R carry 1981-08-26, and nothing compares the two;
   the SPEC_18 precedent reused them.
4. **sigma_T in kind T:** the Voyager 1 files carry 2 bars each, otherwise NaN. I propose NaN in kind
   T, as SPEC_18 did. The figures' sigma_T bands use Lindal's Fig. 4 curve, as SPEC_19 did. Lindal's
   kind N has a finite `refractivity_uncertainty` regardless; the new anchors' will come out of
   `refrac` the same way.
5. **The section's files:** commit each section's control file. I propose ignoring its product, since
   a section is a leaf that nothing reads and it rebuilds in minutes. The alternative is to commit it,
   as `fits/` commits its wind file.
6. **`tests/step04_0/sweep.py`:** its `KIND_OF` table knows only Lindal's names. It does not rebuild
   Kliore's files either. I propose leaving it as it is, with the new anchors rebuilt by their own build
   files, as Kliore's are.

## 4. Order of work after the ruling

1. The ingest module and the three anchors with their closures (acceptance 2).
2. The section tool, checked first against the SPEC_19 runs (acceptance 3).
3. The five sections and the anchor matrix, then the figures, absolute profiles first.

No regression; Lindal's files untouched.
