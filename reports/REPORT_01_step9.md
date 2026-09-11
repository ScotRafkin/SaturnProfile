# REPORT 01, Step 9. Stage two of the Lindal tool, kinds T and D and the manifest

CASSPIAN Saturn atmosphere reference model. Report of the coding agent.

Date: 11 September 2026. Specification: `SPEC_01_Lindal_Tool_Chain.md` v0.17, Step 9, against
`SPEC_00_Architecture_and_Data_Files.md` v0.8. Status: reported, awaiting review.
**Six of six acceptance checks pass. `occul_data/lindal/` is complete.**

**The working tree is not committed**, per the commit on acceptance rule. Four products carry
`-dirty` and are to be rebuilt after the acceptance commit; see finding 4, which is about more
than this step.

---

## 1. What was built

**`src/casspian/tools/lindal/build_inputs.py`**, entry point `casspian-lindal-inputs`, and the
`[stage_two]` section of `lindal_build.toml`. Three products:
`occul_data/lindal/lindal_thermo.nc` (kind T), `lindal_geodesy.nc` (kind D), and
`lindal_reduction.toml` (the manifest of SPEC_00 section 7.1).

The tool reads only the raw bundle and the control file. It separates the transcription into the
two standard kinds, which is the whole of its job: the raw bundle holds the profile and the
fitted geoid together because that is how the paper states them, and SPEC_00 section 6.3 keeps
them apart in the model inputs because the retrieval and the geoid fit are different
measurements with different error budgets. Kind T carries the profile and no geodesy, kind D
carries the surfaces and no profile, and each records which raw bundle it came from.

## 2. Decisions

1. **The four companion files are verified, not rebuilt.** SPEC_01 Step 9 item 3 allows either.
   Verifying is the safer reading: rebuilding them here would mean this tool could silently
   change a file another step is accountable for. The tool refuses with a message naming the
   tool to run if any of G, R, W or C is missing, and it reads each one under its kind, so a
   malformed companion stops the step rather than reaching the manifest.

2. **`height_datum` is a sentence, derived from the Table I footnote.** SPEC_00 section 6.1
   wants the datum stated; the raw bundle carries it only inside the footnote, which says the
   altitude is "relative to the 1-bar pressure level". The attribute says that in those terms.

3. **The manifest is written as text, not serialized from a dictionary.** It is a file a person
   reads and edits, so the alignment, the section order and the two explanatory comment lines
   are written deliberately. The tool refuses before writing it if any named file is missing or
   lacks the directory prefix.

## 3. Findings

### 1. The fit residual is stated as a range, so the per surface variable is NaN

SPEC_00 section 6.3 gives kind D a `fit_residual_m(surface)`, "standard deviation of the fit
where the source gives it". Lindal gives it as a range across his fits, not a value per surface:
p. 1141, "the standard deviations of these fits range from 3 to 4 km". There is no way to say
which surface got 3 and which got 4.

The variable is therefore **NaN at both surfaces**, which SPEC_00 section 5 makes the honest
statement of no information, and the range is carried as `fit_residual_range_m = [3000, 4000]`
with the quotation in `fit_residual_source`. Flagged because a reader expecting a number will
find a NaN, and because SPEC_02's uncertainty accounting will want the range rather than a
fabricated per surface value.

### 2. The gravity citation exists only inside a prose attribute

SPEC_00 section 6.1 wants `source_gravity_citation` on kind T, and its example is
"Null et al. 1981; Campbell 1984 pers. comm." The raw bundle has no key holding that. It appears
only inside `scalars/geodesy.fit_inputs`, a long sentence covering the gravity, the winds, the
pole vector and the rotation together.

Rather than cut a substring out of a sentence, which would break the moment the sentence is
reworded, the attribute carries the **whole `fit_inputs` statement**. It is correct but
overlong, and `source_rotation_system` and `source_wind_citation` beside it are crisp because
the raw bundle has dedicated keys for those. **Proposing `lindal_scalars.toml` gain a
`gravity_citation` key** under `[geodesy]` or a new `[gravity_used_by_source]` table, in the
same shape as `[wind_used_by_source]`.

### 3. The longitude is a swath and section 6.1 asks for a scalar

The source states the ingress swath as running from 186.8 to 185.7 degrees System III
(p. 1138), and the raw bundle carries the pair. SPEC_00 section 6.1 gives `longitude_deg` as a
single value, its example being 186.8. The file carries `longitude_deg = 186.8`, the end where
the ingress data begin, and keeps the pair beside it as `longitude_swath_deg`, in the same shape
as the latitude, which the specification already treats this way. Recorded so the choice of end
is visible rather than implied.

### 4. The rebuild rule does not cover a product a later step rebuilds in passing

Check 1 reports `lindal_gravity.nc` and `lindal_rotation.nc` as **DIRTY**, although Step 3 was
accepted and they were rebuilt clean at the time. They were rerun during Step 7, when the
`control.py` path resolution changed and I reran both tools to confirm the change had not
altered their products. That rerun happened from an uncommitted tree, so the clean commit they
carried was replaced by a dirty one, and nothing since has rebuilt them.

The commit on acceptance rule says a step's tool is rerun after that step's acceptance commit.
It does not cover a product rebuilt incidentally by a later step, which is what happened here.
**Proposing that the rule be stated as a property of the directory rather than of a step**: at
the end of every acceptance commit, any product under `occul_data/` carrying `-dirty` is
rebuilt, whichever step wrote it. The check that would enforce it is the one above, which is
cheap and now exists.

For this step the four dirty products are `lindal_thermo.nc`, `lindal_geodesy.nc`,
`lindal_gravity.nc` and `lindal_rotation.nc`; all four will be rebuilt after the acceptance
commit.

### 5. The manifest tolerance is in degrees while `lib` is radians throughout

SPEC_00 section 7.1 names the manifest key `fixed_point_tolerance_deg`. The Step 6 review ruled
radians everywhere in `lib` with no exceptions, and renamed `tol_deg` to `tol_rad` for exactly
this reason. The manifest is a control file rather than `lib`, so a degree valued tolerance
there is defensible: it is the number a person reading the manifest would rather see, and
`refrac` converts it at one line. It is written as the specification names it. Raised only so
the inconsistency is deliberate rather than inherited, since SPEC_02 is where it lands.

## 4. Acceptance results

Run by `reports/step9/accept_step9.py`; full output in `reports/step9/output.txt`.

| Check | Measured |
|---|---|
| All six input files validate under their kinds | **Pass**, all six. Written by five different tools; four carry `-dirty`, finding 4. |
| Kind T has exactly the three profile variables and their three companions | **Pass**, the six variable names are exactly the expected set. |
| No `x_*` or geodesy content in kind T | **Pass**, none. |
| The three uncertainty companions present and NaN | **Pass**, all three all NaN, the source stating none. |
| `latitude_planetographic_deg = 36.3` with `value_source` | **Pass. 36.3**, `value_source = "Fig. 4 label on the Voyager 2 ingress profile: 36.3 N"`, swath [36.3, 36.7], uncertainty 0.2 (range). |
| Kind D on the 1e4 Pa surface: 60,367 / 54,438 / 0.09822 | **Pass**, exactly, with +- 4 km, +- 10 km and +- 0.00018. |
| Kind D on the 1e5 Pa surface: 60,268 / 54,364 / 0.09796 | **Pass**, exactly, with the same uncertainties. |
| The manifest parses and every path in it exists | **Pass.** Six sections, six inputs, all present, all prefixed `lindal_`, none pointing outside the directory. Anchor isobar 10000 Pa, geoid anchor `radius_polar_m` on the 10000 Pa surface, convergence 1 m. |
| The directory listing matches SPEC_00 section 2.2 | **Pass.** Eight files at the top level and four in `raw/`, none missing and none unexpected. The only file without the `lindal_` prefix is `notes.md`, which section 2.2 lists unprefixed itself. `lindal_refractivity.nc` is absent as expected: `refrac` writes it at SPEC_02. |

Six of six pass. No em dash or en dash appears in any file written in this step.

## 5. What this completes

`occul_data/lindal/` now holds every file SPEC_00 section 2.2 lists except the product, and the
manifest names the six the reduction reads. SPEC_01 is finished at this step; SPEC_02 specifies
`refrac`, whose first acceptance numbers are the wind included `phi_c` and `r0`.

Two numbers from the Step 7 review are queued for that uncertainty accounting and are recorded
here so they are not lost between specifications: the **28.7 km polar asymmetry** the wind
produces, against the roughly 10 km in the same sense that Lindal's Fig. 9 caption allows, and
the **+- 19 km spread** between polar anchor choices, which is the size of the anchoring
uncertainty until the observed radii of his Fig. 9 are brought in. Both live in
`lindal_geodesy.nc` in part: `polar_asymmetry_note` carries the caption, and
`radius_polar_uncertainty_m` carries the +- 10 km.

---

## 6. Addendum: the directory sweep, and SPEC_01 closed

Added 11 September 2026, after `REVIEW_01_step9.md`. Two changes were made, then the sweep the
review asked for was run.

**1. `source_gravity_citation` reads the new table.** The author added
`[gravity_used_by_source]` to `lindal_scalars.toml` in the shape of `[wind_used_by_source]`, so
kind T now carries `"Null et al. 1981; Campbell 1984 (personal communication)"` with its
Appendix quotation in `source_gravity_citation_value_source` and the pole vector references in
`source_pole_vector`. The prose it used to be cut from is still there as `source_fit_inputs`,
the source's own statement of everything that went into the geoid fit. Finding 2 is closed at
the transcription, which is the right place.

**2. The manifest names the anchor rule.** SPEC_00 v0.9 adds `anchor_rule` to the `[geoid]`
vocabulary, and the tool writes it from the control file. The reduction now has, in one place,
the same statement Step 7 built its surface with:

```toml
[geoid]
anchor_surface_Pa = 10000
anchor_quantity   = "radius_polar_m"
anchor_rule       = "mean_polar_radius"   # SPEC_01 v0.16: how the one constant of Eq. B3 is fixed
convergence_m     = 1
```

This is the item the report could not have raised: the `[geoid]` section predates v0.16, and
without the rule the reduction could have frozen `r0` on a differently anchored surface than
the one Step 7 reported, which is a 19 km error and would have been invisible.

**3. The sweep.** The scalars edit changes `lindal_scalars.toml`, which the raw bundle hashes,
so it changes the bundle and everything downstream of it. Every product was rebuilt in
dependency order after the acceptance commit: the raw bundle, then G, R, W, C, then T, D and the
manifest. All now carry the clean commit **`ece58d24733387e8706e2680d26178d435c8f010`**, and the
Step 9 acceptance reports no dirty product.

| product | kind | SHA-256 |
|---|---|---|
| `raw/lindal_raw.nc` | raw | `3054cc6cf0b54905ee3d7d67bbe28ad12bac18288e1fcd53893957d54d9171eb` |
| `lindal_gravity.nc` | gravity | `6db0129cefd01bc9c151020e604cf4563a82381e61c97381d0a7ee383826911e` |
| `lindal_rotation.nc` | rotation | `a5c72017a423c50f2a952ddcbc4ccf777c1e94683fa6965ee96f7f43d949a424` |
| `lindal_wind.nc` | wind | `18d8cdae4dd8faa5cdcea40d25e6a9dfce4391421b3af4205d2b377f1150e7b9` |
| `lindal_composition.nc` | composition | `b882abcac3bf09cb4734fabb2fa6a463f4e9c1bde9428747b7f9b433555165cf` |
| `lindal_thermo.nc` | thermo | `ff9c12991399f94c08caaae80b1a6a91f45a45feaebc4d059ad563b45a720e85` |
| `lindal_geodesy.nc` | geodesy | `16e0c7d551c3026ce5623141b1a2f86187c013a2dbf929f3e330f79f410f787e` |
| `lindal_reduction.toml` | manifest | `d908015605b2293442da8171d580d1ee4b2793cf6b8ac6406d09cfb0500486b9` |

**Regression after the sweep:** ten suites, all passing. 6, 14, 8, 5, 7, 7, 6, 8, 7, 6.

**This closes SPEC_01.** `occul_data/lindal/` holds every file SPEC_00 section 2.2 lists except
`lindal_refractivity.nc`, which `refrac` writes at SPEC_02. The manifest names the six inputs
the reduction reads, every path in it exists, and every product carries a clean commit and the
hash of everything it was made from.
