# REVIEW 03, Step 3. The season identifier, the run directory and its inputs, the run namelist, kind `profile`, and the reader warning

Review of `reports/REPORT_03_step3.md` against SPEC_03 v0.10 Step 3 (deliverables 0 to 4).
Reviewer: S. Rafkin, 15 September 2026. **Disposition: accepted**, with the author's decision
on finding 1 (one date for the whole chain) applied before the acceptance commit and the
specification amended at SPEC_00 v0.19, SPEC_01 v0.28 and SPEC_03 v0.12. Step 4 proceeds after the acceptance commit, the
sweep and the suites of the report's section 5.

## What was verified independently

The SHA-256 of all twelve working-tree files (the seven reduction products, the candidate kind
N, the four run inputs) was recomputed and equals section 7 and check 16 in every digit. The
Step 0 product and the candidate kind N were compared group by group: the 23 root variables and
the radius coordinate are bit-identical, every embedded group's variables are identical, and the
attributes that differ are the writer globals, the hash and path attributes, and the season
attributes the step adds; the wind's `epoch` went from prose to `1981-08-25`; the
`reduction_record` differs only in the commit and the seven input hashes, as check 15 says. A
scan of every path-bearing attribute in every group of the twelve files found no drive letter
and no leading `/`. The closure comparison was repeated with an independent implementation of
the drop rule: all four run inputs are identical to the embedded copies after the drop, and the
attributes that actually differ before it are a subset of the list ruling 2 enumerates
(`decomposition_geometry` does not differ; dropping it anyway is the rule). The season and role
attributes read from the files are those of checks 1 and 5. `asin(sin 26.73° sin 18.2°)` =
8.0759°, 0.016° from the transcribed 8.06°; the bottom-row pressure term 0.005 / 1298.48 =
3.85e-6; 50 m over 31.9 to 57.9 km is 8.6e-4 to 1.57e-3. `lib.control` was read: the loader
refuses in the order decision 7 states, `[grid]`, `[numerics]`, `[estimation]`,
`datum_isobar_Pa` and `transfer` are refused as not implemented with a message distinct from
the unknown-key refusal, the season equality is exact, the gauge match is exact with the nearest
level named, and the closure comparison strips the drop list and every `sha256:` value from
groups and variables alike. No em dash or en dash in any of the fourteen new or changed files
checked. The synthetic F6 was viewed: the residual, the band, the boundary and gauge markers
render as described.

## Rulings on the findings

1. **The epochs of G and R.** Correct, and the contradiction is the specification's: SPEC_00
   v0.17 §5 makes `epoch` an ISO date and SPEC_01 v0.25 Step 3 said "as already recorded".
   The author's decision (15 September 2026): every file of the chain, reduction and closure
   inputs alike, carries `epoch = "1981-08-26"`, the date of the occultation; a file whose
   content has no date of its own (G, R, W, C) takes the date of the observation it serves, and
   its source's own dating goes in `epoch_note` when it is not already stated elsewhere in the
   file. No per-source dates. The `[gravity]` and `[rotation]` build-file sections carry `epoch`
   beside `role` and the tools write it; the prose the files carried becomes `epoch_note`; the
   wind's properties file `[epoch]` says 1981-08-26 with its existing note; the reader parses
   `epoch` as an ISO date or, on a profile, the SPEC_00 §5 string for a season declared without
   a date. Applied before the acceptance commit so that the chain is rebuilt once: the rebuild
   and the acceptance are rerun as scripted and the report's checks 1, 2, 15 and 16 and section
   7 are refreshed. SPEC_00 v0.19, SPEC_01 v0.28 and SPEC_03 v0.12 carry the rule.
2. **The registered kind N refused until the sweep.** As intended; no action.
3. **Kind C's `source_statement` rule declared and not enforced.** A reader defect. Enforced at
   Step 4 in the same commit series (no product changes, no rebuild); SPEC_03 v0.11 Step 4.
4. **The wind's `solar_longitude_source` wording.** Accepted; it matches the raw bundle's and
   names the method and its replacement.
5. **The F6 envelope read column-wise.** Correct: the last printed figure of a value cannot be
   recovered from a float. Restated at v0.11 as the finest precision any value of the column
   shows, 0.01 mbar and 0.1 km for Table I.
6. **F5 presentation.** Deferred to Step 4's figures from the production, with one instruction
   now: the gauge marker is labeled with the declared `gauge_isobar_Pa`, the tabulated value,
   since the gauge is defined by the tabulated level, not by the hydrostatic pressure found
   there.

## On the decisions

Decisions 1 to 11 are accepted as reported. Three are recorded in the specification at v0.11:
the raw bundle's role is always `reduction` (decision 4; the bundle is the transcription of a
source and has no build-file section); the point-latitude marker is
`latitude_planetocentric_absent_meaning`, the SPEC_00 §5 form, and the specification's shorthand
is corrected (decision 5); the closure refusal names every difference (decision 8), which is
better than the first-only wording. Decision 2 (a path on another drive refused rather than
recorded absolute) is the right refusal for a rule that admits no absolute path.

## A finding on the specification, Step 4

The Step 4 acceptance asks that the residual drawn in F6 sit inside the envelope band at every
snapped level. The band is the leading-order budget, ±50 m over the local scale height, 1.3e-3
at 39 km, while the same section's expected values put the residual between 10 and 100 mbar at
±2.6e-3. As written the clause would fail on numbers the specification calls expected. Restated
at v0.11: the band is drawn, the levels outside it are counted and reported with their values,
and the numeric bounds of the acceptance are the acceptance. The bounds do not move.

## Open with the author, not blocking

1. Ruling 5 of §8, committing `forward/*/inputs/*.nc`, remains open; the `.gitignore` rule
   stands until it is made.
2. The `solar_longitude_source` attributes name "the reviewing agent" as the source of the
   computed `Ls`; the values and the method stay whatever the wording, and it is replaced when
   the season tool exists.

## Order of work

1. Commit the author's documents (SPEC_00 v0.19, SPEC_01 v0.28, SPEC_03 v0.12, this review,
   `STATE.md`) in their own commit.
2. Apply finding 1 as ruled (the `epoch` key in the G and R build-file sections of both build
   files, `epoch_note` from the `data_static` prose, the wind's properties file dated
   1981-08-26, the reader parsing `epoch`); rerun `rebuild.sh` and the acceptance; refresh the
   report's checks 1, 2, 15, 16 and section 7. The report keeps its number; the change is
   noted at its head.
3. The acceptance commit for Step 3. Push.
4. The clean rebuild of the chain, the candidate product and the run's inputs from the
   committed tree; the sweep, with the hashes recorded in section 7 in the row format the
   Step 02_1 suite reads and the `-dirty` scan reported.
5. The suites of section 5 on the swept products (02_1, 02_4, 02_5 with its changed checks,
   02_6, 03_1, 03_2), results in section 8. The sweep record commit.
6. Set `STATE.md` to accepted with the commits. Step 4 then proceeds when the author says go.

## Sweep verified (15 September 2026, after the acceptance commit)

The refreshed report was read against the swept products on disk. Every SHA-256 in section 7
(the manifest, the eight reduction files, the four run inputs) was recomputed and matches. Every
one of the twelve netCDF files carries commit `6ae113e` and none `-dirty`. Every file carries
`epoch = "1981-08-26"`; gravity and rotation carry the former prose as `epoch_note`, the wind its
imaging note. The swept kind N is bit-identical to the Step 0 product in every variable of the
root and of every embedded group; its `reduction_record` differs only in the commit and the seven
input hashes. The deferred suites pass on the swept products: Step 02_1 (6 of 6, rerun after its
stale key table was replaced by the gravity tool's own `SECTION_KEYS`), 02_4 (9), 02_5 (7 with
its changed checks), 02_6 (7), SPEC_03 Step 1 (9), Step 2 (8). The committed namelist loads on
the swept products without relaxation, which closes finding 2.

Step 3 is closed. The remaining action for the coding agent is the sweep record commit (the
refreshed report is an uncommitted change in the working tree). Step 4 proceeds when the author
says go.
