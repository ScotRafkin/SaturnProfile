# REVIEW 01, Step 2

Review of `reports/REPORT_01_step2.md` against SPEC_01 v0.5 Step 2. 10 September 2026.

**Disposition: accepted.** Eight of eight acceptance checks pass by exact equality, the
`scalars` groups reproduce every entry of `lindal_scalars.toml` with its `value_source`, the
nine NH3 values sit at the right pressures with the two gaps left as NaN, and the housekeeping
commit and the Step 1 naming correction are done as asked.

**On the decisions.**

1. Decimal parsing before the single rounding to float64: accepted as implemented, and not made
   a rule. The author's position is that a last-bit difference (250.99999999999997 against 251.0)
   is not an issue for this work; floating-point error is unavoidable downstream and
   insignificant at this level. The tool's choice of the double nearest the printed decimal is
   harmless and costs nothing, so it stays; no other tool is required to do the same.
2. Nested TOML tables as nested groups: accepted. The bundle keeps the shape of the
   transcription it came from, which is the point of a raw bundle.
3. `input_hashes` (what was read) distinct from `raw_sources` (what is on record): correct.
4. Attribute order: noted, not a defect.
5. The three extra `table1` group attributes and the coordinate attributes on `pressure_Pa`:
   welcome.
7. `-dirty`: as the rule expects. Rerun after the acceptance commit.

**On the questions.**

1. Nested groups: settled above.
2. Second reader: **withdrawn by the author.** The scalars transcription is taken as checked; if
   an error is found later, the TOML is corrected and the chain rerun from stage one, which is
   what the pipeline is built to allow. `notes.md` is amended to say so, and the raw bundle is
   frozen on this acceptance.

**Then:** commit Step 2, rerun `casspian-lindal-raw` so the bundle carries the clean commit, set
STATE row 2 to `accepted`, begin Step 3.
