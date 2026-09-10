# REVIEW 01, Step 3

Review of `reports/REPORT_01_step3.md` against SPEC_01 v0.5 Step 3. 10 September 2026.

**Disposition: accepted.** Five of five acceptance checks pass by exact equality, the second
harmonic set round-trips with all six degrees, the rotation rate is bit identical to
2 pi / period, and the reader refuses an unknown convention code by name.

**On the decisions.** All six accepted. In particular: keeping `fitted` (how the number was
obtained) apart from `status` (whether the transcription was checked) is right; recording the
transcription's own prose beside the convention code is the only way a human can re-check the
assertion, and it is welcome; the check values stay out of kind G and the Step 4 script reads
them from the static file; the narrow control parser in `tools/gravity` is the correct scope and
is absorbed by `lib.control` at Step 9.

**On the defects and questions.**

1. `J6` uncertainty: agreed, and fixed at the source. `null1981.toml` now carries `nan`, and
   SPEC_00 v0.7 §5 makes it a rule: an uncertainty the source does not state is NaN, never 0.0,
   in transcriptions as much as in netCDF. Rerun the gravity tool after the acceptance commit;
   `J_uncertainty` will then read `[1.8e-05, 3.8e-05, nan]`.
2. Provenance for a label coordinate: `index` is added to the vocabulary (SPEC_00 v0.7 §5) for
   `degree`, `species_name`, `surface`, and the like. Use it here and in kinds C and D.
3. `J` is `derived`. SPEC_00 v0.7 §5 now defines the two: `derived` is obtained from
   measurements by a fit or computation with no added physical assumption (a fitted harmonic,
   a hydrostatic pressure); `inferred` is deduced through a declared model assumption (an
   ammonia abundance from an opacity with an assumed line shape). Every later kind G file
   follows.

**Then:** apply the three changes (they touch the tool's provenance strings and the static
file only), commit Step 3, rerun both tools for clean products, set STATE row 3 to `accepted`,
begin Step 4.
