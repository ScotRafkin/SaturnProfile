# REVIEW 04, Step 1. The reference surface through the anchors, and the wind on the model's geometry

Review of `reports/REPORT_04_step1.md` against SPEC_04 v0.9 Step 1. Reviewer: S. Rafkin,
22 September 2026. **Disposition: accepted with changes.** The two failing checks fail on the
cylinder-extended wind, whose construction the specification did not state and whose stated
values were the reviewing agent's and were wrong; decision P at SPEC_04 v0.10 fixes the
construction and places it at Step 2, where the columns it needs exist, and the two checks move
there with it. Everything else passes and was verified. The rulings are in SPEC_04 §12.

## What was verified independently

The reference surface through the anchor: the report's 60,366,999.9999 m at the equator,
54,420,643.4644 and 54,449,383.8347 m at the poles, and 60,128,612.9664, 59,492,075.6754,
57,053,249.5161 and 55,675,128.9158 m at 10°, 20°, 45° and 60° N are the reviewing agent's
independent march under the linear wind rule to 0.04 m (60,367,000.03, 54,420,643.49,
54,449,383.84, 60,128,613.00, 59,492,075.7, 57,053,249.5, 55,675,128.91). The anchor's wind,
2.16708095 m/s, is the linear value between the file's nodes at 30.5° and 31.0°, and the slope
of that interval, 6.0357 m/s over 0.5°, is the report's −691.638 m/s per radian. The closure
production through `produce` with the wind array being bit-identical to the registered product
is consistent with the reviewing agent's own closure, which the array path cannot change when
every element is the scalar.

On the cylinder-extended wind the reviewing agent rebuilt the construction from scratch after
the report. Pinned on the gauge isobar with the anchor's field-line height, the coding agent's
construction gives 7.636 m/s at the top of the anchor's column; pinned on the gauge isobar with
the radial column, 7.465; pinned on the 1 bar level with the radial column, 9.490. None of
these, and none of the same three with a PCHIP wind, gives the specification's 9.371 or its
−0.133 at the bottom, and the reviewing agent no longer holds the code that produced those
numbers. They are withdrawn. The report's inference that the stated values differ from its
construction by one factor of about 1.30 was the useful clue: it is, to a few percent, the
height above 1 bar over the height above the gauge at the top of the profile, which is what
pinning on the 1 bar level does to the displacement.

## Rulings

1. **The construction (finding 2, decision P).** A kind W file's three parts make
   `u_total(φ, p_ref) = u_reference(φ)` an identity, so the extension is pinned on the file's
   reference level, 1 bar, not on the gauge isobar; the report's construction had the
   reference-level wind at 0.917 m/s at the anchor, which is not "the same reference-level
   wind". The geometry is the library's radial columns under the flat-isobar map, and the
   reference surface and the file are a fixed point (the surface is marched with the file's own
   wind at the gauge), converged in the script. Since the columns are Step 2's, the file is
   built there, and check 9's cylinder half and check 14 move to Step 2 with the values the
   reviewing agent measured under the ruled construction: on the anchor's column 9.490 m/s at
   the top level, 3.717 at the gauge, 2.168 nearest 1 bar, 1.926 at the bottom; the anchor's
   `Φ` lower at the top by 322.3 m²/s² and higher at the bottom by 18.0; at 10° N 427.09,
   345.98 and 326.40 m/s at the top, gauge and bottom against 329.46 on the reference level.
   The Step 5 altitudes under this wind at 10° N are restated (416,300, 99,290, −15,456 m).
   Nothing is loosened: the checks are relocated with the construction they test.
2. **Finding 1** (the march grid overshooting the pole): the snap and the refusal in
   `through_anchor` are right, `wind_geoid` is left alone. Accepted.
3. **Finding 3** (bit-identity measured, not guaranteed): recorded. Accepted.
4. **Finding 4** (check 0 of `accept_step04_0` vacuous on a clean tree): the proposal is
   adopted, a `-dirty` copy under `reports/step04_0/` as the check's subject. Apply at this
   step's acceptance commit and record it in the report.

Decisions 1 to 10 and 12 are accepted as reported. Decision 11 is superseded by decision P.
Section 5b is noted for the Monte Carlo wrapper.

## Order of work

1. Remove check 9's cylinder half and check 14 from the Step 1 acceptance (they are Step 2's
   at v0.10), with a line in the report saying so; apply ruling 4 to `accept_step04_0`; rerun
   both; refresh the report's head, section 3 and section 5a.
2. Commit the author's documents (SPEC_04 v0.10, this review, `STATE.md`) in their own commit.
3. The acceptance commit for Step 1. Push.
4. The sweep: the closure product rebuilt on the clean tree with `u_column_ms`, its hash
   recorded in the report's section 6 in the row format, `step03_4` and the other suites that
   read it rerun on the clean product and recorded. The record commit.
5. Set `STATE.md` Step 1 to accepted with the commits. Step 2 proceeds, and its acceptance
   script builds the cylinder-extended wind by decision P as its first deliverable beyond the
   mesh.
