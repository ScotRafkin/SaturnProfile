# REVIEW 04, Step 3. The kernels

Review of `reports/REPORT_04_step3.md` against SPEC_04 v0.11 Step 3. Reviewer: S. Rafkin,
23 September 2026. **Disposition: accepted with changes**, the regression's table to be filled in
before the acceptance commit. The two failing checks fail on the specification's bounds, which
were written for analytic derivatives and are restated at v0.12 for the interpolant's; the
report's headline conclusion, that decision L and Eq. A15 are in tension at order one and that
the cylinder identity is lost at the 7.6e-2 level, does not survive the measurement the transfer
actually makes, and the reviewing agent's independent kernel puts the effect at 1.1e-3. The
rulings are in SPEC_04 §14 (v0.12).

## What was verified

The kernel under the closure wind is the closed form `2 Ω_abs cos φ u'/g` to 4e-17, with the
largest `|S/g|` 0.0976 per radian to 10° N and 0.0539 to 60° N against the reviewing agent's
0.0977 and 0.0540; `I` is exactly zero at the gauge node of every column; the composition term
meets its closed form to 2e-10 and is exactly zero for the run's uniform composition. Those are
the deliverables of this step and they are right.

On the cylinder-extended wind the reviewing agent rebuilt the kernel independently on the coding
agent's own Step 2 file: the bilinear interpolant of the file, its piecewise-constant partials,
the conversion to fixed radius through the column map, the columns marched under the file's
wind, a 0.05° mesh. With the flat-isobar map continued at the slope of its last interval, as
Step 2 decision 8 states, the largest `|S/g|` on the span is 0.050 per radian; with the map
clamped at the anchor's ends, as the Step 3 script does, it is 0.136, the excess at the clamped
rows and their neighbors. The report's 0.209 is of the second kind. The physically relevant
quantity, the line integral of `S/g` along each of the anchor's isobars from `φ_c` to 10° N,
gives `Δ ln N` of +1.2e-4 at the top level, −5.2e-4 at 10 mbar, −1.7e-4 at the gauge and
−1.1e-3 at the bottom; the vertical integral `I` and the isobar shift it implies reach 224
m²/s². Along the worst isobar the residual changes sign 55 times in 417 nodes, with an rms a
fifth of its maximum: the signature of cell-boundary spikes whose integral is small.

## Rulings

1. **Finding 1, the solid body.** Correct; "to round-off" was for analytic derivatives. The
   bound is `|S/g|` below 1e-4 per radian (measured 3.4e-5).
2. **Finding 2, the map's continuation.** The script clamps the map where Step 2 decision 8
   continues it at the slope of the last interval, which is slope-continuous. The eightfold
   degradation at the ends, the coinciding pressure nodes of finding 4 and most of the cylinder
   maximum are the clamp. The script follows decision 8.
3. **Finding 3.** The difference form for the composition term is the rule. Accepted.
4. **Finding 4.** Synthetic states on the mesh's own geometry: right. The NaN is the clamp's;
   strictly increasing nodes are the schema's business; no guard (decision N). Accepted.
5. **Finding 5, the cylinder-extended wind.** Right about the mechanism, wrong about the
   measure and the size. `max |S/g|` times the span is not `Δ ln N`, and the maximum of a
   sign-alternating residual does not fall under refinement while its integral does. Decision Q
   (v0.12) states the truncation floor of a gridded wind hypothesis: on the 0.5° grid, up to
   1.1e-3 in `ln N` and 224 m²/s² in the isobars over 21° of latitude, against 1.1e-2 and
   31,000 under the closure wind; a property of the data's resolution, reported and carried
   once in the manuscript, not a change to decision L or P. Check 9 is restated as the line
   integrals, bounded at 2e-3 and 300 m²/s², with the maximum reported and the latitude
   refinement reported without an order claimed; the Step 4 and Step 5 cylinder "identity"
   tolerances are restated the same way.

Decisions 1 to 6 are accepted as reported. The report's own correction of its `d ln N` column
is noted; the replacement was still the wrong quantity, and the line integral is the right one.

## Order of work

1. Follow decision 8 for the map in the acceptance script; restate checks 5 and 9 as v0.12
   Step 3 states (the solid-body bound; the line integrals per isobar, the implied shift, the
   maximum reported, the sampling refinement reported); rerun; refresh the report's head,
   section 3 and findings 2, 4 and 5 with the remeasured values, which should close on the
   reviewing agent's within the interpolant's differences.
2. Fill the regression table when the run finishes.
3. Commit the author's documents (SPEC_04 v0.12, this review, `STATE.md`) in their own commit.
4. The acceptance commit for Step 3. Push. No product changes; `STATE.md` to accepted. Step 4
   proceeds, with its cylinder-wind identity at the v0.12 tolerances.

## Refreshed report verified (23 September 2026)

Nine of nine under v0.12, the regression complete at every reference count. The script's clamped
inverse map was the source of the eightfold degradation, the coinciding pressure nodes and most
of the first filing's cylinder maximum, as ruling 2 said; corrected, check 4 measures the whole
mesh. Check 5 at 3.2e-5 per radian and check 9 at 4.7e-4 in `Δ ln N` and 130 m²/s² in the
isobars, inside the v0.12 bounds. The coding agent's aggregate values agree with the reviewing
agent's (largest `|S/g|` 0.054 against 0.050; shift 130 against up to 224); its per-level line
integrals differ by factors of five to nine and one sign, which the reviewing agent accepts as
the sampling of a sign-alternating integrand and not as a defect in either kernel, and SPEC_04
v0.13 says so: the per-level values are a size class and the bounds are the check. The
report's candor on its own `d ln N` column is noted. **Step 3 is accepted.** Order of work from
item 3: the author's documents commit (SPEC_04 v0.13, this review, `STATE.md`), the acceptance
commit, push, `STATE.md` to accepted, Step 4.
