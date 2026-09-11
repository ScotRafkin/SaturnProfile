# REVIEW 01, Step 6. `lib.latitude`

Review of `reports/REPORT_01_step6.md` against SPEC_01 v0.7 Step 6. Reviewer: S. Rafkin,
11 September 2026. **Disposition: accepted**, with one rename. SPEC_01 is revised to v0.8.
Rename, rerun, commit, then proceed to Step 7.

## Rulings on section 3

1. **`tol_rad`.** Agreed; radians everywhere in `lib`, no exceptions. The spec named it in
   degrees and that was the mistake. Rename the argument; the acceptance script converts its
   stated degree tolerances at one line.

2. **Seed flattening: the 100 mbar oblateness, 0.09822.** Agreed on the physics (the Fig. 4
   label is defined on the 100 mbar surface) and now stated in the spec. The default of taking
   the surface's own flattening when none is supplied is fine as a fallback, since the
   converged value does not depend on it.

3. **Iteration counts now carry their tolerances** in the spec: three to 1e-3°, four to 1e-4°,
   no more than eight to 1e-9°.

4. **Null's GM, recorded.** Same pattern as Step 5; both sets of numbers are now in the spec so
   the 1.1 arcsec is never mistaken for a defect.

## Notes

The 0.2° label sensitivity passing through almost one for one (0.187°) is the most useful
number in the report. It means the frozen anchor latitude inherits the uncertainty of a figure
label, and nothing in the conversion reduces it. The spec now says so, and the swath in
`lindal_scalars.toml` is the record the later uncertainty accounting will draw on. This will
matter in SPEC_02.

One property of the returned pair worth stating in the docstring: the returned `psi` is the
tilt at the previous iterate, so `phi_c + psi` equals `phi_g` exactly while `psi` differs from
the tilt at the returned `phi_c` by less than the step size times `dpsi/dphi`, which is below
a tenth of the tolerance. The round trip check (2.3e-11°) confirms it does not matter; it should
still be said so nobody later "fixes" it into an extra evaluation.

The `arctan2` seed at the poles, the per-element stopping, the exact zero at the equator, and
the hemisphere mirror are all correct.

## Actions

- Rename `tol_deg` to `tol_rad`; adjust the docstring per the note above; rerun; commit.
- Step 7 next. It freezes the two numbers Steps 5 and 6 could only bracket: the anchor
  latitude and radius with the Smith wind, and the wind geoid equatorial radius against
  60,367 ± 4 km.
