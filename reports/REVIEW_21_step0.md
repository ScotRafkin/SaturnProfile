# REVIEW 21, Step 0. The general error budget, and hypothesis 4

Reviewing agent, 9 November 2026. Report: `reports/REPORT_21_preexecution.md`. Specification:
`docs/specs/SPEC_21_General_Error_Budget.md` v0.2, with rulings recorded at v0.3.

**Accepted. Proceed with Step 1.**

- **The prior's construction is right:** the Matérn 5/2 temperature covariance differentiated to the
  gradient, conditioned to vanish where the balance and the sphere require it, divided by `w`, with
  the equatorial limit taken explicitly. Draws are finite and smooth at the equator.
- **It does what was wanted toward the poles** (19.5 m/s per unit `ln p` at 74 deg, against SPEC_19's
  40) **and what the physics implies at the equator:** 600 at 110 mbar, 1500 at 10 mbar. That is
  σ_T over a 5 deg length turned into curvature and divided by a small `w'(0)`. It is honest about
  how little temperatures say there; the constraint will come from wind observations (SPEC_23).
  §2's values stand as ruled, to be revisited on results.

**Rulings on §4:**
1. **Equatorial looseness:** keep §2's values; add the wind to acceptance 1(d) (the largest shear and
   the wind change from 10 to 1000 mbar of the 20 draws at the equator). No separate equatorial
   amplitude: that would be a rule for one region.
2. **Conditioning at the poles** as well: accepted, physics.
3. **`w` held** at the base wind and the anchor's `T̄(p)`: accepted.
4. **`σ_rep`'s correlation:** Matérn 3/2 over 0.5 scale height, an input. Accepted.
5. **Huber:** on whitened residuals within each group, threshold 2. Accepted.
6. **The deepest Voyager 2 egress row:** the nearest bar's value. Accepted.
7. **IRIS rows:** the equatorial omission is dropped; the gap rule is kept, since a row where the fit
   only bridges is not a measurement.
8. **IRIS's uncertain level:** one level label serves every row at that level, so its error is fully
   correlated across them, not independent. Enter it as a rank-one term per level: the vector of the
   operator's differences between 110 and 150 mbar, scaled for a uniform range (`d / √12`),
   recomputed each pass.

**The ingest (§3):** accepted as proposed. Each profile's uncertainty at ingest in its kind T;
Lindal's in `occul_data/lindal/lindal_temperature_uncertainty.csv` beside his untouched files, named
in the control file; observations read from kind T.

**Order of work:** Step 1, about 4 to 5 hours of acceptance. No regression; Lindal's files untouched.
