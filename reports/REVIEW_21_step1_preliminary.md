# REVIEW 21, Step 1, preliminary. The general error budget, and hypothesis 4

Reviewing agent, 9 November 2026. Report: `reports/REPORT_21_step1_preliminary.md`. Specification:
`docs/specs/SPEC_21_General_Error_Budget.md` v0.3, amended to v0.4.

**Not accepted; Step 1 reruns under v0.4.** What the run built is accepted as reported (§2 there).

- **What worked.** The fit converged in 4 passes; IRIS 110 mbar fell from 2.32 to 0.79 of its expected
  value; no systematic part was needed and Huber downweighted nothing. Under hypothesis 4 the anchor
  transfers sit inside the stated bands over more of the column, and the warm bump near 10 mbar from
  3.0 S is gone. The temperature draws are physical (1d).
- **What failed, one cause.** Dividing the temperature-gradient covariance by `w` leaves the shear prior
  without a physical limit at the equator: draws change the equatorial wind by thousands of m/s
  between 1000 and 10 mbar. The same breadth is why 1b cannot recover the equatorial layer, why the
  synthetic fits show a spurious equatorial shear column, and much of why the overall ratio is 0.14.
  The real fit is not harmed, but its stated uncertainty and any Monte Carlo would be.
- **Why a temperature bound alone does not cure it.** The draws' temperatures are already physical. At
  the equator the shear follows the curvature of temperature in latitude, and that curvature, set by
  `σ_T / ℓ²`, is what runs away. A clip on the shear (150 m/s per unit `ln p` was considered) would stop
  it without being right (author).

**Rulings (author, 9 November):**
1. **A deformation-radius floor on the latitude correlation length.** A balanced structure cannot be
   narrower than the deformation radius; narrower features are waves. The length used at each latitude
   is the larger of the input length and the local deformation radius, `L_R = c / |f|` away from the
   equator and `√(c a / 2Ω)` at it (a standard form joining the two), with `c = N × (the prior's
   vertical correlation length in height)` and `N` from the anchor's own profile. General physics; it
   binds near the equator, where it is about 8 to 9 deg for a structure one scale height deep. The
   covariance becomes non-stationary in latitude; the gradient covariance may be formed numerically on
   the fine grid before conditioning and division by `w`.
2. **No superadiabatic profile.** A one-sided constraint on every delivered profile: a penalty that is
   active only where the lapse rate would exceed the dry adiabat for normal hydrogen (the stricter), with
   `c_p` from the run's composition. Applied in the fit at each pass and reported where active; the
   Monte Carlo inherits it, since each member rebuilds its wind by the same procedure.
3. **No shear clip** unless 1 and 2 still leave absurd equatorial winds; then it returns as a safety
   bound, labeled as one.
4. **The budget's size.** If the overall ratio stays well below 1 after the rerun, report it with the
   radio and IRIS ratios; no tuning on four profiles. Equal scaling of prior and errors would not move
   the best guess, only its uncertainty; Cassini data will set it.

**The two diagnoses (§4.1, §4.2): approved,** in parallel with the change.
- **The southern march:** find the stalled latitude and what switches at pass 38; fix the transfer or the
  hand-off of the 2.5 deg wind. No tolerance loosened: the tool must work for any wind.
- **IRIS 290 mbar in 1c:** 500 draws for the ratios at the final linearization, each group's effective
  degrees of freedom, and the condition number of S.

**Departures (§5): accepted.** `voyager_fit_1` refits at its own commit; the generalized
`build_matrix.py`; a resume point after the matrix runs.

**Order of work:** the diagnoses and the two rules; then Step 1 rerun, with 1d also reporting the
equatorial deformation radius used and where the stability constraint was active. No regression;
Lindal's files untouched.
