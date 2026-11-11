# REVIEW 21, Step 1, preliminary under v0.4. The stability constraint's form

Reviewing agent, 10 November 2026. Report: `reports/REPORT_21_step1_preliminary_v04.md`. Specification:
`docs/specs/SPEC_21_General_Error_Budget.md` v0.4, amended to v0.5.

**Not accepted; Step 1 reruns under v0.5.** The diagnosis is right and clearly shown.

- **The deformation-radius floor works:** with the stability constraint off the fit converges in 3
  passes; the equatorial length is 7 to 9.3 deg in the stratosphere.
- **The two diagnoses are closed.** The southern march was duplicate latitude nodes in the section mesh,
  fixed in `lib.mesh.build_mesh` for every mesh; no forward run had such a pair. IRIS 290 mbar in 1c was
  sampling: about 5 effective degrees of freedom, 0.98 over 500 draws. Both accepted, with the added
  lines in 1c.
- **The stability constraint fails on its form, not its physics.**
  1. Its pairs are delivered levels 0.047 apart in `ln p`, inside one interval of the forward run's wind
     grid (0.115), where the forward run holds the wind linear and moves both levels together. The rows
     predict a slope change about eight times what the forward run makes, so the iteration cannot settle.
  2. The pair it pins (1045 to 1095 mbar) is slightly superadiabatic in the anchor's own profile, within
     its error. No wind can remove that, and none should: it is part of the observation.

**Rulings (author, 10 November):**
1. **Stability at the forward model's resolution.** Layer slopes between the forward run's wind-grid
   nodes, with the rows formed from the wind as the forward run sees it (linear in `ln p` between nodes).
2. **Relative to the source.** A delivered profile may not be less stable than its source allows: the
   constraint is on the excess over the dry adiabat (normal hydrogen, `c_p` from composition) beyond the
   anchor's own excess in the same layer. Excess already in the source is reported, not penalized.
   General: no profile or level is named.
3. **A smooth one-sided penalty** (a softened hinge, its width an input with its basis), so layers do not
   switch the iteration on and off.
4. **No finer wind grid.** It would change every forward run and need a regression.
5. **Fallback, only if 1 to 3 still do not converge:** stability as a reported check on the fit, the
   penalty kept for the Monte Carlo draws. Report before using it.

**Order of work:** the three changes; Step 1 from its start (the synthetic fits refit too). Report where
the constraint is active and by how much. No regression; Lindal's files untouched.
