# REVIEW 18, Step 1. Kliore et al. (1980) as a second registered profile

Reviewing agent, 5 November 2026. Report: `reports/REPORT_18_step1.md`. Specification:
`docs/specs/SPEC_18_Kliore_Kind_N.md` v0.4.

**Accepted. SPEC_18 closes.**

- **Lindal's answer is unchanged.** No registered file changed. His kinds C, T and D and the kind C of
  all 23 runs, rebuilt from scratch in a copy by the split tools, differ in no value, and the
  manifest is byte-equal. No other suite is rerun, as SPEC_18 §0 directs.
- **Kliore is a working kind N.** Refractivity to RMS 0.031 n-units against the table's 0.05 half
  step; radius to +1.27 km with no spread, so the height ruling keeps his spacing exactly; the label
  returns 9.5 S; kinds D, G, R and W equal Lindal's; none of Kliore's files records a Lindal file
  among its inputs.
- **The closure's drift is the source's, explained.**
  - The table integrated outside the model reproduces the closure to 8.4e-4: the model is not
    involved.
  - From 10 mbar down the residual is a constant, +6.0e-3, plus an offset of about -10 Pa from the
    top that fades with pressure. With both removed, 108 of 112 levels are inside the printed
    precision.
  - The constant is the product of molar mass and gravity, 0.6 percent below ours. The table cannot
    split the two. The molar-mass reading (2.2 amu, the 90/10 mixture with whole-number masses)
    fits with no fitted value and brings 110 of 112 levels inside. It is likely, and recorded as an
    inference: neither Kliore paper states the number.
  - This is enough to close. It is a property of Kliore's printed pressures and temperatures, not
    of his measurement or of our model.
- **For when Kliore becomes an anchor.** His measured quantity is refractivity against radius;
  his temperatures and pressures carry his constants. That ruling should consider re-deriving them
  from his refractivity with the exact molar mass, as the measurement, rather than inheriting the
  0.6 percent. It belongs with ruling 3 (his uncertainties). The inference of 2.2 amu goes in a
  comment in `kliore1980_scalars.toml` now, so it is not lost.

**On the report's three points:**
1. **The closure run is not registered.** Its inputs and product rebuild from its two control
   files, so by the repository rule only those are committed. No `.gitignore` change.
2. **Kliore's products are rebuilt clean** after the acceptance commit, as earlier steps did, and
   committed with the hash commit.
3. **`pip install -e .`** before the acceptance commit, so the two new commands work as the runbook
   says.

**Order of work:** commit and push, with no regression.
1. The acceptance commit: the tools, `occul_data/kliore1980/`, the closure's two control files,
   `pyproject.toml`, `docs/RUNBOOK.md`, `tests/step18_1/`, both SPEC_18 reports and this review,
   SPEC_18 v0.4, `STATE.md`, the status note, and
   `docs/CASSPIAN_Wind_Estimation_Design_Note_2026-11-05.md` (accepted at v0.4).
2. Kliore's products rebuilt clean at that commit; with them, the hash in SPEC_18's row of
   `STATE.md`. Push.
