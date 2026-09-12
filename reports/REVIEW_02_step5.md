# REVIEW 02, Step 5. The standard diagnostics, `casspian.tools.plots`

Review of `reports/REPORT_02_step5.md` against SPEC_02 v0.7 Step 5 and SPEC_00 v0.14 §3.5
and §7.1. Reviewer: S. Rafkin, 12 September 2026. **Disposition: accepted with two small
figure changes**, to be made before the acceptance commit. SPEC_02 closes at this step and
reopens at v0.8 with a new Step 6 (equatorial anchoring); see the order of work.

## The figures

F1 to F4 were read against the numbers the earlier steps established, and they agree: at
`phi_c` the rigid centrifugal part of `g_N - g_eff` is 1.158 m/s² at 100 mbar (Ω² r cos²φ
with r = 58,520 km gives 1.158) and the wind part 6.1e-4 m/s² (2Ωu cosφ with u = 2.17 m/s
gives 6.1e-4); `G_phi` = -0.958 m/s² against `g_eff` = 9.96 gives ψ = 5.49°, which is the
frozen 5.495°; `g_N` falls by 0.15 m/s² over the 400 km from 1 bar to 0.2 mbar, which is
2Δr/r; F3 shows the 28.7 km asymmetry with the south pole higher, the 127 km equatorial
dynamical height, and 65 km at the anchor; F4's radius panel puts the anchor line at
58,519.883 km and the temperature recovery is at machine precision. Every profile panel has
the anchor isobar, every latitude panel the profile latitude, every figure the footer with the
commit and hash. The style is plain and consistent. The layout fixes described in decision 8
worked; nothing collides.

Two changes:

1. **F4, the refractivity panel.** The ±1σ band is 2.3 percent of N and is invisible on the
   log axis, so the one uncertainty the product states does not show. Draw the fractional
   uncertainty `refractivity_uncertainty / refractivity` as a second curve on a top axis of
   that panel (linear, in percent), labeled as such, keeping the band. When the band later
   becomes visible (a Monte Carlo product), both will still read correctly.
2. **F1, the wind inset.** The inset's horizontal axis carries no label. Label it
   `planetocentric latitude (deg)`, short form acceptable.

Both are code only in `tools/plots`; rerender, replace the two copies under
`reports/figures/`, and note the change in the report.

## Rulings on the decisions

1. **The acceptance runs on a copy of the manifest.** Correct for this step; changing the
   committed manifest would have changed a recorded hash mid-step. The question for the
   author is answered: yes, the Lindal manifest is to carry `[diagnostics] figures = true`,
   and it is done in Step 6, which rebuilds the manifest and the product for its own reason,
   so the hashes change once, not twice.
2. **Time band masking.** Sound. The same-timestamp byte identity is the stronger check and
   it is good that it was run.
3. **Names, and the combined PDF for kind N only.** Accepted.
4. **The wind view's `u(p)` at five latitudes.** Accepted; a file offered alone has no
   `phi_c`.
5. **`--format` and `--dpi` on the CLI.** Accepted.
6. **F5 and F6 on the named per-level fields; `Phi` on the anchor isobar and the boundary
   pressure not implemented.** Accepted. SPEC_03 will fix the field names; `geopotential_m2s2`,
   `pressure_hydrostatic_Pa` and `boundary_pressure_Pa` are recorded there as the placeholders
   this code already understands, so SPEC_03 either adopts them or renames them in both
   places.
7. **Wind interpolated linearly in log pressure for F2.** Accepted; a rendering choice that
   reaches no product, and for Lindal a no-op.
8. **Layout checked by eye.** This is the check the acceptance cannot make and the report
   was right to say so.

## Housekeeping

`lindal_refractivity.nc` on disk carries `4a9aca2...-dirty` because the Step 4 regression
suite drives `casspian-refrac`. That is the sweep rule working as intended: after the Step 5
acceptance commit the sweep rebuilds it clean, and the new hash supersedes `327f8141...` in
REPORT_02_step4 section 7. Step 6 will rebuild it again with the new manifest; record both.

## Order of work

1. The two figure changes; rerender; update `reports/figures/` and the report; commit; sweep;
   record the hash.
2. SPEC_02 v0.8 Step 6, equatorial anchoring, with the SPEC_01 v0.20 amendments to Steps 5
   and 9 and SPEC_00 v0.15 §7.1. Read the status lines and STATE.md first. Step 6 is a
   comparison step: its acceptance is that the comparison is reported completely and the
   author then decides the default rule.
