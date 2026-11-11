# REPORT 21, pre-execution: the general error budget, and hypothesis 4

SPEC_21 v0.2, Step 0, before any code changes. 9 November 2026. For the author and the reviewing agent.
`main` is at `a24b2f1`, pushed. No code has changed. The check of §1 is a read-only scratch computation on
the base wind of `voyager_fit_1`.

**In short.**
- **The suggested construction works, with one addition:** condition the gradient to vanish at the poles
  as well as the equator.
- **The prior's draws are finite and smooth at the equator,** as required (§1).
- **With §2's first values the equatorial prior is very loose:** about 600 m/s per unit `ln p` at 110 mbar.
  Toward the poles it is tighter than SPEC_19's. Not a block, since §2 is to be revisited on results, but
  the author should see it before the fit (§4 item 1).
- **The 2.5 deg basis costs almost nothing.** A fit's time is its forward runs, as before (§2).
- **Lindal's uncertainty needs a file of its own** beside his registered ones (§3).
- **Nothing blocks.** §4 lists eight choices for a ruling.

## 1. The prior: from temperature to shear

**The construction.** SPEC_17's linearized balance on an isobar is `d ln T / dφ = w s`, with

    w = 2 Ω_abs r m_bar (sin φ - (cos φ / r)(dr/dφ)_p) / (R T)        (lindal_wind.balance_weight)

So `w` has a simple zero at the equator, where the bracket changes sign, and is finite elsewhere.

1. **The input** is a temperature anomaly on isobars, `T'(φ, ln p)`. Its amplitude is `σ_T(p)`, and its
   correlation is separable: Matérn 5/2 in latitude (planetocentric, radians) and the ruled form in `ln p`.
2. **The gradient** `g = d ln T / dφ = (1 / T̄) dT'/dφ` is then a stationary Gaussian process. Its covariance
   in latitude is `-ρ''`, and for Matérn 5/2 with length `ℓ` and `q = √5 |Δφ| / ℓ`:
   `(5 / 3ℓ²)(1 + q - q²) e^(-q)`.
   - This is the spec's "stationary prior on `d ln T / dφ`", derived from the temperature covariance
     rather than stated beside it.
   - Matérn 5/2 is the least smooth form whose gradient is itself differentiable, which the limit in
     item 4 needs.
3. **Conditioned on zero gradient at the equator** (`K - k(φ,0) k(0,0)^-1 k(0,φ')`), as the balance demands
   for any finite shear.
   - **Added: also at both poles.** The basis drops the polar functions, so `s = 0` there. With `w` finite
     at the poles, `g` must vanish there too. A field smooth on the sphere has zero gradient at a pole, so
     this is physics and not a choice. It changes the shear sd by 5 percent at 87.5 deg and by nothing
     equatorward of 70 deg.
4. **Divided by `w`:** `s = g / w`, so `B_s(φ, φ') = K_g(φ, φ') / (w(φ) w(φ'))`, scaled by
   `σ_T(p) σ_T(p') / (T̄ T̄')` and the vertical correlation.
   - **At the equator** both vanish, and the limit is `s(0) = g'(0) / w'(0)`. Its variance comes from the
     conditioned kernel's mixed second derivative, so it is finite whenever `T'` has a curvature.
   - **In code** the equator node takes this limit explicitly. Dividing there would not work, because the
     computed `w` at the node is `-1e-8`, not 0.
5. **Onto the basis:** `B_c = P B_s Pᵀ`, with `P` the basis's least-squares projection from a fine grid,
   0.5 deg and 0.1 in `ln p`. `B_c` is dense: the `w` scaling breaks the separable form `prior.apply`
   uses.
6. **What `w` is evaluated with:**
   - **Once**, with the base run's wind (Sanchez-Lavega at 398 mbar with the IRIS shear) and the anchor's
     temperature `T̄(p)`, uniform in latitude.
   - **Held for the whole fit**, since the prior must not move with the iterate.
   - **The uniform temperature** is a first-order choice: `T` varies by about 10 percent along an isobar,
     and enters `w` as `1/T`.

**The check, at 110 mbar** (σ_T 6.0 K; Matérn 5/2, 5 deg; conditioned at the equator and poles). Shear sd in
m/s per unit `ln p`, planetocentric:

| Latitude | 0 (limit) | 0.5 | 1 | 2 | 5 | 10 | 20 | 40 | 60 | 74 | 87.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 110 mbar | **600** | 535 | 467 | 345 | 152 | 77 | 43 | 25 | 21 | 19.5 | 18.5 |
| 1000 mbar | 197 | 176 | 154 | 114 | 50 | 26 | 14 | 8 | 7 | 6.5 | 6.1 |
| 10 mbar | 1508 | 1346 | 1175 | 869 | 383 | 192 | 107 | 63 | 52 | 48 | 46 |

- **Finite at the equator:** the sd rises smoothly to its limit.
- **Smooth at the equator:** in 20 draws at 110 mbar, the equator node differs from the mean of its two
  neighbors by at most 0.07. Draws change by 130 to 650 between neighboring 0.5 deg nodes.
- **Every draw is finite** everywhere.
- **What the draws carry:** a temperature span along the isobar of 19 to 33 K, about ±2.5 σ_T over
  36 correlation lengths, as a Gaussian field should.

**What this means for the fit.**
- **Poleward of 20 deg** the prior is tighter than SPEC_19's uniform 40: 19.5 at 74 deg. That is what was
  wanted against the unphysical draws at 74.2 S.
- **Within 5 deg of the equator** it is far looser than SPEC_19's: 600 at 110 mbar against hypothesis 1's
  jet shear of about 60.
  - **The cause:** σ_T 6 K over a 5 deg length means curvatures of `σ_T / ℓ²`, and the balance turns
    curvature into shear divided by a small `w'(0)`.
  - **This is the physics the design note asks for,** at a first value that may be generous. The radio
    paths through the equator and the robust weighting will carry the fit there.
  - **Acceptance 1(d) checks only temperature draws.** It would not show a wind draw of 1000 m/s per unit
    `ln p`, so §4 item 1 proposes adding the wind.

## 2. The cost

| | 5 deg knots (SPEC_19) | 2.5 deg knots |
|---|---|---|
| Latitude functions (polar dropped) | 37 | 73 |
| `ln p` functions (0.5, half of one scale height) | 21 | 21 |
| Unknowns | 777 | 1533 |
| Dense `B` and its Cholesky | 0.06 s | 0.10 s |

- **The basis is cheap.** The solver works in observation space (289 rows), so its cost hardly depends on
  the unknowns.
- **The rest scales mildly:**
  - The posterior wind sd at every level is about 50 GFLOP, seconds.
  - The radio operator's path integrals and IRIS's balance rows scale with the latitude functions, small
    beside the runs.
  - Building `B_c` once from the fine grid (361 × 60 points, a few seconds).
- **A pass** is three forward runs and the balance weights, about 5 minutes alone, unchanged. The SPEC_19 fit
  with the floor took 5 passes.
- **The second fit for IRIS's systematic part** starts from the first fit's state with the term added. I
  expect 2 to 4 passes, so about 10 to 20 minutes.
- **One full fit is therefore about 40 to 50 minutes.**
- **Step 1's acceptance** needs about seven fits:
  - 1(a) and 1(b) nonlinear;
  - the real fit and its second fit;
  - the start from the prior mean;
  - the reference swap of acceptance 5, with its own second fit.

  With the four sections and the runs, the whole acceptance is about 4 to 5 hours.
- **The equatorial resolution.** The knots resolve the 5 deg correlation length, half of it as ruled. But the
  prior's sd falls by half within 2 deg of the equator, so 2.5 deg knots will smooth that peak somewhat. I
  will report the projected sd against the exact one there.

## 3. The ingest of each profile's uncertainty

**The Voyager profiles** (source-specific layer, `tools/lindal/build_profile_raw.py`, rebuilt):
- **Where:** each profile's scalars transcription gains an `[uncertainty]` table. Stage two writes kind T's
  `temperature_uncertainty_K`, NaN since SPEC_18, with the method in its attributes.
- **Voyager 2 egress:** its own traced Fig. 4 sigma_T, 390 of 391 rows. The deepest row, 1337.7 mbar, is
  beyond the traced bars and stays NaN unless a value is ruled (§4 item 6).
- **Voyager 1 ingress and egress:** the Fig. 4 curve at their pressures, times 2 (ruled), plus 2 K in
  quadrature from 600 to 1200 mbar at 74.2 S. The factor, the band and the source sentence live in the
  transcription, not in code.
- **What leaves the control file:** the hand-set 15 K at 3.0 S and 2 K at 74.2 S. Representativeness is not
  written at ingest; the fit adds it, by observation type (§1 item 4).
- **Rebuild:** the three anchors' kind T, D and N files rebuild, and Lindal's are not touched. The kind N
  refractivity uncertainty picks up the temperature uncertainty through `refrac`'s companions, as it would
  for any source.

**Lindal's anchor.** His registered `lindal_thermo.nc` carries NaN and must not change. I propose:
- **A new file:** `occul_data/lindal/lindal_temperature_uncertainty.csv`, the Fig. 4 curve on his 66 levels,
  written by the same ingest module from his kind T's pressures and the egress trace's sigma column. Its
  method goes in its header.
- **How the fit reads it:** each profile's errors from its kind T, or from a file named beside the profile in
  the control file when its kind T has none.
- **Why this is not a rule for one dataset:** only Lindal's entry names such a file, and the fit holds no
  value from it.
- **How the fit reads observations:** the observed profiles move from the traced CSVs to each profile's kind
  T, so the fit sees what the reduction saw. For the Voyager profiles they are the same numbers.

## 4. For a ruling

None of these blocks the step.

1. **The equatorial looseness (§1).** I propose:
   - keeping §2's values, as ruled;
   - adding the wind to acceptance 1(d), as the largest shear and the wind change over 10 to 1000 mbar of
     the 20 draws at the equator;
   - and, if the author prefers, a separate `σ_T` amplitude within 10 deg of the equator. That is the same
     rule with a different input.
2. **Conditioning at the poles** as well as the equator (§1 item 3). I propose it.
3. **`w` held at the base wind and the anchor's `T̄(p)`** (§1 item 6). I propose it.
4. **The correlation form of `σ_rep`.** §2 gives its amplitude, not its correlation. I propose Matérn 3/2 over
   0.5 scale height, the layering's own scale and the radio's form, as an input. The alternative is
   uncorrelated, which would weigh `σ_rep` more lightly.
5. **Huber with correlated errors.** A residual is "normalized" against its own sd or against the whitened
   residual. I propose the whitened residual within each group: the Cholesky factor of that group's `R`,
   with the threshold, 2, applied per whitened row and the weights carried back.
6. **The deepest Voyager 2 egress row** (1337.7 mbar, no traced bar): the nearest bar's value, or the row
   dropped. I propose the nearest bar.
7. **IRIS rows.** "No exclusions" (§1 item 9) drops the 5 deg equatorial omission, and I propose that.
   - **The 2 deg gap rule I propose to keep:** it is not an exclusion but absence of data. Across −27.7 to
     −38.5 deg the IRIS fit only bridges, so a row there would be an interpolation presented as a
     measurement.
   - **Near the equator** IRIS's rows come in with their operator's own weight, which `w` drives toward 0.
8. **IRIS's uncertain level** (§1 item 6): its operator at 110 and at 150 mbar, recomputed each pass, since the
   operator depends on the state. The difference enters `R` as an uncorrelated term on each 110 mbar row. A
   correlated alternative is possible if the author prefers.

No regression is due. The model does not change.
