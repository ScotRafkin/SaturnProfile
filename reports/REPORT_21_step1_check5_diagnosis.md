# REPORT_21 Step 1, check 5: diagnosis of the Voyager 2 egress-referenced fit

**Status: reported, not accepted. Check 5 fails at 3.45 K after the swap fit converges (section 9). Under the instruction, SPEC_21 closes on hypothesis 4 with check 5 carried to SPEC_22. Nothing is committed.**

Log: `reports/step21_1/diag_swap.log`. Script: `reports/step21_1/diag_swap.py`. Both fits start from
`project(base_shear)` and are linearized with the production forward runs. The swap control is
`fits/voyager_fit_2/voyager_fit_2_reference_31s.toml`: reference 31.2 S (26.1 S), targets 36.3 N (30.8 N),
74.2 S (71.1 S) and 3.0 S (2.4 S). The control fit is `fits/voyager_fit_2/voyager_fit_2.toml`: reference 36.3 N
(30.8 N), targets 31.2 S (26.1 S), 74.2 S (71.1 S) and 3.0 S (2.4 S).

## 1. What prompted it

In the first resume (`reports/step21_1/run_resume_1.log`), the swap fit's solve asked for about 15 K every pass,
whatever the forward runs did:

| Pass | Solve's largest linear change (K) | Actual change in delivered T, 36.3 N / 74.2 S / 3.0 S (K) |
|---|---|---|
| 1 | 15.0 | 9.9 / 13.4 / 13.7 |
| 2 | 15.5 | 0.02 to 0.07 |
| 3 | 15.6 | 14.0 / 10.4 / 30.6 |
| 4 | 15.5 | 13.4 / 8.1 / 24.7 |

The Anderson weights were +1.004, then −135.36, then +1.633. In pass 5 the forward run reached 0.997 Pa, the
wind file's top, and crashed. Since then every wind file the fit writes extends one decade above 1 Pa
(REVIEW_21_step1 item 1).

## 2. Item 1: one-row steps at the swap fit's start

For each target, the step goes along its row nearest 110 mbar: `d = 0.01 B r / (r^T B r)`, a predicted change of
+0.0100 in ln T at that row (about 1 K). The forward runs are then rerun. The table gives the ratio of actual to
predicted change in ln T at every target. Each value is the median over rows with |predicted| > 1e-4, with the
correlation between actual and predicted in parentheses.

| Step along | 36.3 N (30.8 N) | 74.2 S (71.1 S) | 3.0 S (2.4 S) | Actual at the stepped row |
|---|---|---|---|---|
| 36.3 N, 108.7 mbar, \|d\| 7.57 | 1.037 (1.000) | 1.011 (1.000) | 1.012 (1.000) | +0.01042 |
| 74.2 S, 108.7 mbar, \|d\| 6.70 | 1.053 (1.000) | 1.014 (1.000) | 1.020 (1.000) | +0.01012 |
| 3.0 S, 108.7 mbar, \|d\| 22.1 | 1.058 (1.000) | 1.013 (1.000) | 0.973 (0.999) | +0.01021 |

The ratios over rows span about 0.85 to 1.48. The 3.0 S group in its own step reaches −2.1, at a row where the
predicted change crosses zero.

For comparison, the control fit, stepped along its 31.2 S (26.1 S) row at 114.7 mbar with \|d\| 8.03, gives
0.981, 0.969 and 0.978, with correlations of 0.999. The actual change at the stepped row was +0.00956.

**Finding.** At the start, the swap fit's operator describes its forward runs as well as the Lindal-referenced
one describes its own, to within 2 to 6 percent. The failure is not in the operator at the starting point.

## 3. Item 2: assumptions that the reference is Lindal's

**In the code: none remain.** The following are all general:

- The path integrals start at `anchor_latitude_deg`, taken from `[radio.reference]`.
- The anchor temperature behind the weight and the stability anchor are the reference's own.
- The shared anchor error uses `self.reference`, with its stated uncertainty and representativeness on the
  anchor's labels.
- Variable and file names are general. The earlier fix in `radio.temperature_on_wind_grid` reads the single
  `ln_refractivity_*` variable.

The swap control's target latitudes match its runs' namelists, and those runs use the Voyager 2 egress
refractivity as anchor. The orientation along each path is correct: the swap fit's starting misfit at 36.3 N is
close to the negative of the control fit's at 31.2 S (section 4).

**By the control file's choice.** The prior's temperature and deformation speed are Lindal's. The problem's wind,
composition, gravity and rotation, which `unit_weight` uses, come from base run `lindal_iris_v_sl_ii_110`. The
swap control sets both on purpose, so the prior is shared between the two fits.

**In the error budget: the floor.** This is one assumption tied to Lindal's anchor:

- Each delivered product sits on its anchor's levels. Between 1 Pa and 1.3 bar:

  | Profile | Levels | Median spacing in ln p |
  |---|---|---|
  | Lindal, 36.3 N (30.8 N) | 58 | 0.092 |
  | Voyager 2 egress, 31.2 S (26.1 S) | 307 | 0.023 |
  | Voyager 1 ingress, 74.2 S (71.1 S) | 300 | 0.023 |
  | Voyager 1 egress, 3.0 S (2.4 S) | 306 | 0.023 |

- The rows per target are therefore 390, 376 and 348 in the swap fit, against 65, 60 and 59 in the control.
- The 0.3 K floor (SPEC_19 v0.5 section 2 item 2) is added to each target's diagonal uncorrelated, row by row,
  in `radio.Radio.block`.
- The floor was sized on Lindal's levels to absorb the operator's rough departure of 0.1 to 0.3 K, which the
  Matérn matrix alone would have weighted as heavily as the whole radio noise.
- On levels 4 times as dense, an uncorrelated floor weighs 4 times as much per unit ln p, an effective floor near
  0.15 K. The departure it was meant to absorb then carries weight again.

I have not fixed this. Any general rule changes the error budget:

- Scaling the floor's variance by the row spacing needs a reference spacing.
- Delivering every target on a common grid would change the Lindal-referenced values.

The choice is for the review.

## 4. Item 3: starting misfits and the first solve

Starting misfit, observed minus delivered (K):

| Fit | Target | Rows | RMS | Largest | 730 | 290 | 110 | 30 | 10 | 1 mbar |
|---|---|---|---|---|---|---|---|---|---|---|
| swap | 36.3 N (30.8 N) | 390 | 5.84 | 15.29 | +0.84 | +2.11 | +5.05 | −0.25 | −1.77 | +5.25 |
| swap | 74.2 S (71.1 S) | 376 | 5.01 | 11.97 | −2.08 | +2.20 | +4.47 | +3.06 | +6.58 | +1.00 |
| swap | 3.0 S (2.4 S) | 348 | 7.30 | 13.85 | +0.32 | −0.04 | +4.13 | −5.12 | +8.82 | −11.29 |
| control | 31.2 S (26.1 S) | 65 | 4.29 | 15.01 | −0.56 | −2.01 | −4.88 | +0.41 | +1.49 | −5.22 |
| control | 74.2 S (71.1 S) | 60 | 3.84 | 11.53 | −2.69 | +0.15 | −0.38 | +3.52 | +8.09 | −4.21 |
| control | 3.0 S (2.4 S) | 59 | 7.52 | 19.82 | −0.28 | −2.16 | −0.75 | −4.84 | +10.32 | −16.45 |

The first robust solve, with Huber's rule:

| Fit | \|step\| | Rows downweighted | Largest predicted change per target (K) | RMS misfit after (K) |
|---|---|---|---|---|
| swap | 445 | 5 | 9.66 / 13.30 / 13.60 | 2.44 / 3.59 / 6.65 |
| control | 453 | 0 | 6.73 / 6.89 / 11.14 | 1.98 / 2.42 / 3.72 |

The largest predicted changes for IRIS are 0.178, 0.107 and 0.100 in the swap fit, against 0.178, 0.107 and
0.096 in the control. For stability they are 0.006 and 0.045.

**Why the first solve asks for about 15 K.** The misfits are that size, in both fits; the control fit's first
request is 11 K. The swap fit's request is larger at 74.2 S (13.3 against 6.9 K), yet its misfit after is larger.

The 36.3 N path, from 26.1 S to 30.8 N planetocentric, contains the 3.0 S path, from 26.1 S to 2.4 S. With 390
rows weighted under the tighter floor, 36.3 N pulls on the shared shear, and the misfit at 3.0 S stays at 6.65 K
RMS.

## 5. What is not yet separated

A first request of 15 K is expected. What failed is that the later passes kept asking for it. The operator holds
for small steps at the start (section 2), which leaves two candidates:

1. The full step (\|step\| 445, about 60 times the row steps) leaves the region where the operator holds.
2. The floor, tightened by the denser rows, makes the solve chase the operator's rough departure, so each
   relinearization's small changes move the solve far.

The control fit takes a step of the same size and converges, which points to the second. That is inferred, not
tested.

## 6. Proposed test (needs approval)

Linearize the swap fit at its start again (about 15 minutes), with the floor's variance scaled by each target's
row spacing over 0.092, the Lindal spacing. Compare the first solve with section 4, then run two passes and see
whether the request falls. The scaling is for the test only, not a ruling on the floor.

No Anderson safeguard has been added. The wind extension, the change to 1c and the addition to check 2 stand.
Check 5 will be rerun from the resume point after the review's ruling.

## 7. The correlated floor, tested

**The change (ruled after section 6).**

- The 0.3 K floor is now correlated in ln p over 0.115, the forward wind grid's spacing. It is an input,
  `floor_correlation_ln_p`, with its form `floor_correlation_form`, both in `[radio]` of both control files with
  their basis in a comment. With no length given it stays uncorrelated, so `voyager_fit_1` is unchanged.
- The form is exponential (Matérn 1/2), added to `prior.correlation`. A Matérn 3/2 floor would remove the diagonal
  SPEC_19 added for conditioning. On a test matrix of 300 levels 0.023 apart, with 1 K Matérn 3/2 over 0.5 beside
  the floor, the smallest mode's sd is:

  | Floor | 0.023 spacing | 0.092 spacing |
  |---|---|---|
  | uncorrelated (SPEC_19) | 0.300 K | 0.304 K |
  | exponential over 0.115 | 0.095 K | 0.192 K |
  | Matérn 3/2 over 0.115 | 0.019 K | 0.129 K |

**The test.** The swap fit from its start: two passes as `driver.iterate` makes them, then the solve on the third
linearization (`reports/step21_1/floor_test.log`, script `floor_test.py` beside it). The table gives, per target,
the solve's largest request, then the RMS misfit before and after the solve's linear change (K).

| Pass | \|step\| | Huber downweighted | 36.3 N (30.8 N) | 74.2 S (71.1 S) | 3.0 S (2.4 S) | Actual change after the pass, 36.3 N / 74.2 S / 3.0 S |
|---|---|---|---|---|---|---|
| 1 (plain) | 583 | 48 | 14.3 / 6.07 / 5.60 | 25.4 / 5.15 / 8.73 | 21.0 / 7.11 / 9.79 | 13.7 / 24.1 / 19.5 |
| 2 (Anderson +1.020) | 7.67e3 | 78 | 16.8 / 5.53 / 14.57 | 41.5 / 7.99 / 30.72 | 42.1 / 9.17 / 33.34 | 0.36 / 0.96 / 0.68 |
| 3 (solve only) | 7.5e3 | 76 | 16.9 / 5.37 / 14.68 | 42.9 / 7.52 / 30.87 | 43.3 / 8.87 / 33.65 | |

**Finding: the request does not fall; it grows.** Under the correlated floor the first request is 25 K, against
13.6 K under the uncorrelated one, and by pass 2 it is 42 K.

**What the passes show.**

- **Pass 1 is linear.** The actual change follows the request to within 5 percent (24.1 against 25.4 K).
- **The solve makes the radio misfit worse.** In every solve the RMS misfit after the solve's own linear change is
  larger than before at 74.2 S and 3.0 S, and from pass 2 at all three targets: 8 K before, 31 K after. A
  Gauss-Newton solve of this problem should not raise the misfit of the data that drive it.
  - So the solve is moving the state along directions its error budget makes cheap, not toward the radio rows.
  - The step grows thirteenfold from pass 1 to pass 2, and Huber's rule downweights 48, then 78 rows.
- **The Anderson weight near +1 explains the near-still passes.** In this run's pass 2 and the first resume's pass
  2, the weight was near +1 (+1.020 and +1.004). That mixes the new solve almost entirely back to the previous one,
  so the state barely moves (0.4 to 1.0 K here, 0.02 to 0.07 K before). It is a symptom, not the cause.
- **What the correlated floor does to the problem.** It removes most of the uncorrelated diagonal on 300 closely
  spaced rows: the smallest mode's sd falls from 0.30 to 0.095 K. Large, smooth departures then grow cheaper in R
  relative to rough ones. Under the swap, the shared anchor term and the representativeness, up to 15 K aloft and
  correlated across all three targets, sit on the Voyager 2 egress labels, also 300 rows. In the Lindal-referenced
  fit they sit on 58.

**Not done, per the instruction.** No step limit or Anderson safeguard has been added. Hypothesis 4 is not refit and
check 5 is not rerun, since the condition for both was that the request fall.

## 8. Proposed next diagnosis (needs approval)

The failure is now in the solve, not the operator: the robust solve raises the misfit of every radio target. I
propose to save the swap fit's pass 1 and pass 2 linearizations (about 30 minutes), then, with no forward runs:

1. Break the pass 2 solve's cost into its parts: the prior, each radio target, the shared anchor term, IRIS and
   stability, before and after. This shows what the 7.5e3 step buys.
2. Solve the same linearization without Huber's rule, and with the floor uncorrelated, to see which one turns the
   solve against the radio rows.
3. Find R's smallest modes on the swap's rows, and the solve's projection on them.

The `--refit` mode for the acceptance script (hypothesis 4 refit, its change at the six levels, then check 5) is
written and ready for when the swap converges.

## 9. Observations at the forward model's resolution

Ruled after section 8: the fit had been asked to use structure below the forward model's resolution.

### 9.1 The change

- **Layer averaging.** With `average_layers = true` in `[radio]`, input and basis in both control files, every
  radio profile is averaged in ln p into the layers between the forward wind grid's nodes. These are the 121 levels,
  0.115 apart, that the stability constraint uses. The averaging is the exact mean of the piecewise-linear profile
  over each layer, `radio.layer_mean`.
  - The observed profile is averaged on its own levels.
  - The delivered profile, its rows `H` and every error part are averaged with the same matrix: the target's own
    error and representativeness, and the reference's shared term on the anchor's labels.
  - One row per layer, only for layers both profiles span.
- **The floor.** It is uncorrelated again, added once per layer. The correlated form stays in code
  (`floor_correlation_ln_p`), unused.
- **The cost breakdown.** `solver.cost_parts` breaks each solve's cost into the prior, the shared anchor term,
  each target, each IRIS level and stability. `driver.iterate` prints it every pass.
  - Before the solve, the cost is taken at the linearization's state; after it, at the solve's linear prediction.
  - The radio part `v^T R v` (`v = R^-1 r`) splits exactly into the shared term and each target's own block.
- **Unchanged.** Without `average_layers`, `voyager_fit_1` is unchanged. No step limit or Anderson safeguard was
  added.
- **Item 3, the trust region.** It was not needed, since the request fell. I wrote it before the test finished and
  took it out after, so the code does not carry it.

### 9.2 The test: swap fit, three passes from its start

Log: `reports/step21_1/layers_test.log`.

| Pass | Solve's request (K) | Actual change, 36.3 N / 74.2 S / 3.0 S (K) | Cost before -> after the solve |
|---|---|---|---|
| 1 (plain) | 9.07 | 6.88 / 6.04 / 8.25 | 208.7 -> 34.2 |
| 2 (Anderson +0.029) | 0.42 | 0.41 / 0.25 / 0.26 | 37.2 -> 35.1 |
| 3 (Anderson −0.029) | 0.38 | 0.39 / 0.16 / 0.15 | 34.9 -> 34.7 |

The request falls from 9 K to 0.4 K in one pass. Under the rows of section 4 it stayed near 15 K, and under the
correlated floor it grew to 42 K. In every solve the cost falls at every target as well as in total. The full swap
fit, in the check 5 run below, converged in 4 and 2 passes.

### 9.3 The cost breakdown of each solve

Log: `reports/step21_1/run_refit.log`.

**Hypothesis 4, the Lindal-referenced fit, first fit:**

| Pass | Prior | Shared anchor | 31.2 S | 74.2 S | 3.0 S | IRIS 110 | IRIS 290 | IRIS 730 | Stability | Total |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 110.1 -> 15.8 | 3.54 -> 1.69 | 15.6 -> 13.7 | 13.8 -> 13.2 | 22.7 -> 18.1 | 48.9 -> 2.66 | 13.0 -> 1.20 | 10.8 -> 0.86 | 205.1 -> 2.27 | 443.5 -> 69.5 |
| 2 | 15.8 -> 14.9 | 1.71 -> 1.71 | 13.8 -> 13.7 | 13.2 -> 13.2 | 21.5 -> 21.4 | 4.72 -> 2.82 | 1.28 -> 1.18 | 0.96 -> 1.29 | 1.37 -> 0.07 | 74.3 -> 70.3 |
| 4 | 14.9 -> 14.9 | 1.65 -> 1.64 | 13.7 -> 13.7 | 13.2 -> 13.2 | 21.4 -> 21.4 | 3.00 -> 3.00 | 1.24 -> 1.24 | 1.30 -> 1.28 | 0.13 -> 0.17 | 70.5 -> 70.5 |

The second fit, after the IRIS 110 mbar part of 0.0514 per radian, ends at 69.5 in 3 passes.

**The swap fit, 31.2 S-referenced, first fit:**

| Pass | Prior | Shared anchor | 36.3 N | 74.2 S | 3.0 S | IRIS 110 | IRIS 290 | IRIS 730 | Stability | Total |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 110.1 -> 15.3 | 4.50 -> 2.51 | 3.75 -> 1.80 | 4.12 -> 3.49 | 10.7 -> 6.51 | 48.9 -> 2.66 | 13.0 -> 1.08 | 10.8 -> 0.72 | 2.95 -> 0.07 | 208.7 -> 34.2 |
| 4 | 14.5 -> 14.5 | 2.57 -> 2.57 | 1.90 -> 1.90 | 3.54 -> 3.54 | 6.99 -> 6.98 | 2.73 -> 2.70 | 1.23 -> 1.22 | 0.98 -> 0.98 | 0.09 -> 0.10 | 34.5 -> 34.5 |

The second fit, after the IRIS 110 mbar part of 0.0458 per radian, ends at 33.8 in 2 passes.

What the breakdown shows:

- **The prior and IRIS buy most of the first solve.** The prior part falls from 110 to 15 and IRIS from 73 to 5. The
  radio terms change little.
- **The two fits end at different costs.** Each radio target's term is 4 to 7 times larger in the Lindal-referenced
  fit than in the swap fit at the same targets: 74.2 S 13.2 against 3.5, 3.0 S 21.4 against 7.0. Its stability
  term at the start is 205, against 3.
- **The radio data fit easily within their errors in both.** The misfit over expected after IRIS's parts is 0.20,
  0.19 and 0.33 at the Lindal-referenced fit's targets, and 0.04, 0.06 and 0.12 at the swap fit's.

### 9.4 Hypothesis 4 refit (`--refit`)

**Passes: under the 1 K bound.** Delivered T, the refit minus the previous fit (K):

| Target | 730 | 290 | 110 | 30 | 10 | 1 mbar | Largest, 1 mbar to 1.3 bar |
|---|---|---|---|---|---|---|---|
| 31.2 S (26.1 S) | −0.04 | −0.03 | +0.01 | −0.05 | +0.11 | +0.08 | 0.12 |
| 74.2 S (71.1 S) | −0.18 | −0.15 | −0.05 | −0.09 | +0.02 | +0.10 | 0.23 |
| 3.0 S (2.4 S) | −0.02 | −0.01 | +0.02 | −0.01 | +0.33 | +0.75 | 0.75 |

The refit's wind replaces `fits/voyager_fit_2/voyager_fit_2_wind.nc`. The previous file is kept outside the
repository.

### 9.5 Check 5

**Fails: largest difference 3.45 K, bound 0.5 K.** The swap fit converged (4 and 2 passes). Delivered T, the
31.2 S-referenced wind minus the Lindal-referenced wind, through the same runs (K):

| Profile, from anchor | Largest | RMS | 730 | 290 | 110 | 30 | 10 | 1 mbar |
|---|---|---|---|---|---|---|---|---|
| 31.2 S (26.1 S), from Lindal | 1.90 | 0.48 | −0.49 | −0.24 | +0.13 | −0.16 | −0.05 | −0.29 |
| 74.2 S (71.1 S), from Lindal | 0.79 | 0.21 | −0.19 | −0.02 | +0.08 | −0.14 | −0.08 | −0.06 |
| 3.0 S (2.4 S), from Lindal | 3.45 | 0.87 | −0.25 | −0.10 | +0.13 | −0.21 | −0.68 | −0.40 |
| 36.3 N (30.8 N), from 31.2 S | 1.93 | 0.44 | +0.49 | +0.27 | −0.13 | +0.15 | −0.02 | +0.27 |

- **At the six levels.** Every difference is within 0.5 K except 3.0 S at 10 mbar (−0.68 K). The largest values lie
  elsewhere in the profiles.
- **Where the largest lie.** Check 5 deletes its runs after comparing them (`work_reference`), so I can't place
  them without rerunning.
- **What check 5 compares.** It compares delivered profiles at full resolution, while the fit now constrains them
  only as layer means.
- **The sign pattern.** The two directions of the 36.3 N and 31.2 S pair are close to mirror images (+0.49 against
  −0.49 at 730 mbar). That is what a difference in the wind between the two latitudes, rather than in either
  profile, would give.

Under the instruction, check 5 failing after the swap converges means SPEC_21 closes on hypothesis 4, with check 5
carried to SPEC_22. Nothing is committed. Since nothing is committed, the following are left for the closure's
ruling:

- the layer averaging and cost breakdown in the code;
- the refit wind file;
- `--refit` in `tests/step21_1/`;
- REPORT_21_step1's section 7, which still misreports check 5.

Two questions for SPEC_22, neither tested:

- whether check 5 should be judged on layer means, the resolution the fit now uses;
- why the Lindal-referenced fit's radio terms are 4 to 7 times the swap fit's at the same targets.
