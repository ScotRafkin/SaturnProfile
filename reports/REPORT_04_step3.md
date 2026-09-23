# REPORT 04, Step 3. The kernels

Coding agent, 23 September 2026; refreshed under REVIEW_04_step3 and SPEC_04 v0.12.
Specification: `docs/specs/SPEC_04_Transfer.md` v0.12, section 5, decision Q and section 14.
Acceptance script `reports/step04_3/accept_step04_3.py`, output `reports/step04_3/output.txt`.
Regression `reports/step04_3/run_regression.sh`, output `reports/step04_3/regression.txt`.

**Nine of nine checks pass.**

The first filing was 7 of 9 and its headline conclusion was wrong. It reported that decision L and
Eq. A15 are in tension at order one and that the cylinder identity is lost at the 7.6e-2 level.
The mechanism was right and the measure was not: the transfer applies the **line integral** of the
kernel along an isobar, not its maximum, and the residual of a gridded wind alternates in sign
inside the cells holding a slope jump. Bounding `d ln N` by the maximum times the span overstates
a sign-alternating integrand by two orders of magnitude. Measured as the line integral, the
largest `d ln N` over the span is **4.746e-04** against a bound of 2e-3, and the isobar shift
**130.3 m2/s2** against 300. Decision Q states that as the truncation floor of a wind hypothesis on
a grid; neither decision L nor decision P changes.

Two of the three failures were mine and one was the specification's:

* **The slope discontinuity was my clamped inverse map, not Step 2 decision 8.** `Phi_of_p`
  continued at the last interval's slope as decision 8 states, but `p_of_Phi` used `np.interp`,
  which clamps. The mesh reaches one geopotential spacing past the anchor's levels at each end, so
  the clamped inverse put several nodes on one pressure: measured, **2 duplicate pressures** and a
  non-monotonic grid, where the continued inverse gives none. That single fault produced the
  eightfold degradation at the ends, the NaN of the first filing's finding 4, and most of check
  9's maximum. Corrected; check 4 now measures the whole mesh with no node set held aside.
* **Check 5 asked for round-off and measures truncation**, which v0.12 restates as `|S/g|` below
  1e-4 per radian. Measured **3.180e-05**.

Findings 1 to 5 carry the measurements and their rulings.

## 0. Before this step

HEAD `8aef692`, Step 2 accepted at `562c769`. Every product on disk is the swept one: the reduction
chain clean at `1c9b310` and the closure product at `780a2de`. This step registers nothing, changes
no existing module, and writes nothing inside the repository outside `reports/step04_3/`, which git
ignores. `lib/kernel.py` is new.

## 1. What was built

**`lib.kernel.shear_kernel`.** `S = 2 Omega_abs r (du/dZ)_R` of Eq. A15, with `Omega_abs` from `u`
at the node by Eq. A2. `u` is never differentiated on the mesh: the wind file is the only thing
that knows how `u` varies, and `lib.windfield` returns its interpolant's own partials in the file's
coordinates. Those are at fixed pressure and A15 needs them at fixed radius, so they are converted
by `isobar_slopes`.

**`lib.kernel.isobar_slopes`.** The change of variables, and the only place a mesh difference is
taken. The mesh is in `(phi, Phi)` and both slopes are wanted in `(phi, r)`, so with
`dr/dPhi = 1 / g` from Step 2:

    (dln p/dr)_phi  = g (dln p/dPhi)_phi
    (dln p/dphi)_r  = (dln p/dphi)_Phi - g (dr/dphi)_Phi (dln p/dPhi)_phi

**`lib.kernel.shear_integral`.** `I` of Eq. A16, the trapezoid of `S / g` in `Phi` taken outward
from the `Phi = 0` node in both directions, as `lib.geopotential` takes the geopotential, so `I` is
exactly zero there and no node near the gauge carries the round-off of a long cumulative sum.

**`lib.kernel.composition_term`** and **`_three_point_derivative`.** The second term of Eq. A27 on
the composition file's own latitude grid at each isobar label, by the three-point formula for
unequal spacing written in difference form. **`lib.kernel.transfer_kernel`.** `K` at a point of a
curve, `S / g` bilinear on the mesh and the composition term linear in latitude.

## 2. Decisions

1. **`shear_kernel` takes the mesh, the node arrays and the field, not `mesh` alone.** The
   deliverable writes `shear_kernel(mesh)`, but `S` needs the radius, the pressure, `g` and `u` at
   every node, and those are the columns' and the isobar map's, not the mesh's. The mesh carries
   the nodes and the difference rule; the caller carries the state on them.
2. **`shear_kernel` returns `(S, du_dZ)`.** The axial derivative is what the synthetic states have
   a closed form for, and returning it lets them check the quantity the identity is about rather
   than inferring it back through `2 Omega_abs r`.
3. **`isobar_slopes` is its own function.** It is the whole of the conversion from the file's
   coordinates to fixed radius, it is the only place a mesh difference enters the kernel, and the
   acceptance measures it on its own.
4. **The composition term uses a difference-form three-point derivative, not `np.gradient`.**
   Forced by finding 3.
5. **`transfer_kernel` takes the composition arguments optionally.** Omitting them gives the shear
   term alone, which is what a run whose composition is uniform in latitude computes anyway, and
   saves the caller assembling a zero array to add.
6. **The synthetic wind states are posed on the mesh's own columns and held in memory.** Forced by
   finding 4. They are never written: a field that is nonzero at the poles and does not span them
   is not a kind W file.

## 3. Acceptance results

| # | Check | Result |
|---|---|---|
| 1 | closure wind: `S/g` is the closed form, no difference taken | PASS, 4.163e-17 |
| 2 | largest `\|S/g\|` 0.0977 to 10 N, 0.0540 to 60 N | PASS, 0.097565 and 0.053882 |
| 3 | `I = 0` at the `Phi = 0` node of every column exactly | PASS, exactly 0.0 |
| 4 | synthetic `u = beta r sin(phi)`: `(du/dZ)_R = beta` | PASS, whole mesh |
| 5 | solid body: `\|S/g\|` below 1e-4 per radian (v0.12) | PASS, 3.180e-05 |
| 6 | composition term against its closed form to 1e-8 | PASS, 1.922e-10 |
| 7 | the run's composition gives exactly zero | PASS, exactly 0.0 |
| 8 | beyond: `transfer_kernel` at a node and off it | PASS, exactly 0 at a node |
| 9 | cylinder wind: line integrals below 2e-3 and 300 m2/s2 (v0.12) | PASS, 4.746e-04 and 130.3 |

**The closure wind (checks 1, 2, 3).** `du/dln p` is `0.0e+00` at every one of the 329,643 nodes,
so both conversion terms multiply zero and `S` reduces to `2 Omega_abs cos(phi) u'(phi)` with no
difference taken anywhere. The measured departure from that closed form is 4.163e-17, which is the
`r (cos phi / r)` round trip and not an approximation. The largest `|S/g|` is 0.097565 per radian
from `phi_c` to 10 N against the specification's 0.0977, and 0.053882 to 60 N against 0.0540. `I`
is exactly zero at the gauge node of all 421 columns.

**The cylinder-extended wind (check 9).** The fixed point converges in five passes. On the file's
own 0.5 degree grid the line integral of `S/g` along the anchor's isobars gives a largest
`d ln N` of **4.746e-04** against a bound of 2e-3, and an isobar shift of **130.3 m2/s2** against
300. The largest `|S/g|` is 5.408e-02 per radian, reported and not bounded, against the reviewing
agent's 0.050. Per level, against the v0.12 instance values: +1.30e-05 at the top against +1.2e-4,
-9.33e-05 at 10 mbar against -5.2e-4, **-1.78e-04 at the gauge against -1.7e-4**, and +2.21e-04 at
the bottom against -1.1e-3. Finding 5 carries what that agreement and disagreement means.

Resampling the cylinder construction at 0.5, 0.25 and 0.1 degrees with the mesh unchanged gives
4.746e-04, 3.415e-04 and 3.867e-04, ratios 1.39 and 0.88. No order is claimed, as v0.12 directs,
and the sequence is not monotonic.

**The composition (checks 6, 7, 8).** The synthetic field `x_He = 0.06 + 0.02 sin(phi)` with H2 as
the closure reproduces the closed form to 1.922e-10 against a bound of 1e-8, on a 0.01 degree grid.
The run's own composition, whose columns are identical, gives **exactly** zero. `transfer_kernel`
returns a node's own `S / g` exactly and adds the composition term off it.

## 4. Findings

**1. The solid-body check asked for round-off and measures truncation. Ruled and restated.**
`u = dOmega r cos(phi)` is `dOmega R`, a function of the cylindrical radius alone, so `(du/dZ)_R`
vanishes for any `r` and `phi`. v0.11 said `S = 0` "to round-off"; the field is read through the
file's linear interpolant, whose derivatives are secants, so what is left is that interpolant's
truncation. Measured: `(du/dZ)_R` is 1.637e-08 against `dOmega` of 1e-6, a ratio of 1.637e-02, and
`|S/g|` is **3.180e-05** per radian. **Ruled (section 14.1):** correct, "to round-off" was written
for analytic derivatives; the bound is `|S/g|` below 1e-4 per radian, which check 5 now carries and
meets.

**2. The slope discontinuity was my clamped inverse map, not the specification. Withdrawn.** The
first filing attributed an eightfold degradation at the ends of the anchor's levels to Step 2
decision 8, which continues `Phi(p)` at the slope of its last interval outside them. That is
wrong. Decision 8's continuation is continuous in slope and leaves no kink. What had a kink was my
inverse, `p_of_Phi`, which used `np.interp` and therefore clamps: the mesh reaches one geopotential
spacing past the anchor's levels at each end, so the clamped inverse put several nodes on one
pressure. Measured on the 783 geopotential nodes: the clamped inverse gives **2 duplicate
pressures** and a grid that is not strictly monotonic; the continued inverse gives none and spans
19.68 to 131,400 Pa against the anchor's 19.95 to 129,848.

**Ruled (section 14.2):** the script follows decision 8; no change to the specification's text. The
reviewing agent measures the largest `|S/g|` on the span at 0.050 per radian with the map continued
and 0.136 with it clamped. Corrected here, check 4 measures the whole mesh with no node set held
aside, and the first filing's 4.203e-01 at the ends is gone.

**3. `np.gradient` cannot return an exact zero, and the specification asks for one. Fixed.** Given
a coordinate array `np.gradient` takes its non-uniform path, whose weighted coefficients
`a + b + c` do not cancel in floating point; on a constant field of about -63 with coefficients of
about 57 it returned **9.095e-13** rather than zero. The three-point formula is now written in
difference form in `lib.kernel._three_point_derivative`, algebraically identical and carrying
differences of neighbouring values, so equal columns give an exact zero. Verified: exactly `0.0` on
a constant and 3.608e-15 on a quadratic with uneven spacing. **Ruled (section 14.3):** accepted,
and it is the rule for the composition term. `lib.mesh.d_dlatitude` keeps `np.gradient`, where no
exact zero is required.

**4. A synthetic wind state has to be posed on the geometry the kernel differentiates on.** The
identities `u = beta Z` and `u = dOmega R` are kinematic: they hold for `u` as a function of
position. Posing them through a file whose radius came from one integration while the mesh's radius
came from another makes `u` at the node a different function and the closed form does not apply;
measured before the correction, `(du/dZ)_R` was 42 percent off `beta`. **This was never a property
of `lib.kernel`**, which reproduces every term of the decomposition at a well-posed node to 0.68
percent. The NaN reported with it was the clamp of finding 2, not a separate fault. **Ruled
(section 14.4):** the posing is right; the NaN is the clamp's; strictly increasing nodes are the
schema's business and no guard goes in `lib.windfield` (decision N). Accepted.

**5. The cylinder-extended wind: right about the mechanism, wrong about the measure and the size.**
The first filing reported that decision L and Eq. A15 are in tension at order one and that the
cylinder identity is lost at 7.6e-2. **That conclusion is withdrawn.**

What was right: `u = U(s)` gives `(du/dZ)_R = U' sin cos - U' cos sin = 0`, two terms of about
7e-4 per second cancelling; `U` is built by linear interpolation through the file's `(s,
u_reference)` pairs and read back by the decision L interpolant, so it carries slope
discontinuities; the cancellation fails inside the cells holding a jump; and the resulting `|S/g|`
is spiky, alternates in sign, and its maximum does not fall under refinement. The first filing
measured that maximum correctly and the reviewing agent's independent kernel agrees on it, 5.408e-02
per radian here against 0.050.

What was wrong: the measure. The transfer applies the line integral of the kernel along an isobar,
not its maximum, and `max |S/g|` times the span is not `d ln N`. For an integrand that alternates
in sign the two differ by two orders of magnitude. Measured as the line integral along each of the
anchor's isobars, the largest `d ln N` over the span is **4.746e-04** against a bound of 2e-3 and
the isobar shift **130.3 m2/s2** against 300.

**Ruled (section 14.5, decision Q):** a wind hypothesis on a grid carries a truncation floor and
the model reports it; the cylinder checks of Steps 3, 4 and 5 are stated at that floor as line
integrals with the maximum reported. Neither decision L nor decision P changes.

**What the two kernels agree and disagree on, measured.** On the aggregate quantities they agree:
the largest `|S/g|` 5.408e-02 against 0.050, the isobar shift 130.3 against "up to 224", and both
well inside the 2e-3 and 300 bounds. Per level they do not: `d ln N` is +1.30e-05 at the top against
the stated +1.2e-4, -9.33e-05 at 10 mbar against -5.2e-4, **-1.78e-04 at the gauge against -1.7e-4**,
and +2.21e-04 at the bottom against -1.1e-3, which differs in sign. The gauge level, where both
kernels sit on the reference surface and the residual is largest relative to the cell structure,
agrees to two percent; the others differ by factors of five to nine and one in sign. That is what a
sign-alternating residual does: its integral over a span depends on where the cells fall relative
to the slope jumps, and two implementations that agree on the field will not agree on the
cancellation cell by cell. Both are inside the floor decision Q sets, which is the quantity the
specification bounds. It is recorded rather than reconciled, and no value was adjusted toward the
other.

**A correction to this report's own record.** The first filing reported a `|d ln N|` column computed
as `max |I (phi_a - phi_target)|`, where `I` is the vertical integral of `S/g` in `Phi`; that is not
`d ln N` and gave 1.4e4. I corrected the formula to `max |S/g|` times the span, which was still the
wrong quantity, and reported it three times. The line integral is the right one, and it is what
check 9 now measures. Correcting arithmetic without re-examining whether the quantity was the right
one is the error worth recording.

## 5. Regression

| Suite | Result | Reference |
|---|---|---|
| `step1/accept_step1` | 6 of 6 | 6 |
| `step1/verify_review_changes` | 14 of 14 | 14 |
| `step2/accept_step2` | 8 of 8 | 8 |
| `step3/accept_step3` | 5 of 5 | 5 |
| `step4/accept_step4` | 7 of 7 | 7 |
| `step5/accept_step5` | 7 of 7 | 7 |
| `step6/accept_step6` | 6 of 6 | 6 |
| `step7/accept_step7` | 8 of 8 | 8 |
| `step8/accept_step8` | 9 of 9 | 9 |
| `step9/accept_step9` | 6 of 6 | 6 |
| `step02_1/accept_step02_1` | 6 of 6 | 6 |
| `step02_2/accept_step02_2` | 7 of 7 | 7 |
| `step02_3/accept_step02_3` | 9 of 9 | 9 |
| `step02_4/accept_step02_4` | 9 of 9 | 9 |
| `step02_5/accept_step02_5` | 7 of 7 | 7 |
| `step02_6/accept_step02_6` | 7 of 7 | 7 |
| `step03_1/accept_step03_1` | 9 of 9 | 9 |
| `step03_2/accept_step03_2` | 8 of 8 | 8 |
| `step03_3/accept_step03_3` | 16 of 16 | 16 |
| `step03_4/accept_step03_4` | 13 of 13 | 13 |
| `step04_0/accept_step04_0` | 15 of 15 | 15 |
| `step04_1/accept_step04_1` | 13 of 13 | 13 |
| `step04_2/accept_step04_2` | 13 of 13 | 13 |
| `step04_3/accept_step04_3` | 9 of 9 | this step |

All 23 accepted suites return their reference counts. The registered kind N restored to
SHA-256 `920304440ee4554648c3aa6321b713741ee26eb0eb1ff223f5a81a8f9b20fd37`, the swept value, and the
closure run's inputs and product and the transfer run's inputs all restored equal. Porcelain at the
end holds this step's `lib/kernel.py`, this report, and the author's SPEC_04 v0.12 and
REVIEW_04_step3.

This step adds a new module and changes no existing one, so no accepted behaviour can move; the run
is the procedure's confirmation of that. It was run twice: the first run executed the acceptance
script while it was still being edited under the rulings and recorded a broken row for
`step04_3`, so it was discarded and rerun against the stable script rather than reported.

## 6. Hashes

Nothing was registered. No product of the Step 1 or Step 2 sweeps changes, so
`reports/REPORT_04_step1.md` section 6 and `reports/REPORT_04_step0.md` section 6 remain the record.

## 7. Next step

Step 4: `forward.transfer`, `forward.estimate`, the outer loop, the M = 2 identity test.

Decision Q sets the tolerances the cylinder state is tested at there: the transfer under that wind
is an identity to the truncation floor of the grid the wind is given on, not to round-off, and the
floor on the Lindal state at 0.5 degrees is 4.7e-4 in `d ln N` and 130 m2/s2 in the isobars,
against 1.1e-2 and 31,000 under the closure wind. The closure wind is exempt: its vertical
derivative is exactly zero and the kernel is the closed form to 4.163e-17.

Nothing in this step is left open. The five findings are ruled on in section 14 and applied.

For the record, three defects in this step's own scaffolding were caught by `lib.windfield`
refusing to extrapolate rather than clamping: two coinciding pressure nodes, a latitude subset that
stopped one ulp short of the pole, and a pressure one ulp below the file's grid. A field that
clamped silently would have returned plausible numbers in all three.
