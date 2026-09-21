# CASSPIAN design note: occultations, thermal retrievals, and the seasonal propagator

Version 0.5, 16 September 2026. Author of record: S. Rafkin, Southwest Research Institute.
Status: draft for the author's markup; v0.4 applied the author's decisions of 16 September 2026
(§10, decisions 9 to 12): no preference by origin among anchors, the wind as a hypothesis whose
shear the anchors correct in the estimate, the temperature field never an input. v0.5 (decision
13): a retrieval is registered in geopotential from its own gauge isobar and carries no absolute
radius; the "constructed heights" of v0.3 are withdrawn. A design document, separate from the manuscript and from
the code specifications; it records how the three sources of information are defined, what each
supplies, how they are combined, and what the decisions of 14 and 15 September 2026 require of
the code now and later. Equation numbers here are S1 onward; A and B numbers refer to the
manuscript's appendices.

---

## 1. The principle

CASSPIAN produces a temperature and pressure profile against altitude at a chosen latitude and a
chosen season, with an uncertainty that is honest about where the information came from. Three
sources feed it, and they are not interchangeable:

- **Occultation profiles** (Lindal et al. 1985; the Cassini radio occultations) are measurements
  of refractivity against absolute planetocentric radius at one latitude and one season. They are
  the only source with absolute altitude registration, and expressed as refractivity they carry
  no composition assumption. They are sparse in latitude and in time.
- **Thermal retrievals** (Cassini CIRS, nadir and limb) are inversions of infrared radiance into
  temperature against pressure over most latitudes and a third of a Saturn year. They carry a
  prior, a vertical smoothing of about a scale height, and the composition assumed in the
  retrieval, and they carry no altitude. They are dense in latitude and time within their range.
- **A seasonal propagator** is a model of how a column's temperature changes between two seasons.
  It is needed wherever the run's season is not the season of the data.

The principle that orders them: **the observations are the state; the model propagates them.**
In latitude the propagator is the transfer kernel of Appendix A, which needs the wind field and
nothing about the thermal state. In season the propagator is the subject of §5, and it is built
the same way: it moves an observed column in time and needs only the seasonal tendency, never a
complete state of its own. The seasonal model is therefore never the starting state to be
corrected by data. A season in the middle of a gap between observations still gets an answer; its
envelope says how far the nearest observations are.

**No preference by origin (author, 16 September 2026).** Every source is reduced to
refractivity under its own assumptions, and from that point the pipeline sees anchors with
declared uncertainties and nothing else. An occultation and a retrieval differ in what they
measured and how well, and those differences are numbers in their uncertainties (the retrieval's temperature carries
its prior and its resolution; the occultation alone carries a measured radius, which feeds the
reference surface and not the anchor's weight); they are not ranks. The weight of an anchor at a target is
its uncertainty as it arrives there: its own, grown by the transfer in latitude along the path and
by the propagation in season across the seasonal circle (§6). A sharp profile far away in
latitude and season pushes weakly; a broader one nearby pushes hard. The model never
differentiates an anchor and never computes a gradient from the profiles.

**The wind is a hypothesis, and the anchors correct its shear.** The wind field the transfer
needs is supplied as a prior: a reference-level wind against latitude from any source (cloud
tracking where it exists, a published table, an invented profile), and a shear along the local
vertical declared in any form the user chooses (none; constant on cylinders; a parameterization
with parameters), each with a declared uncertainty. Within one run the delivered temperature field
is in balance with that wind by construction, because the transfer builds `N` from the kernel and
the kernel is the shear (Eq. A39 holds on every isobar of the output). Where several anchors are
present, their disagreement after transfer to the gauge (Eq. A33) is the data saying what the
integrated shear between their latitudes is, isobar by isobar, and the estimate of Eq. A35 corrects
the kernel within its declared uncertainty until the anchors agree; the corrected kernel is a
corrected shear, the posterior wind, which the run writes back in the user's form. The
reference-level wind is invisible to the kernel (A9, A11) and is never corrected by the anchors.
The temperature field is never an input: the wind is the prior, the anchors are the data, the
estimate reconciles them, and B8's deferred decision is settled this way. Where anchors are dense
the data set the shear and the hypothesis only fills the shape between them; where they are
sparse or absent, below the depth they reach or in an unobserved hemisphere, the hypothesis and
its declared uncertainty stand.

---

## 2. The season identifier

Every file, input or output, carries the epoch of what it describes and the season it belongs to.

The season is the planetocentric solar longitude `Ls`, in degrees, zero at northern spring
equinox. It is defined through the sub-solar planetocentric latitude `δ_s`, which the ephemeris
gives directly:

    sin δ_s = sin ε sin Ls,        ε = 26.73°,                                          (S1)

with the quadrant fixed by the sign of `dδ_s / dt` (Ls in the first or fourth quadrant when the
sub-solar point moves north). `Ls` is a property of the date, so it is computed by a tool from
the date and the ephemeris and transcribed with its source, never typed from memory. For the
Voyager 2 ingress of 26 August 1981 it is about 18° (northern spring, one and a half Earth years
after the equinox of early 1980); the value that goes into the file comes from the ephemeris.

Every file is valid at an instant. Nothing the model reads or writes spans an interval over
which anything varies for the model's purposes: an occultation is a date, a model output is its
run's season, and a retrieval compiled over months is assigned the date its authors assign it,
with the span recorded in a note and carried into its uncertainty only if it ever matters (author,
15 September 2026). Attributes, on every kind:

- `epoch`: an ISO 8601 date, the instant the file is valid at; `epoch_note` (free text) when the
  date stands for a span, saying how it was chosen.
- `solar_longitude_deg` and `solar_longitude_source`, for every file that belongs to a season.
- `season_absent_meaning = "uniform"` on the kinds that have no season by their nature (a
  harmonic set, a rotation system, molecular properties, a composition declared
  season-independent), in the SPEC_00 vocabulary for an absent dimension, so that the absence
  of `solar_longitude_deg` is never silent.
- On a forward product: the run's declared `solar_longitude_deg`, and the season of every anchor
  it used.

A run declares its season in the namelist. In closure mode the run's season equals the anchor's,
because a closure at another season is not a closure. No other logic acts on the season until the
propagator exists; the identifier is put on every file now so that nothing later has to guess.

---

## 3. The occultation: definition and role

An occultation profile arrives as tabulated pressure, temperature and altitude at one latitude and
one date (Appendix B2). The reduction (SPEC_01, SPEC_02) turns it into refractivity on absolute
radius using only the source's own assumptions:

    n = p / (k_B T),      r = r0 + (h − h_ref) cos ψ,      N = n ℛ̄_source,                (S2)

with `r0` the radius of the anchor isobar at the planetocentric anchor latitude, constructed on
the wind-included geoid from the source's own gravity, rotation and wind (Eq. B3 marched from the
declared anchor rule), and `ψ` the tilt of the local vertical there. The altitude `h` is measured
along the local vertical, so its radial projection carries `cos ψ` (SPEC_03 decision 1).

The forward production (SPEC_03) recovers pressure and temperature from `N` under any chosen
composition:

    Φ(h) = ∫ |g_eff| dh  along the local vertical, Φ = 0 on the gauge isobar,             (S3a)
    p(Φ) = p_b + ∫_Φ^{Φ_b} (m̄ / ℛ̄) N dΦ',      T = p ℛ̄ / (k_B N).                        (S3b)

What the occultation supplies that nothing else does: the absolute radius of an isobar, from
which every altitude in the model descends; the structure below the retrieval floor (to 1.3 bar
for Lindal); and a measured column at a known season against which the whole construction is
checked (§7). In the estimate it is one anchor among the `M`, weighted by its declared
uncertainty as it arrives at the target (v0.4, §1): sharp where the target is near its latitude
and season, broadened by the transfer and the propagator elsewhere; the part of it that is the
mean state persists at any season (§6.3).

The Cassini radio occultations are occultations in exactly this sense, at other latitudes and
seasons, and enter the same way once transcribed and reduced.

---

## 4. Thermal retrievals: definition and the four uses

A retrieval is a temperature field `T_R(φ, p)` on a latitude and pressure grid for a date or a
range of dates, with a stated uncertainty, a vertical resolution, and the composition and
ortho-para assumptions of the retrieval. It enters the repository as a kind T file of instance
`retrieval`, transcribed with those assumptions the way Lindal's were, its resolution declared as
part of its uncertainty, and its season by §2. The CIRS record spans 2004 to 2017, from northern
winter through the equinox of August 2009 to about northern summer solstice, in both hemispheres
at every date.

**What a retrieval is to the pipeline (v0.4).** One more anchor, reduced under its own
assumptions, weighted by its declared uncertainty as it arrives at the target, and nothing else.
The two facts about retrievals that the files record are facts, not preferences: a retrieval
carries no absolute radius and is registered in geopotential from its own gauge isobar (§4.2), and
the limb retrievals reach above any occultation's top level (§4.4).

### 4.1 The retrieval as an anchor: the reduction

A retrieval column at latitude `φ_i` is reduced to refractivity exactly as Lindal's column was,
under the retrieval's own composition assumption, which plays the part Lindal's 94 percent
hydrogen played (author, 15 September 2026):

    N_i = p ℛ̄_R / (k_B T_R(p)),                                                            (S6)

with `ℛ̄_R` the mean molecular refractivity of the composition the retrieval assumed.

### 4.2 Registration in geopotential; no absolute radius (v0.5)

A retrieval has no altitude information, and the model lends it none. Its levels are placed in
the model's `(φ, Φ)` plane from its own temperature and composition by Eq. A8,

    Φ(p) − Φ(p_g) = −∫_{p_g}^{p} R̄_R T_R d ln p',      Φ(p_g) = 0 on its own gauge isobar,     (S6a)

wherever that isobar sits in radius. This is the gauge every anchor uses, so the retrieval enters
the transfer as the set of points `(φ_i, Φ, N)` asserting only the thermal structure above and
below its gauge isobar, which is what it measured, and nothing about where that surface was. If
its epoch differs from the run's, the propagator of §5 moves `N` at fixed pressure and carries the
hydrostatic shift of its gauge isobar (the paragraph on the anchor radius), which is why the anchor
carries `Φ` and not `r`. No radius enters the reduction of a retrieval: (S6a) needs the
retrieval's temperature and composition and nothing else, since gravity is absorbed in the
coordinate. Gravity and radius enter only on the run's side, in the kernel through the run's own
column and reference surface, the geometry every anchor is transferred through. The kind T
`retrieval` instance therefore carries `Φ(p)` in place of a height column, and the kind N product
for a retrieval carries geopotential as its level coordinate in place of radius and
`altitude_registration = "none"` against `"measured"` for an occultation (a schema decision the
retrieval-leg specification makes explicitly).

Absolute altitude is a property of the delivered state, not of any anchor. It descends from the
run's reference surface, whose constant of integration comes from the occultations' measured
radii: one now, Lindal's; with several, their radii at their latitudes and epochs over-determine
the constant and it is estimated with weights, the same structure as the anchor constant on an
isobar. Retrievals contribute nothing to it, and the delivered altitude far from any occultation
carries the geoid's uncertainty accordingly.

Reduced this way, retrievals are more anchors of the same kind and the transfer in latitude
handles them without a second code path; their uncertainty (temperature error, resolution, the
composition dependence of (S6)) is declared in the product and is what the estimate weighs.

### 4.3 The record as the calibration of the seasonal propagator

The retrievals over the record calibrate the propagator of §5: they supply the measured seasonal
differences where both seasons were observed, and the amplitude and lag that fix the response
model elsewhere. An occultation at a season the record covers enters that calibration on the same
footing, as one more observation at its season.

### 4.4 The limb retrievals as the upper boundary

The limb retrievals extend above the top of any occultation. The column above the occultation's
top level is then observed rather than assumed: the limb column is an anchor in the overlap and
the sole source above it, and the boundary pressure `p_b` of (S3b) moves to the retrieval's top,
where by Eq. B7 its influence on everything below is smaller in proportion to the pressure ratio.
Above the limb ceiling nothing is measured and the state is a parameterization with the seasonal
envelope as its uncertainty.

### 4.5 Validation

At the delivered latitude and season, a retrieval not used as an anchor in that run is compared
with the delivered temperature (SPEC_05). The same-season comparison of §7 is the first such test.

---

## 5. The seasonal propagator

### 5.1 Definition

Write the thermal field at a latitude and pressure as a mean and a seasonal anomaly,

    T(φ, p; s) = T̄(φ, p) + A(φ, p; s),                                                     (S7)

with `s` the season (`Ls`, or the date). The propagator from season `s_a` to season `s` is the
difference of the anomaly,

    𝒫[s_a → s]: T(s) = T(s_a) + A(s) − A(s_a) ≡ T(s_a) + ΔA.                                (S8)

The mean `T̄` never enters the propagation. It is whatever the data say it is: the occultation's
column, or a retrieval's. This is the sense in which the seasonal model is a propagator and not a
state. A model that supplied `T̄` as well would be a full radiative model and would replace the
data with its own biases; that is not built.

Applied to an anchor, the propagator acts on `N` at fixed pressure. Since `N = p ℛ̄ / (k_B T)`
with `ℛ̄` fixed for the anchor, the seasonal change of `N` at a given pressure is the log ratio of
the temperatures, and when the change is taken from the retrieval record the composition cancels
exactly:

    Δ ln N(φ, p; s_a → s) = −ln [ T(s) / T(s_a) ] = −ln [ 1 + ΔA / T(s_a) ],                (S9)

so `N` at fixed pressure is the transportable quantity in time as well as in latitude (author,
15 September 2026). The propagator is a field `𝒫(φ, p; s_a, s) = Δ ln N` on latitude, pressure and
season, and everything about its construction and its application can be written without leaving
`N`:

**Construction.** Seasons are cyclic, and the propagator is one function on the seasonal
circle, not two cases (author, 15 September 2026). Every observed season contributes a point:
the retrieval-derived refractivity (§4.2) `ln N_R(φ, p; s)` at each date of the record, and the
occultations at theirs. The anomaly `a(φ, p; s) = ln N_R − ⟨ln N_R⟩` is fitted over the whole
circle by one of the forms of §5.2, and `𝒫 = a(s) − a(s_a)` from the fit. At an observed season
the fit reproduces the observation within its residual, so the propagator between two observed
seasons is the measured difference of `ln N_R`, in which the retrieval's composition assumption
cancels (the same retrieval family at both seasons; otherwise the ratio `ℛ̄_R(s) / ℛ̄_R(s_a)` is
divided out). Between observed seasons the fit supplies the shape, which is what a cyclic
response has and a straight interpolation across a gap does not; the bounding observations fix
its amplitude and phase. The fit is made in `ln N` rather than in temperature: radiative
relaxation is linear in the temperature anomaly, and `ln T` and `ln N` differ by a constant at
fixed pressure, so `a` obeys (S10) with its own gain and the same time constant. No temperature
conversion enters the propagator at any point.

**Application to an anchor.** At each tabulated pressure of the anchor's kind N product, the
propagator is interpolated in `ln p` and latitude to the anchor's levels and applied,

    N_s(p) = N_a(p) exp[ 𝒫(φ_a, p; s_a, s) ],                                                (S9a)

and the geopotential is re-formed from the hydrostatic relation in its `N` form (from `dp = −ρ dΦ`
with `ρ = (m̄ / ℛ̄) N`),

    Φ_s(p) = Φ_g − ∫_{p_g}^{p} (ℛ̄ / m̄) (p' / N_s(p')) d ln p',                                (S9b)

integrated from the gauge isobar with the run's composition, exactly as the production does. So
the run's composition enters the season-propagated anchor only where it enters every anchor,
through the geopotential mapping and the production, and the column expands or contracts on `Φ`
as it warms or cools. The propagated column is then `N_s(Φ)`, ready for the transfer in latitude.

**The anchor radius.** `r0` is the radius of the anchor isobar at `s_a`. At another season that
isobar sits at `r0 + δr0`, where `δr0` is the hydrostatic expansion of the column between the
seasonally invariant depth (below about a bar, where the record shows no change within its
noise) and the anchor isobar, computed from the same `𝒫` field with the composition by (S9b) in
altitude. For an anomaly of a few kelvin over a scale height it is of order a kilometer, small
against the anchor radius uncertainty but not negligible against the altitude precision the
model aims at, and it is carried into the propagated anchor and its uncertainty.

**The shear at the new season (v0.4).** No shear is derived from the propagated thermal field.
The wind hypothesis (§1) is the prior at every season (its reference-level wind may itself be
declared seasonal), and the season-propagated anchors correct its shear in the estimate of §6
exactly as same-season anchors do. Temperature and wind at the new season are consistent by
construction because the delivered field is built from the corrected kernel. The physics of §5.2
is written in temperature because that is where the radiative relaxation lives; its fit and its
application are in `ln N`. The radius of the gauge isobar shifts by the hydrostatic change of the
column below it, as (S9) and the paragraph on the anchor radius describe.

### 5.2 Forms of the anomaly, each a declared rule

The anomaly over the seasonal circle is fitted to every observation by one of the following
forms; the choice is a declared rule in the propagator's control file with the others as
alternatives for the Monte Carlo:

1. **Linear response to the insolation.** The anomaly is the response of a radiatively relaxing
   layer to the seasonal forcing:

       τ dA/dt + A = G F'(φ, t),                                                            (S10)

   with `F'` the diurnal-mean insolation anomaly about its annual mean at that latitude (§5.3),
   and two parameters per latitude, pressure and hemisphere: a gain `G(φ, p)` and a time
   constant `τ(φ, p)`. For a forcing written as a Fourier series in the Saturn year,
   `F' = Σ_n F_n cos(nω(t − t_n))`, the periodic solution is

       A = Σ_n G F_n (1 + (nωτ)²)^{−1/2} cos(nω(t − t_n) − arctan(nωτ)).                     (S11)

   `G` and `τ` are fitted from the observations: the amplitude and the lag of the observed cycle
   at each latitude and pressure fix them. This requires no radiative transfer; the only
   physics computed is the forcing, which is geometry. In the limits the model is transparent:
   `τ` much longer than a Saturn year gives `A → 0` (the deep troposphere, no seasonal change);
   `τ` much shorter than a season gives `A → G F'` (the upper stratosphere tracks the sun). The
   band where `τ` is comparable to a season, the lower and middle stratosphere, is where the lag
   is large and the fit matters, and it is the band the record constrains best.
2. **Harmonic fit with the period fixed.** The empirical alternative to 1: the anomaly fitted as
   a Fourier series in `Ls` with the Saturn year as the period, no forcing and no time constant.
   It is the simplest rule near the observed seasons and it does worst in the middle of the
   largest gap; the disagreement between 1 and 2 is a measured uncertainty.
3. **No change**, the limiting case, for pressures where the observations show none within
   their noise.

No hemispheric symmetry is assumed in any of these. The seasonal forcing is nearly symmetric but
the response is not: the thermal structure at the scale of the belts and zones is tied to the jets
through A19, and the jets differ between hemispheres. `G`, `τ` and the residual are fitted per
hemisphere and per latitude. The record contains both hemispheres at every date, and at the
equinox of 2009 both sat under the same sun; the difference between them at that date is the
measured asymmetry of lag and dynamics, and it is the first thing to be examined before any
symmetry is declared even for the smooth part of the response.

### 5.3 The forcing

The diurnal-mean insolation at the top of the atmosphere at planetocentric latitude `φ` is

    F(φ, Ls) = (S_0 / d(Ls)²) (1/π) [H sin φ sin δ_s + cos φ cos δ_s sin H],
    cos H = −tan φ tan δ_s   (clamped to [−1, 1]),                                          (S12)

with `δ_s` from (S1), `d` the heliocentric distance from the orbit (eccentricity about 0.056,
perihelion near `Ls` 280°), and `H` the half-day. Ring shadowing multiplies the direct beam by a
transmission factor evaluated along the ray to the sun: from a point at radius `r` and latitude
`φ` with the sun at hour angle `h`, the ray crosses the ring plane at a cylindrical radius

    ρ² = (r cos φ + t cos δ_s cos h)² + (t cos δ_s sin h)²,    t = −r sin φ / sin δ_s,        (S13)

for `t > 0` (the winter hemisphere only); if `ρ` lies within a ring of normal optical depth
`τ_ring`, the transmission is `exp(−τ_ring / |sin δ_s|)`. The diurnal mean of the shadowed beam
is a numerical integral over `h`. Ring radii and optical depths are declared data with their
citation. The forcing is what makes the two hemispheres' seasons differ even before the dynamics
do: perihelion falls in southern summer, and the ring shadow falls on the winter hemisphere at
latitudes set by `δ_s`.

### 5.4 Uncertainty of the propagator

Three terms, all declared and all drawn by the Monte Carlo wrapper:

- the parameter uncertainty of the rule (the covariance of `G` and `τ` from the fit, or of the
  harmonic coefficients), which grows with the angular distance on the seasonal circle to the
  nearest observations on either side, smallest at an observed season and largest at the middle
  of the largest gap, because the observations fix amplitude and phase only through the shape
  the rule supplies;
- the interannual floor: two observations of the same season in different years differ, and that
  difference (Lindal against the CIRS record at the same `Ls`, and the CIRS record against
  itself where it repeats a season) is the uncertainty that remains even at an observed season;
- the residual `ε(φ, p; s) = T_R − T̄ − A` of the rule against the record, which is everything
  the radiative picture leaves out: the equatorial oscillation, meridional overturning, storm
  aftermaths. It is carried as a stochastic term with the amplitude and the correlation time the
  record shows, the way the gravity wave statistics are carried;
- the choice of rule itself, a discrete alternative.

Where the residual is large the envelope is wide, and the fit says where that is, latitude by
latitude and level by level. Because the circle is closed, no target season is unbounded: the
question is only how far it lies from the observations on each side, and the weighting toward
the observations is strongest where a target is closely bounded by them.

---

## 6. The estimate at a latitude and season

### 6.1 Assembly

For anchors `i` at `(φ_i, s_i)` with refractivity columns `N_i(Φ)` and uncertainties
`σ_i^meas(Φ)`:

1. **Propagate in season**, (S8) and (S9), to the run's season `s`, adding the propagator's
   uncertainty `σ_i^season(Φ; s, s_i)`.
2. **Transfer in latitude** along isobars by Eq. A28 with the kernel of Eq. A27 formed from the
   wind hypothesis, adding the transfer uncertainty `σ_i^transfer(Φ; φ, φ_i)` accumulated along
   the path from the hypothesis's declared shear uncertainty (Eq. A31).
3. **Combine and correct** at the gauge latitude. With one anchor the constant is fixed and the
   hypothesis stands. With several, the inverse-variance mean of Eq. A30,

       C(Φ) = Σ_i w_i C_i(Φ) / Σ_i w_i,     w_i = 1 / [ (σ_i^meas)² + (σ_i^season)² + (σ_i^transfer)² ],   (S14)

   with the consistency diagnostic `D_ij` of Eq. A33 between every pair, and the estimate of
   Eq. A35: the kernel used is `K + η`, with `η(φ, Φ)` a correction whose prior covariance is
   the declared shear uncertainty with a declared correlation length in latitude, solved per
   isobar as the linear least-squares problem in `C` and `η` given the `M` values `ln N_i`. The
   posterior field's gradient on every isobar is `K + η`, and `η` integrated from the reference
   level is the posterior shear, written back as the wind the run actually used. A misfit that
   exceeds the declared uncertainties is reported and inflates them, as A10 says.
4. **Produce** pressure, temperature and altitude at the target by (S3), with the boundary
   pressure from the highest anchor (§4.4) and the altitude from the wind-included surface
   through the occultation's `r0`.

The weights do what the principle asks. An anchor near the target in latitude and season
arrives sharp and pushes hard, whatever its origin; one far away in either arrives broadened by
the transfer and the propagator and pushes weakly, however precise it was where it was measured;
at a season in the middle of the largest gap between observations the propagator's uncertainty is
the envelope. The two things that cannot be constrained by any number of anchors are the
reference-level wind, which the kernel does not see, and the attribution of a correction between
shear and a composition gradient the run did not supply, since A39 constrains the product `R̄ T`;
the record says under which composition the correction was made.

### 6.2 What persists and what decays

The propagated anchor at a far season is `T̄ + A(s)` with `T̄` from the data and `A` from the
rule. So the occultation's contribution at a far season is its mean, which is the persistent part
of what it measured; the anomaly is the model's. In the language of assimilation, the observation's
correction to the model splits into a persistent bias correction, which carries to any season, and
a transient part, which decays with the dynamical correlation time. The propagated-observation
form used here and the corrected-background form give the same estimate for a linear anomaly model;
they differ only in that the corrected-background form needs the model to supply a mean state and
this one does not.

### 6.3 The limits

- **Deep troposphere.** `τ` is longer than a Saturn year; `A` is small at any season; the
  occultation is the state everywhere in time. This is where the probe's numbers mostly come from
  and where the seasonal machinery has almost nothing to do.
- **Middle stratosphere.** The anomaly is of order ten kelvin; a response model calibrated on a
  third of a cycle reaches a far season with an error of a third to a half of that. Not decisive,
  not useless.
- **Upper stratosphere and above the limb ceiling.** The anomaly is largest and least calibrated;
  above the ceiling nothing is measured and the guess is the state. The envelope there is the
  amplitude itself.
- **Arrival season.** The next northern spring equinox is 2039. A probe arriving in the early
  2040s arrives in the season the CIRS record covers from 2009 onward and the season Lindal
  measured, one and two Saturn years earlier. For that case the target is closely bounded by
  observations on the circle, the weighting is almost entirely theirs, and the honest uncertainty
  is the interannual floor, the difference between the two observations of the same season. The
  width of the delivered envelope depends on how closely the arrival season is bounded by
  observations, and that should be said to the mission as such.

---

## 7. The same-season check

Voyager 2 measured at about `Ls` 18°; Cassini observed `Ls` 0° to about 25° in 2009 to 2011,
one Saturn year later. Lindal's column transferred in latitude on a shear field from cloud winds,
without any seasonal propagation, compared with the CIRS field of 2010 to 2011 at 36° N, tests the
transfer, the production and the retrieval-to-anchor conversion with season out of the way.
Whatever residual remains is interannual variability and retrieval bias together, with the 2010
storm as the known intruder. It is the first end-to-end test of SPEC_05 and the one that
calibrates how much to trust the rest.

---

## 8. What this requires of the code, now and later

### 8.1 Now, in the occultation leg (SPEC_00 v0.17, applied at SPEC_03 Step 3)

The occultation leg and its transfer in latitude (SPEC_03, SPEC_04) do not change. One thing is
added now because it touches every file and is cheap now and expensive later: the season
identifier of §2. Concretely:

- SPEC_00 §5 gains the required global `epoch` (with the optional `epoch_note`), and
  `solar_longitude_deg` with `solar_longitude_source` on every kind that belongs to a season;
  kinds G and R carry `season_absent_meaning = "uniform"` instead; kinds T, D, W and N carry the
  observation's season; kind C carries a season or `uniform` as the composition declares; kind
  `profile` carries the run's season and the seasons of its anchors.
- `lindal_scalars.toml` gains the observation date (26 August 1981, Fig. 2 caption) and the
  season, the latter computed from the ephemeris's sub-solar latitude by (S1) and transcribed with
  its source. The wind's season is that of the Voyager imaging (1980 to 1981); the harmonic set
  and the rotation are `uniform`.
- The namelist `[run]` gains `solar_longitude_deg`; in closure mode it must equal the anchor's.
- The reduction chain is rebuilt once to carry the attributes, in the same in-memory candidate
  procedure as Step 0, at Step 3 where the schema is already being touched for kind `profile`.

No logic acts on the season in the occultation leg beyond the closure equality. Nothing else in
this note changes SPEC_03 or SPEC_04.

### 8.2 Later, as their own code bases, in order

1. **The retrieval leg.** Transcription of the CIRS record (§9) into kind T `retrieval` files
   with their assumptions and resolution; the kind T `retrieval` instance carrying `Φ(p)` by
   (S6a) from its own gauge isobar, so that `refrac` reduces it to kind N by (S6) with no
   radius; kind N with geopotential as its coordinate for that instance; the
   `altitude_registration` attribute on kind N (`"measured"`, `"none"`). Specified after SPEC_04,
   with the forward-inputs specification. **The combination and the estimate** (§6.1 step 3, per
   isobar, with the posterior wind written back; kind W carrying the declared shear uncertainty
   the prior covariance needs) is its own specification, after SPEC_04 and before the
   propagator, since it is what makes the profiles and the wind hypothesis consistent.
2. **The seasonal propagator.** The forcing function (S12, S13) with the ring data; a fitting tool
   that reads the transcribed record and writes the propagator file (the rule, `G`, `τ` or the
   harmonic coefficients, the residual, per hemisphere and latitude and pressure); a new input
   kind for that file with a season dimension; the application of (S8), (S9) and the shear
   propagation inside the transfer; the anchor weights of (S14). Specified after the retrieval leg.
3. **SPEC_05.** The same-season check (§7), the identity tests, the composition scaling, and
   runs under chosen inputs, including runs at seasons outside the record with their envelopes.

Each of these is a set of tools that write files the model reads, plus the weighting in the
transfer. The model itself gains the season-propagation step and the weights and nothing else;
every assumption in §5 lives in a control file where it can be replaced without touching the code.

---

## 9. The retrieval record to transcribe

To be transcribed from the papers, each with its coverage in pressure, latitude, dates and `Ls`,
its resolution, its uncertainties, and its composition and ortho-para assumptions; the list is
what to read, not a statement of what each contains:

- Fletcher et al. (2010), Icarus 208, 337: seasonal change 2004 to 2009 from nadir CIRS.
- Sinclair et al. (2013), Icarus 225, 257: temperatures and hydrocarbons 2005 to 2010.
- Sylvestre et al. (2015), Icarus 258, 224: stratospheric temperatures and hydrocarbons from limb
  observations across the record.
- Fletcher et al. (2016), Icarus 264, 137: tropospheric temperatures, winds and para-hydrogen.
- Fletcher et al. (2018), in Saturn in the 21st Century (Cambridge): the seasonally changing
  atmosphere over 2004 to 2017, the synthesis.
- Guerlet et al. (2009), Icarus 203, 214: limb temperatures and hydrocarbons; Guerlet et al.
  (2014), Icarus 238, 110: a seasonal radiative model, the published alternative to §5.2 rule 2;
  Guerlet et al. (2018) on the equatorial oscillation from limb data.
- Fletcher et al. (2017), Nature Astronomy 1, 765: the 2010 storm's disruption of the equatorial
  oscillation, which bounds the residual of §5.4.
- Friedson and Moses (2012), Icarus 218, 861, and Greathouse et al. (2008): dynamical and
  radiative seasonal models, the reference for what §5 leaves out.

And the Cassini radio occultations (Schinder et al. 2011, 2015 and later), which are anchors of
the occultation kind at other latitudes and seasons.

---

## 10. Decisions recorded (author, 14 and 15 September 2026)

1. Every input and anchor carries its epoch (one date; everything is valid at an instant) and its
   season; a run declares its season; every output carries it.
2. The model never refuses an anchor for seasonal distance. It produces its best estimate for any
   selected season from the data it has, with the seasonal guesses carried as uncertainty in the
   Monte Carlo.
3. Temperature and wind are consistent at every season by construction: the delivered field is
   built from the kernel, and the kernel is the wind hypothesis's shear corrected by the anchors
   (superseded wording of 15 September, "thermal field primary and shear derived", withdrawn on
   16 September).
4. The seasonal propagator supplies tendencies and is never the starting state. Seasons are
   cyclic: the propagator is one function on the seasonal circle fitted to every
   observation, every target season is bounded by observations on both sides, and the weighting
   toward the observations is strongest where the target is closely bounded by them.
5. No hemispheric symmetry is assumed. The response is fitted per hemisphere; the equinox
   comparison is examined before any symmetry is declared even for the smooth part.
6. The retrievals' limb data serve the upper boundary.
7. A retrieval is reduced to kind N under its own composition assumption and enters the
   occultation pipeline as one more anchor (its registration per decision 13); `N` at fixed pressure is the transportable quantity in time as well as in latitude,
   and the propagator is constructed, fitted and applied in `ln N` (S9a, S9b).
8. The occultation leg of the code proceeds unchanged apart from the season identifier; the
   retrieval leg and the propagator are separate code bases, specified after SPEC_04.
9. (16 September 2026) No preference by origin. After reduction to refractivity an anchor is an
   anchor; its weight is its declared uncertainty as it arrives at the target, grown by the
   distance in latitude and in season. The model never computes a gradient from the profiles. The
   "differential principle" of v0.3 is withdrawn.
10. (16 September 2026) The wind is a hypothesis: a reference-level wind from any source and a
    shear along the local vertical declared in any form, each with a declared uncertainty, in the
    user's coordinates; the model builds the field on its own geometry and forms the axial shear;
    for a shear declared constant on cylinders the model writes back what that implies along the
    local vertical. The formulation imposes no constraint on the shear; a vanishing cylindrical
    shear at depth is a hypothesis a user may declare, not a rule.
11. (16 September 2026) The anchors correct the shear in the estimate (A35), per isobar, within
    the declared uncertainty; the reference-level wind is never corrected; the posterior wind is a
    product of every run with more than one anchor; the temperature field is never an input.
12. (16 September 2026) The combination and the estimate are a specification of their own, after
    SPEC_04.
13. (16 September 2026) A retrieval carries no absolute altitude and is lent none: it is placed in
    geopotential from its own gauge isobar by Eq. A8 with no radius entering, and its product
    carries geopotential as its coordinate and `altitude_registration = "none"`. Absolute altitude belongs to
    the delivered state through the reference surface, whose constant is set by the occultations'
    measured radii (estimated with weights when there are several). "Constructed heights" (v0.3
    to v0.4) withdrawn.
