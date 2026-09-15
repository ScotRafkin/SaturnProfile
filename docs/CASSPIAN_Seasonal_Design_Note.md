# CASSPIAN design note: occultations, thermal retrievals, and the seasonal propagator

Version 0.3, 15 September 2026. Author of record: S. Rafkin, Southwest Research Institute.
Status: draft for the author's markup. A design document, separate from the manuscript and from
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
corrected by data. The occultation is the absolute reference at its own season; the retrievals
supply differences and gradients; the propagator supplies tendencies; and the estimate at any
latitude and season is the observed state moved there, with its uncertainty grown by every step
of the move. A season in the middle of a gap between observations still gets an answer; its
envelope says how far the nearest observations are.

The consistency requirement that binds temperature and wind is Eq. A19: on a geopotential
surface the meridional gradient of temperature (and composition) and the axial shear of the wind
determine one another. Whatever the source of the thermal field, the shear is derived from it
through A19, or the thermal field is derived from the shear; never both independently. The
sections below always take the thermal field as primary and derive the shear, because cloud
tracking gives the wind at one level and nothing about its change with season.

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
which every altitude in the model descends; a refractivity free of composition; the structure
below the retrieval floor (to 1.3 bar for Lindal); and a measured column at a known season
against which the whole construction is checked. Its role in the estimate is the anchor of
highest weight, at its own latitude and season. Its weight decreases with distance in latitude,
through the uncertainty of the transfer kernel, and with distance in season, through the
uncertainty of the propagator (§6), but the part of it that is the mean state persists at any
season (§6.3).

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

**The differential principle.** The retrievals are used for what they measure best, which is
change: the gradient with latitude and the change with season. Both are differences in which the
retrieval's absolute biases (its prior and its composition assumption) cancel to first order. The
absolute thermal state remains the occultation's wherever an occultation reaches.

### 4.1 As the source of the shear

Given `T_R(φ, p)` and a composition, the field is placed on geopotential at each latitude by the
hydrostatic relation of Eq. A8,

    Φ(p) = Φ_g − ∫_{p_g}^{p} R̄ T d ln p',                                                  (S4)

with `Φ = 0` on the gauge isobar `p_g`, and the shear kernel follows from Eq. A20 (uniform
composition; the composition terms of A19 are added when the composition varies):

    S / g = −(1/T) (∂T/∂φ)_Φ + (I / T) (∂T/∂Φ),      I(φ, Φ) = ∫_0^Φ (S/g) dΦ',            (S5)

a fixed point in `I` that converges quickly because the second term is small. `S` is the axial
shear `2 Ω_abs r (∂u/∂Z)_R`, from which the sheared component of the wind is integrated given the
barotropic component at the reference level from cloud tracking (§B4, Eq. A39). This is the use
the manuscript already plans; it uses the latitudinal gradient of the retrieval and nothing of its
absolute level.

### 4.2 As an anchor of lower weight: the retrieval reduced to kind N as if it were an occultation

A retrieval column at latitude `φ_i` is reduced to refractivity exactly as Lindal's column was,
under the retrieval's own composition assumption, which plays the part Lindal's 94 percent
hydrogen played (author, 15 September 2026):

    N_i = p ℛ̄_R / (k_B T_R(p)),                                                            (S6)

with `ℛ̄_R` the mean molecular refractivity of the composition the retrieval assumed. The one
thing the retrieval lacks that the occultation has is the altitude column. It is constructed: the
radius of the anchor isobar at `φ_i` comes from the wind-included reference surface of Eq. B3
marched from the occultation's `r0` (the model's own geodesy at that latitude), and the height of
each level above it from the hydrostatic integration of the retrieval's own temperature under its
own composition,

    h(p) − h_ref = −∫_{p_ref}^{p} (R̄_R T_R / |g_eff|) d ln p',                                (S6a)

a fixed point in `|g_eff|` through the radius that converges at once because gravity varies by a
part in a thousand over the column. The retrieval then enters the occultation pipeline as a kind T
`retrieval` instance with constructed heights, and `refrac` reduces it to a kind N product with
everything embedded, as it does for an occultation. The radius enters the transfer only through
`|g_eff|` along the column, so the constructed geodesy is adequate for the physics; what it is not
is a measurement, and the product says so: the heights carry `provenance = "derived"` and the
product carries `altitude_registration = "constructed"`, against `"measured"` for an occultation,
and the radius uncertainty of the reference surface (the geoid's wind sensitivity) is in its
companions.

Reduced this way, retrievals are `M` more anchors of the occultation kind and the transfer in
latitude handles them without a second code path. Their weight is lower: the retrieval's
temperature error, its resolution and the composition dependence of (S6) together are of order a
percent in `N`, against a part in a thousand for an occultation, and the altitude registration is
the model's rather than the retrieval's. The inverse-variance weighting then leans
on the occultation wherever the target is near one in latitude and season and on the retrievals
wherever it is far. A retrieval used both as shear source and as anchor in the same run is not a
contradiction; it is A19 enforced. The independent checks are the same-season occultation
comparison (§7) and the cloud-tracked wind's change between epochs, not the consistency diagnostic
between a retrieval and a shear field derived from it.

### 4.3 As the calibration of the seasonal propagator

The retrievals over the record calibrate the propagator of §5: they supply the measured seasonal
differences where both seasons were observed, and the amplitude and lag that fix the response
model elsewhere. This is the differential principle applied in time.

### 4.4 As the upper boundary

The limb retrievals extend above the top of any occultation. The column above the occultation's
top level is then observed rather than assumed: the limb column is an anchor in the overlap and
the sole source above it, and the boundary pressure `p_b` of (S3b) moves to the retrieval's top,
where by Eq. B7 its influence on everything below is smaller in proportion to the pressure ratio.
Above the limb ceiling nothing is measured and the state is a parameterization with the seasonal
envelope as its uncertainty.

### 4.5 As validation

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

**The shear.** Eq. A20 has an exact form in `N`. With `ln T = ln p − ln N + const` on a surface of
uniform composition and Eq. A17 for the pressure gradient along the surface,

    S / g = (∂ ln N / ∂φ)_Φ + I [ 1 / (R̄ T) + ∂ ln T / ∂Φ ],      I = ∫_0^Φ (S/g) dΦ',          (S9c)

so the shear kernel at the new season is formed from the season-propagated `N` field on `Φ` at
neighboring latitudes, with the temperature in the bracket taken from the production (it
multiplies `I`, which is small, so its composition dependence is second order). This is what
"the shear is propagated by the same field" means in practice: (S9a) and (S9b) applied at every
latitude, then (S9c); temperature and wind are consistent at the new season by construction and
no separate propagation of the wind exists. The physics of §5.2 is written in temperature because
that is where the radiative relaxation lives; its fit and its application are in `ln N`. The radius of the gauge isobar shifts by the
hydrostatic change of the column below it, as (S9) and the paragraph on the anchor radius
describe.

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
2. **Transfer in latitude** along isobars by Eq. A28 with the kernel of Eq. A27, whose shear is
   derived from the season-propagated thermal field by (S5), adding the transfer uncertainty
   `σ_i^transfer(Φ; φ, φ_i)` from the kernel error.
3. **Combine** by the inverse-variance mean of Eq. A30 at the gauge latitude,

       C(Φ) = Σ_i w_i C_i(Φ) / Σ_i w_i,     w_i = 1 / [ (σ_i^meas)² + (σ_i^season)² + (σ_i^transfer)² ],   (S14)

   with the consistency diagnostic `D_ij` of Eq. A33 between every pair.
4. **Produce** pressure, temperature and altitude at the target by (S3), with the boundary
   pressure from the highest anchor (§4.4) and the altitude from the wind-included surface
   through the occultation's `r0`.

The weights do what the principle asks. An occultation near the run in latitude and season is
sharp and dominates; a retrieval is broader but present everywhere; an anchor far in season is
broadened by the propagator until, at a season in the middle of the largest gap between
observations, the propagator's uncertainty is the envelope.

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
   with their assumptions and resolution; a tool that constructs the heights by (S6a) on the
   model's reference surface so that `refrac` reduces the retrieval to kind N by (S6) as it does
   an occultation; the `altitude_registration` attribute on kind N; the thermal-wind path of the
   wind tool by (S5). Specified after SPEC_04, with
   the forward-inputs specification.
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
3. Temperature and wind are consistent at every season by construction, through A19, with the
   thermal field primary and the shear derived.
4. The occultation is the absolute state at its own season; the retrievals are used
   differentially; the seasonal propagator supplies tendencies and is never the starting state.
   Seasons are cyclic: the propagator is one function on the seasonal circle fitted to every
   observation, every target season is bounded by observations on both sides, and the weighting
   toward the observations is strongest where the target is closely bounded by them.
5. No hemispheric symmetry is assumed. The response is fitted per hemisphere; the equinox
   comparison is examined before any symmetry is declared even for the smooth part.
6. The retrievals' limb data serve the upper boundary.
7. A retrieval is reduced to kind N under its own composition assumption, with constructed
   heights on the model's reference surface, and enters the occultation pipeline as one more
   anchor; `N` at fixed pressure is the transportable quantity in time as well as in latitude,
   and the propagator is constructed, fitted and applied in `ln N` (S9a to S9c).
8. The occultation leg of the code proceeds unchanged apart from the season identifier; the
   retrieval leg and the propagator are separate code bases, specified after SPEC_04.
