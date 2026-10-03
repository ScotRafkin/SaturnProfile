# CASSPIAN design note: where the variables live, and a staggered arrangement

Reviewing agent, 3 October 2026. **Tabled by the author on 3 October 2026**: the realistic shear
case and the first comparison with CIRS come first, on the current numerics. This note is kept for
when the comparison or the Monte Carlo reopens the question. SPEC_07 was closed, not adopted, at
v0.5 and amended at v0.6.

Read from the working tree on the author's machine: `forward/production.py`, `forward/transfer.py`,
`forward/estimate.py`, `lib/mesh.py`, `lib/kernel.py`, `lib/hydrostatic.py`, `lib/geopotential.py`,
SPEC_00 §7.3 and SPEC_04 §4. Nothing was run. Every claim about behaviour below is a reading of
the code and of the equations it implements, and is for the coding agent to confirm before any
specification is written from it.

---

## 1. What the specifications asked for, and what was built

SPEC_00 §7.3 fixed two numerical rules before any integral was coded:

- (i) a quantity near exponential in the coordinate (N, n, rho, p) is log-linear across a layer;
  T, u and composition are linear;
- (ii) **every variable and every derivative is assigned a location on a staggered (phi, Phi)
  grid, at levels or half-levels and at latitudes or half-latitudes, in a table written and
  agreed before the first integral is coded.**

Rule (ii) was answered by SPEC_04 §4 in one sentence: "Every quantity lives on the `(phi_i, Phi_j)`
nodes", with `S`, `I`, the composition term and the curves listed as exceptions only in how they
are formed, not where they live. That is a collocated arrangement, Arakawa A in latitude and
Lorenz in the vertical. The table that rule (ii) asked for, with a deliberate choice of where each
variable sits, was not written. SPEC_03 Step 1 did stagger one thing correctly (the geopotential
increment on the layer, the gravity on the levels), and that is the only place it was done.

## 2. Where everything lives today

The model works on four grids at once.

| Grid | Coordinates | Set by |
|---|---|---|
| the wind file | latitude nodes (0.5 deg) by pressure nodes (ten per decade) | the hypothesis's author |
| the anchor's levels | 66 tabulated levels, each an isobar with label `p_k` | the occultation |
| the mesh | latitude (0.05 deg) by geopotential (5e4 m2/s2) | the namelist |
| the composition file | its own levels by latitude | the composition tool |

| Quantity | Lives at | Formed by |
|---|---|---|
| `u`, `du/dphi`, `du/dln p` | anywhere, read from the wind file's interpolant | decision L (or L2) |
| `r`, `z`, `g`, `|g_eff|` | mesh nodes | the column march, RK4 in `Phi` |
| `ln p` on the mesh (the isobar map) | mesh latitudes; linear in `Phi` between the curves' knots | the traced curves |
| isobar slopes | mesh nodes | centred differences of the map |
| `S/g` | anywhere; at mesh nodes for the record | Eq. A15, wind at the point, geometry bilinear |
| `I = int S/g dPhi` | any `Phi` in a mesh column | Gauss with the wind's pressure nodes as breakpoints (SPEC_06) |
| `Phi_k(phi)`, the isobar of level `k` | mesh latitudes, one per anchor level | RK4 of `dPhi/dphi = -I` (A24) |
| `ln N_k(phi)` | mesh latitudes, one per anchor level | trapezoid of `K` on the curve, `S/g` read at the label (A27, A28) |
| `rho`, `p`, `T` at the target | the anchor's levels | the production: log-linear density, `p` summed from the top (B4 to B6) |

The isobars are a Lagrangian vertical coordinate: each anchor level is a surface of fixed pressure
that the model follows in latitude. Two quantities are carried on every one of them, `Phi_k` and
`ln N_k`, at the same points.

## 3. The three balances, and why the identity exists

The model carries three relations.

1. **The isobar's height** (A16, A24). Along isobar `k`, `dPhi_k/dphi = -I(phi, Phi_k)`, with `I`
   the vertical integral of `S/g` from the reference surface. This is the gradient wind balance.
2. **The density on the isobar** (A27, A28). Along isobar `k`, `d ln N_k/dphi = K`, with `K` the
   shear term at the isobar's label plus the composition term. This is the thermal wind balance.
3. **The hydrostatic column at the target** (B4 to B6). From `ln N_k` and `Phi_k`, the density is
   taken log-linear across each layer, the layer masses are summed from the top, and that gives
   `p` and `T` on the levels.

In the continuum, 2 and 3 follow from 1. The thickness of layer `k+1/2` is
`Phi_k - Phi_(k+1)`, and by 1 it changes along the isobars as

    d(Phi_k - Phi_(k+1))/dphi = -[I(Phi_k) - I(Phi_(k+1))] = -int over the layer of S/g dPhi

the shear integrated through the layer. Relation 2 with the hydrostatic equation gives the same
change only if the shear inside the layer is what its two level values imply.

Discretely they are not the same.

- **The tracing** integrates `S/g` through each layer at the wind file's resolution, with breakpoints
  and Gauss points (SPEC_06).
- **The transfer** samples `S/g` at the two labels only.
- **The production** then rebuilds the layer from those two samples by a log-linear rule.

The pressure identity, `p/p_label - 1`, is the difference between the two discretizations of one
layer integral, summed from the top. The check 7 law is that difference for a step inside a layer:
the layer mass, times the jump, times where in the layer it falls. Two consequences that SPEC_06
and SPEC_07 ran into come out of this directly:

- **Placement luck.** The size depends on where the change falls within the anchor layer, so run 6
  was small under L by luck and large under L2.
- **The offset carried down.** Pressure is summed from the top, so an error in one layer's mass is
  carried to every level below it. That is the uniform offset below `p_s` in the dense runs.

One more observation. The isobar map is `ln p` linear in `Phi` between isobars, and by hydrostatics
that is an isothermal layer. So the map already holds a layered state: one temperature per layer,
`T = -dPhi / (R_specific d ln p)`, consistent with the traced heights and the fixed labels. The
production at the target computes a second, level-based state from `ln N`. The identity is the
disagreement between the two.

## 4. A staggered arrangement for this model

The Arakawa C idea, carried into this model's coordinates: place each variable where the balance
that defines it is centred, so that no balance needs an interpolation.

**Vertical: Charney and Phillips on the isobars.**

| Where | What |
|---|---|
| interfaces, the isobars `k` | the label `p_k` (exact by definition), the height `Phi_k(phi)` (traced), `u` read there for the column |
| layers `k+1/2` | the thickness `Phi_k - Phi_(k+1)`, the layer temperature `T_(k+1/2) = m_(k+1/2) (Phi_k - Phi_(k+1)) / (R ln(p_(k+1)/p_k))`, the layer composition |

- **Hydrostatic.** The layer's thickness and its mean temperature are the same quantity in
  different units. A change confined to part of a layer moves the layer mean and nothing else:
  there is no "where in the layer".
- **Thermal wind.** The change of a layer's thickness along the isobars is the shear integrated
  through that layer, exactly the integral the tracing already computes. The thermodynamic state is
  transported by the tracing itself. No separate transfer of `ln N` is needed for it.
- **The pressure identity is exact by construction.** The pressures at the target are the labels,
  since nothing integrates them from the top; there is no identity to fail and no offset to carry
  down.

**Latitude.** `Phi_k` sits on the mesh latitudes. The slope between them, `-I`, is evaluated by RK4
at the half steps, where the C grid would put it, and the wind enters `I` read at the point. The
horizontal is already arranged this way in effect, and nothing here changes it. The F8 staircase
was `S/g` drawn at the nodes; its integral across each cell was always smooth.

**What is delivered on the levels.** Temperature and refractivity at the levels are diagnosed from
the layers. This is the one place an interpolation remains, the "none perfect" step, and it is
local, stated and bounded by the neighbouring layer values. Two rules are worth weighing:

- **A.** The level temperature is the mean of its two layers'. Simple and second order, but it
  replaces the anchor's own level structure with layer means even at the anchor's latitude.
- **B.** The level temperature is the anchor's own level value plus the transferred change,
  interpolated to the level from the changes of its two layers. At the anchor's latitude it returns
  the data exactly, and only the change carried by the hypothesis is reconstructed.

I would recommend B. In both, `N` at a level follows from `T`, the label and the composition, and
`p` at a level is its label.

**At the anchors.** Nothing changes. The placement already produces each anchor's labels and
heights from its own levels by the log-linear rule of SPEC_00 §7.3 (i). That is the one point where
data on levels becomes layers, and it is where the data lives.

## 5. What it does to the problems of SPEC_05 to SPEC_07

- **The identity:** gone, by construction rather than by tolerance.
- **A sharp change in the wind.** A layer receives the hypothesis's own shear, integrated through
  it at the wind file's resolution.
  - A step inside a layer is delivered as the layer mean.
  - There is no overshoot, no luck of placement and no offset carried below.
  - The limit is stated and general: structure thinner than a layer is averaged over the layer.
- **The interpolant matters much less.** A layer's thickness change is the integral of the shear.
  For every wind interval that lies wholly inside a layer, that integral is the difference of the
  node values, the same under L and L2. The two readings differ only for intervals cut by a layer
  boundary. L2 still helps there and in latitude, but SPEC_07's question shrinks to a refinement,
  to be judged on the general test of section 6.
- **A general error estimate, for nothing.** Keep the present transfer of `ln N` along each isobar,
  demoted to a diagnostic. Its disagreement with the layered state is exactly the sub-layer
  structure in the hypothesis that the layers cannot carry. Computed on every run, for any
  hypothesis, it is the resolvability measure discussed on 3 October:
  - one physical tolerance;
  - levels beyond it flagged in the product;
  - no bound set per case.
- **Refinement, which absorbs SPEC_08.**
  - Sub-isobars can be inserted between the anchor's levels, with the anchor's state inside its own
    layers given by its log-linear rule. The data carry no information inside a layer; the
    hypothesis may, and the sub-isobars carry it.
  - The diagnostic says where refinement is needed, so it can be added where the hypothesis is
    sharp rather than everywhere.

## 6. Testing on general terms

Every acceptance below is stated for any hypothesis, and the SPEC_05 experiments, the dense
references and randomly drawn winds are samples of it, not owners of their own bounds.

1. **The layered state is consistent.** Thickness, layer temperature and labels agree with the
   hydrostatic relation to round-off on every layer of every run.
2. **The diagnostic predicts the error.** On a set of winds of increasing sharpness, the
   disagreement the diagnostic reports matches the change delivered when the isobars are refined.
3. **Convergence under refinement.** Wherever the diagnostic is below the tolerance, refining the
   wind nodes, the mesh or the sub-isobars changes the delivered temperature by less than the
   tolerance.
4. **The flag is honest.** Wherever it is above, the level is flagged, on every run.

## 7. What changes, and what does not

**Unchanged:** the reduction chain, the closure production, the placement of anchors, the mesh and
its columns, the kernel, the tracing, the outer loop and the reading of the wind.

**Changed:**

- `produce_at_target`: the target's state from the traced heights and the labels, not from a
  production of `ln N`.
- The transfer of `ln N`: kept, as the diagnostic.
- Kind P: layer variables added (the thickness, the layer temperature, the diagnostic and the flag);
  the levels' `p` becomes the label.
- F6, F7 and F8.
- The registered transfer product rebuilt. The closure product and the reduction chain are not.

**The estimate, the largest consequence.** A29 to A35 combine the anchors' `ln N` on the isobars at
the gauge latitude, and `D_ij` is a difference of `ln N` there. In this arrangement the natural
quantities are each anchor's `Phi_k` at the gauge (per isobar) or its thicknesses (per layer). This
is a change of formulation, not of code alone, and it reaches the manuscript (A27 to A35 and B7 at
a transferred latitude). It is the author's.

**Accepted checks restated or retired:**

- the identity checks of SPEC_04 Steps 4 and 5;
- `step05_4` check 4;
- SPEC_06 checks 6 and 7;
- SPEC_07's check 7.

Most become exact, or measure the diagnostic instead.

## 8. For the author

1. **The arrangement.** Charney and Phillips on the isobars as in section 4, the horizontal as it
   stands?
2. **The level rule.** A or B (B recommended)?
3. **The estimate.** Anchors combined in `Phi_k` per isobar at the gauge, or in thickness per
   layer? This decides how the manuscript's A29 to A35 are restated.
4. **The order.**
   - My suggestion: this as the next specification, ahead of any commit of SPEC_07, since it
     changes how much the interpolant matters.
   - SPEC_07's L2 is then judged on the general test of section 6 and committed or dropped on that
     evidence.
   - SPEC_08 becomes the sub-isobar refinement inside this design.
5. **The tolerance** for the flag: one physical number for every run.

Before any of this becomes a specification, the coding agent's reading would confirm the claims of
sections 2 and 3 against the code: in particular that the identity is the disagreement of the two
discretizations as stated, which a short measurement can show on run 6 and on a dense run.
