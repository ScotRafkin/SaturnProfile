# SPEC 20. Latitude-pressure sections, and the Voyager profiles as anchors

Version 0.5, 9 November 2026. Author of record: S. Rafkin. **Status: closed. Step 1 accepted
(REVIEW_20_step1), with heat capacity and SPEC_15's gravity fix at the acceptance commit.**
Depends on SPEC_15, SPEC_17, SPEC_18 and SPEC_19.

## 0. Purpose

Two looks at the model, with no fitting and no tuning.

1. **The latitude-pressure section.** Transfer an anchor to every latitude and draw the result as a
   section. It shows at once whether the delivered atmosphere looks right: the tropopause, the
   stratosphere, the stability, and what happens away from the data.
2. **Every Voyager profile as an anchor.** The design note (§6) says no radio profile is special.
   Here the three other Voyager profiles become anchors beside Lindal's, and each is transferred to
   every latitude and to every other profile. In a perfect world the sections would be identical;
   where they differ shows what the wind and the model leave unexplained. One case is wanted on its
   own: Voyager 1's egress (3.0 S) transferred to its own ingress (74.2 S), the same occultation,
   hours apart, with no Lindal anchor and no nine months between.

Lean: one new tool for the section, the existing chains for the anchors, figures. No new physics,
no wind change. Protocol as SPEC_04 §0. No regression: Lindal's files are untouched, and his answer
must come out the same.

## 1. Step 0: pre-execution report

Before building, a short report (`reports/REPORT_20_preexecution.md`) proposing:
1. **The anchors' chain.** The three profiles are traced temperature against pressure (SPEC_15),
   on Lindal's 94/6 composition. The proposal follows SPEC_18's three layers: source-specific
   ingest, then the general tools from the standard kinds on. Lindal's reduction set (kinds G, R and
   W) is used as is. Heights: the Voyager 2 egress carries Lindal's altitudes; for the two Voyager 1
   profiles, the proposal says how the radius is fixed (Lindal's Fig. 9 100 mbar point at each,
   `occul_data/lindal/isobar_100mbar.csv`, with hydrostatic integration from it, is the obvious
   choice) and what that costs. **Ruled (REVIEW_20_step0):** heights hydrostatic from 100 mbar, with
   gravity at kind D's radius, which also sets the anchor's absolute radius; the Fig. 9 point is a
   check only. Composition Lindal's C with his ammonia; Lindal's G and R reused, each file's own date
   in its comment.
2. **The section.** Whether one march per hemisphere already carries the delivered state through
   every intermediate latitude, so the section is a byproduct; otherwise one transfer per latitude.
   Where the section's product and control file live. **Ruled:** two marches per anchor, the
   production redone at each latitude, 87.5 S to 87.5 N; `casspian-section` in
   `src/casspian/tools/section/`; `sections/<anchor>_<wind>/`, control file committed, product
   ignored.
3. **Anything that blocks.**

## 2. Step 1

**Deliverables.**
1. **Three new anchors, kind N, with a closure run each:** Voyager 2 egress (31.2 S), Voyager 1
   ingress (74.2 S), Voyager 1 egress (3.0 S), under `occul_data/lindal/` beside Lindal's files,
   none of which changes. Epochs and the reliability notes of Lindal (scintillation, C-ring
   attenuation, lost records, the ray's sweep from 1.5 S to 8.9 S at egress) go in the file comments.
2. **The section tool** (the command named in the pre-execution report): an anchor, a wind and a
   latitude grid in; temperature, density, the wind, and the lapse rate on that grid out, 1.3 bar
   to 0.2 mbar. Grid: every 2.5 deg planetocentric, pole to pole, as far as the run's wind extends.
3. **The sections:** each of the four anchors under hypothesis 3 (`voyager_fit_1`), and the Lindal
   anchor under hypothesis 1.
4. **The anchor matrix:** each of the four anchors transferred to the latitude of each of the other
   three, under hypothesis 3.

**Acceptance (`tests/step20_1/accept_step20_1.py`, under `reports/step20_1/`).**
1. **Lindal's answer is unchanged.** No registered Lindal file changes; any shared tool that changes
   rebuilds his files in a copy, equal in value (SPEC_18's rule).
2. **Each new anchor's closure:** delivered at its own latitude against its own profile, within the
   tracing's precision, with the pressure identity.
3. **The section against single runs:** the Lindal section under hypothesis 3, at 31.2 S, 74.2 S and
   3.0 S, against the delivered profiles of the SPEC_19 runs. They should agree to about 0.1 K; any
   larger difference is explained.
4. **For the author, absolute first:**
   - **The anchor matrix as absolute profiles:** at each profile's latitude, the observed profile
     with its sigma_T band and the profiles delivered there from each of the other anchors.
     Then a table of delivered minus observed at 730, 290, 110, 30, 10 and 1 mbar, every pair.
     **Voyager 1 egress to Voyager 1 ingress** is drawn and tabled on its own.
   - **Reference adiabats on every temperature profile panel (author, 9 November):** the dry
     adiabat `d ln T / d ln p = R / c_p`, with `c_p(T)` computed from the run's composition
     (REVIEW_20_step1 ruling 1), for normal and for equilibrium hydrogen, both drawn and labeled,
     started from the observed profile's deepest point and drawn up to 100 mbar. They show how
     close the deep profile is to adiabatic; Lindal found it about 6 percent below near 1 bar.
   - **The temperature sections,** one per anchor and the hypothesis 1 section, isotherms on log
     pressure, the same color scale throughout. Overlaid: the zonal wind of the section's wind as
     dashed contours labeled in m/s (author, 9 November); the IRIS points at their three levels in
     the same colors; the latitudes of the four Voyager profiles and of Pioneer 11 marked. No
     hatching (author, 9 November).
   - **The sections' differences** from the Lindal anchor's, one per other anchor, on one diverging
     scale.
   - **The zonal wind section** of hypothesis 3, and its posterior sd.
   - **The lapse rate** against the same adiabat, per anchor; any
     superadiabatic region flagged and measured.
   - **Density** at 1 bar, 100 mbar and 10 mbar against latitude, one line per anchor.

## 3. Notes (accepted with the specification)

1. **Not independent.** Hypothesis 3 was fitted to these profiles as targets of the Lindal anchor, so
   the anchor matrix is partly constrained by construction. Where it still disagrees is what the
   fit could not absorb. The Voyager 1 pair is the least constrained case.
2. **Expected, not a fault.** The 3.0 S anchor carries its 10 to 20 mbar inversion and its layering to
   every latitude, and its section will differ most above 30 mbar. Lindal rates both Voyager 1
   profiles less reliable. These are looks, not tests with bounds.
3. **Pioneer 11** is marked on the sections and enters nothing (epoch, SPEC_17).
4. **No tuning on what the sections show.** Findings go to the status note; changes to the wind
   wait for new data (author, 9 November).

## 4. Revision history

- v0.1, 9 November 2026: first draft.
- v0.2, 9 November 2026: accepted.
- v0.3, 9 November 2026: Step 0 accepted; its rulings recorded in §1.
- v0.4, 9 November 2026: a reference adiabat on every temperature profile panel.
- v0.5, 9 November 2026: closed; heat capacity from composition, both hydrogen adiabats; SPEC_15's
  gravity fixed at the source; wind contours on the sections, hatching removed.
