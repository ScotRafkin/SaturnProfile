# REVIEW 20, Step 1. Sections, and the Voyager profiles as anchors

Reviewing agent, 9 November 2026. Report: `reports/REPORT_20_step1.md`. Specification:
`docs/specs/SPEC_20_Sections_and_Anchors.md` v0.4, closed at v0.5.

**Accepted, with two items at the acceptance commit. SPEC_20 closes.**

- **The build is clean.** Lindal's files are byte-identical; the three new anchors close on their own
  profiles within the tracing's precision; the section reproduces SPEC_19's runs to 6.4e-3 K. No
  shared tool or model code changed. Departures 1 to 7 of §2 are accepted.
- **What the looks show.** Each anchor's offset from Lindal's is nearly uniform in latitude: a
  single-anchor transfer carries everything its anchor holds, including what no zonal wind can
  produce (the 3.0 S inversion appears at every latitude). This is expected of one anchor. The
  product blends all anchors, so an observation dominates where and when it was made (author).
- **The weighting favors IRIS.** At 110 mbar the two southern profiles miss each other by 6.5 to
  6.8 K, about 3 sigma_T, while the radio group ratios are 0.2 to 0.3 and IRIS's is 2.3. IRIS's errors
  are random only and too small, and its top level is labeled 110 mbar by Lindal and 150 by Conrath
  and Pirraglia. The IRIS systematic term deferred in SPEC_19 is worth 3 to 7 K here; that deferral
  was wrong. It goes into SPEC_21.
- **The north pole** under hypothesis 3 (63 K at 110 mbar, 87.5 N, poleward of IRIS's last row) is
  doubtful; a temperature-based prior would bound it. For the status note.

**Rulings (author, 9 November):**
1. **Heat capacity from composition.** No run specifies it. The mixture's `c_p(T)` comes from the
   run's composition, with each species' `c_p(T)` from physics: H2 by its rotational partition
   function (rotational temperature about 85 K), both normal (ortho to para 3:1, frozen) and
   equilibrium; He 5/2 R; minor species as rigid rotors, immaterial at their abundance. A table in
   `data_static`, computed once with its basis in a comment, is equally acceptable. Every temperature
   panel and the lapse-rate figure draw both adiabats, labeled.
2. **SPEC_15's gravity:** an error in our code, so fixed at the source though under materiality:
   `egress_profile.to_pressure` evaluates gravity at the geoid radius. Rebuild
   `voyager2_egress.csv` and the Voyager 2 egress anchor from it. No refit now; the next fit uses it.
3. **SPEC_21's scope** widens (status note, Next).
4. **The section figure (author):** zonal wind of the same wind overlaid on each temperature section
   as dashed contours, labeled in m/s; the shear hatching removed from every figure, as it gives the
   eye nothing usable. Standing for all section figures.
5. **The final single answer** blending all anchors is SPEC_22's product, with its own profile matrix
   and sections.

**Order of work:** the two items and the section figure, the products rebuilt clean, the adiabat figures and the rebuilt
egress closure reported briefly, then commit and push. No regression: Lindal's files untouched.
