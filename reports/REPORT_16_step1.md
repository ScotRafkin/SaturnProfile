# REPORT 16, Step 1. Before acceptance

Coding agent, 30 October 2026. Specification: `docs/specs/SPEC_16_Kliore_Refractivity.md` v0.2.
Working tree on `main` at `0128d65`, uncommitted. No regression (SPEC_16 §0); no change to the
model, the tools, a data file or a control file.

Acceptance: `tests/step16_1/accept_step16_1.py` **2 of 2**; output in `reports/step16_1/output.txt`,
figures, `self_check.csv` and `summary.json` beside it. The two SPEC_15 runs at 11.9 S (9.51 S
planetocentric) are read as they are.

**In short.**
- **The self-check agrees.** Kliore's T and p, as published, give his refractivity column through
  `reduction.refractivity` (RMS 0.034 n-units). With the `lindal1985` H2 and He values, the implied
  mean refractivity is 125.87 n-units at STP against our 125.90, so there is no systematic offset.
  One level misses the rounding, and it is a misprint in Kliore's table: 17.1 mbar should read
  about 17.5 (§1).
- **Against radius, with no composition,** the delivered refractivity matches Kliore's after a
  radius offset of +9.2 km on the planetocentric reading. From 10 to 182 mbar the residual is 2.3
  percent RMS. Above 10 mbar it reaches 10 to 30 percent, the same wave structure SPEC_15 found in
  temperature (§2).
- **The two readings of the radius mapping** fit equally well, because the offset absorbs the
  frame. The planetocentric reading needs only +9.2 km; the planetographic one needs +96 km.
- **Against pressure,** the comparison restates SPEC_15's temperature comparison exactly. The
  rescaling leaves Kliore's refractivity unchanged, so at equal pressure the refractivity ratio is
  the inverse temperature ratio. Delivered over observed: 0.97 at 100 mbar, 0.84 (hypothesis 1)
  and 0.82 (hypothesis 2) at 10 mbar (§3).

## 1. The self-check (deliverable 1, acceptance 1)

Kliore's temperature and pressure as published, for 90 percent H2 and 10 percent He by number,
were run through `reduction.number_density` and `reduction.refractivity`. The mean refractivity
per molecule was converted as `tools/composition/build_composition.py` converts it (e6 x 1e-6 /
n_Loschmidt). The species values are the `lindal1985` set: H2 136.0 and He 35.0, a mean of 125.90
n-units at STP.

- 128 levels: RMS 0.034 n-units, largest -0.154 n-units at 17.1 mbar. In percent: RMS 0.74, largest
  -3.50 at 2.7 mbar.
- **The largest percentages are rounding.** Above 3 mbar the column is 0.5 to 0.9 n-units, so the
  table's 0.05 n-unit half step is 6 to 10 percent there. Ours lies within it at every such level.
- **No systematic offset.** The least-squares ratio of Kliore's column to ours, from 10 n-units
  down, is 0.99980. That implies a mean refractivity of 125.87 n-units against our 125.90, a
  difference of 0.02 percent. Kliore's species refractivities are the ones Lindal used, or agree
  with them at this mixture to that level. The mean of ours minus his is +0.005 n-units.
- **One level beyond rounding, 17.1 mbar, is a misprint in the original.**
  - That row reads 6.2 n-units, 96.0 K and 17.1 mbar. The scan (journal page 5864) prints 17.1, so
    the transcription is faithful.
  - The neighbors run 15.5, 16.5, [17.1], 18.4 and 19.3 mbar, and 17.1 breaks the spacing.
  - The pressure that gives 6.2 n-units at 96.0 K is 17.54 mbar. Interpolating in ln p against
    refractivity between 16.5 and 18.4 mbar also gives 17.55.
  - Without this row, the largest difference is 0.086 n-units, inside the rounding.
  - The data file is left as it is: this step changes no data file. **For the author:** should the
    row be corrected to 17.5, with an inline comment as for SPEC_15's two misprints? It is one
    level of 128, and it moves no result in this report by more than its last digit.

Every level is in `reports/step16_1/self_check.csv`.

## 2. Refractivity against radius (deliverables 2 and 3)

**Kliore's radius, back to the local radius.** His column is the equatorial radius from his
measured radius, mapped with a constant oblateness, f = 0.096:
- p. 5861: "a constant oblateness value based on an ellipsoidal approximation";
- p. 5857: 0.096 at the 100 mbar level.

Undone with r_local = r_eq (1 - f) / sqrt((1 - f)^2 cos^2 phi + sin^2 phi), there are two readings:
- **planetocentric:** phi = 9.5 S, his table's latitude;
- **planetographic:** the same latitude converted on his own ellipse, 11.57 S. Lindal's label of
  11.9 S comes from Lindal's figure of Saturn, not Kliore's.

The two readings differ by 86.8 km in r_local, close to the specification's "about 80 km". At 99.6
mbar, Kliore's equatorial 60,333.0 km becomes:
- 60,150.0 km (planetocentric) or 60,063.3 km (planetographic);
- against 60,157.0 km delivered by both runs.

**The fit.** The offset is delivered radius minus Kliore's r_local at equal refractivity. It is
found by least squares in ln N over all 128 of Kliore's levels, unweighted as specified. It is
required that every level lie inside the delivered column.

| run | reading | offset (km) | RMS in ln N before (%) | RMS after (%) |
|---|---|---:|---:|---:|
| hypothesis 1 (case ii) | planetocentric | +9.20 | 23.7 | 4.58 |
| hypothesis 1 (case ii) | planetographic | +96.00 | 213.7 | 4.64 |
| hypothesis 2 (case iii) | planetocentric | +9.37 | 24.4 | 6.22 |
| hypothesis 2 (case iii) | planetographic | +96.16 | 213.6 | 6.32 |

**What the fit shows.**
- **The readings are indistinguishable in shape.** The residual after the offset is the same to
  0.1 percent. The 87 km difference is nearly a rigid shift over the 1.8 to 182 mbar span, so the
  fit cannot choose between them. Only the size of the offset does: +9 km against +96 km.
- **The offset is robust.** Fitted on 10 to 182 mbar only, it is +9.05 km (hypothesis 1) and +8.86
  km (hypothesis 2), against +9.20 and +9.37 on all levels. The upper levels do not drive it.
- **By range, residual in ln N x 100, planetocentric reading:**

  | pressure (mbar) | levels | hypothesis 1 RMS (mean) | hypothesis 2 RMS (mean) |
  |---|---:|---:|---:|
  | 10 to 182 | 112 | 2.31 (-0.48) | 2.23 (-1.37) |
  | 2.5 to 10 | 13 | 8.97 (+1.31) | 12.99 (+8.44) |
  | 1.8 to 2.5 | 3 | 18.7 (+18.2) | 27.1 (+26.8) |

  Below 10 mbar, the delivered refractivity follows Kliore's to about 2 percent once the frame is
  aligned. That is to first order what the author's working assumption, refractivity invariant,
  expects.
- **Above 10 mbar the departure is structure, not offset.** Delivered refractivity is low by up to
  10 percent near 60,250 km (about 9 mbar), then high by 10 to 30 percent above 60,280 km (above
  about 4 mbar). This is the 10 to 30 K warm and cold layering SPEC_15 found in Kliore's
  temperatures at 11.9 S, seen in the quantity the occultation measures. At the top levels the
  table's rounding contributes up to 10 percent.

**The offset combines** the radius mapping (Kliore's constant f against the model's geoid at 9.51 S),
the Pioneer 11 orbit, and the model's reference surface. This step does not separate them.

Figure: `reports/step16_1/refractivity_radius.png`.
- Left: ln N against radius, both runs and both readings, before (open) and after (filled) the
  offset fitted to hypothesis 1.
- Right: the residual after each fit.

## 3. Refractivity against pressure (deliverable 4)

Kliore's pressure was rescaled to the anchor's mixture by SPEC_15's rule (x0.9341); his refractivity
was not. Delivered N, interpolated in ln p, was compared with his N at the rescaled pressures.

| run | RMS (%) | mean (%) | delivered over observed at 10 / 30 / 100 / 150 mbar |
|---|---:|---:|---|
| hypothesis 1 (case ii) | 7.65 | -4.32 | 0.8428 / 1.0049 / 0.9714 / 0.9886 |
| hypothesis 2 (case iii) | 8.16 | -5.81 | 0.8180 / 0.9636 / 0.9654 / 0.9890 |

**This comparison is SPEC_15's temperature comparison, not a new one.**
- Lindal's rule rescales T by m and p by m / r, where m is the ratio of mean molecular masses and r
  the ratio of mean refractivities.
- The rescaled profile therefore has p_r R_94 / (k T_r) = p R_90 / (k T), Kliore's own N, exactly.
- So at a rescaled pressure, delivered N over Kliore's N is T_rescaled over T_delivered, the inverse
  of SPEC_15's ratio.

For example, at 10 mbar, 0.84 is the +16.9 K of SPEC_15's table on about 90 K. The pressure columns
also do not depend on the reading, which is why §2 ruling 1 calls this the comparison independent
of the frame.

Figure: `reports/step16_1/refractivity_pressure.png`.
- Left: N against pressure, the two runs and Kliore rescaled, with his pressure ranges rescaled.
- Right: delivered over observed minus one.

## 4. For the next specification (§2 ruling 2)

The self-check agrees and the comparison is sensible, so ruling 2's condition is met. Three points
bear on reducing Kliore's table to a kind N:
- **Radius.** Kliore's radius column is equatorial. A kind N needs the local radius, so the reducer
  needs the mapping, f and the reading.
  - The planetocentric reading is the one whose offset is small (9 km).
  - The planetographic reading would put Pioneer 11's levels about 87 km lower, against the same
    model.
  - Choosing it is the author's ruling. This step can only say that the shape does not decide it.
- **The 17.1 mbar row** (§1).
- **Kliore's error ranges** are extremes, not one sigma. As validation at weight 0 they enter no
  fit. They will matter when the next step assimilates.

## 5. Files

New, uncommitted:
- `tests/step16_1/accept_step16_1.py`;
- `reports/step16_1/`:
  - `output.txt`, `summary.json` and `self_check.csv`;
  - `refractivity_radius.png` and `refractivity_pressure.png`;
- this report.

Untracked from the author: `docs/specs/SPEC_16_Kliore_Refractivity.md`; modified,
`docs/specs/STATE.md`.
