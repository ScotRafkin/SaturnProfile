# Transcription notes, Lindal et al. (1985) Voyager 2 ingress

Documentation only. Not read by any tool; its hash is recorded in the raw bundle's
`raw_sources` beside the two data files.

## lindal_table1.csv

Transcribed from Table I of Lindal, Sweetnam and Eshleman (1985), AJ 90, 1136, during Phase 1
(August 2026). Sixty-six rows, columns `pressure_mbar`, `temperature_K`, `nh3_ppm` (blank where
Table I has no entry), `altitude_km` (above the source's 1 bar level). Checked in Phase 1:
66 rows; the 1 bar row exact at 1000.00 mbar, 134.8 K, 0.0 km; nine finite NH3 entries from
831.76 mbar (2.6 ppm) to 1258.93 mbar (66.9 ppm), with the interior gap at 1047.13 mbar and
the end gap at 1298.48 mbar left blank. The profile values are trusted (S. Rafkin, 10 September
2026). Units are converted by the stage-one tool, never in this file.

## lindal_scalars.toml

Every scalar the paper states that the reduction needs, each with the page, table, figure, or
caption it was read from. Read from the PDF in `docs/` on 9 and 10 September 2026 by the app
agent and accepted by S. Rafkin on 10 September 2026 (REVIEW_01_step2). If a value is found
wrong, correct it here and rerun the chain from stage one; nothing downstream is edited by hand.

Two items worth knowing when checking it. The paper states no numerical top boundary condition
for its hydrostatic integration (only "from the top of the detectable atmosphere and
downward", p. 1138); the `[top_boundary]` entry records that as a statement, not a number. The
Fig. 9 caption names a third wind source, Garneau (1984), beyond the two named in the Appendix.

## The pressure grid

Table I prints its pressures to two decimals of a millibar, and every row from 0.20 to 1258.93
mbar lies on the grid `10^(k/100)` mbar within that precision. Of the 65 rows above the last, 62
admit exactly one grid value (48 of them printed to four, five or six figures, 10 to three and 4
to two). Three rows printed to two figures admit two: 0.20 (`k` = −70 or −69), 0.25 (−61 or
−60) and 0.32 (−50 or −49). They are adjacent, so they are resolved one at a time from below by
continuing the spacing of the two nearest fixed rows: 0.32 takes −50 from 0.40 and 0.50, then
0.25 takes −60, then 0.20 takes −70. The spacing along the table is then 10 hundredths of a
decade to 20 mbar, and 8, 6, 4 and 2 below that. The last row, 1298.48 mbar, is the end of the
data rather than a grid level and stays as printed. This is an inference by the reviewing agent
(14 September 2026), not a statement of the source, adopted by the author; it is declared in the
`[pressure_grid]` table of `lindal_scalars.toml`, and the CSV is not changed. Both the printed and
the grid values are carried in the raw bundle and in kind T. If the source's level list is ever
found, it replaces the inference.

## The anchor latitude

Carried as 36.5 degrees through the project record until 9 September 2026, when the Fig. 4
label was read directly: 36.3 degrees N planetographic, the latitude at which the lowest point
on the ray touched the 100 mbar surface (Fig. 4 caption). The swath was 36.3 to 36.7 degrees
(p. 1138). The superseded value survives only in the Phase 1 netCDF preserved under the
`phase1-prototype` tag.
