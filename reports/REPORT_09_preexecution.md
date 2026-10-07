# REPORT 09, pre-execution. The coding agent's reading of SPEC_09 v0.3

Coding agent, 7 October 2026. `main` at `2fbcf1d`.

**Nothing blocks Step 1.** Fig. 1 is on page 2 of `docs/conrath_etal_IRIS_1983.pdf` (journal page
287). It is a clean line scan:
- one temperature axis, 80 to 120 K, ticked every 10 K;
- latitude from 90 on the left to -90 on the right, ticked every 30 deg;
- the three levels well apart in temperature.

At 600 dpi a dot is about 11 px across and a stroke of the broken curve about 6 px thick, so an
opening wider than a stroke and narrower than a dot removes the curve and keeps the dots. The
method is the one `data_static/winds/digitize_smith1982_fig4.py` used for Smith et al. (1982):
threshold, open, label, and split each merged component into `round(area / median single-dot
area)` dots by k-means.

Points stated, none needing a ruling:

1. **The renderer.** No Python PDF renderer is installed, and none is a dependency. The script calls
   Ghostscript (`gswin64c` or `gs`, found on the path, overridable by argument). The output
   depends on it only through the rendered page.
2. **The PDF is not committed.** It is untracked in `docs/`. The script takes its path, and the CSV
   is the committed product. Committing the paper is the author's choice.
3. **Decisions by hand, coded** (deliverable 1):
   - the figure's crop on the page;
   - the three label boxes ("730 mbar", "290 mbar", "150 mbar") blanked;
   - the boundaries that assign a dot to its level. The boundary between 730 and 290 is 105 K. The
     boundary between 290 and 150 is a piecewise-linear line in latitude, because the two levels
     come within about 4 K near 75 N (a 290 mbar dot at about 90 K) and near -90 (a 150 mbar dot
     at about 90 K).
4. **The broken curve** touches 150 mbar dots between about 0 and -15 deg. The opening separates
   them, and the overlay shows the result.
