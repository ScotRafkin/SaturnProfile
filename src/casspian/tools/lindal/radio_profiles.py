"""Three more radio profiles: Lindal et al. (1985) Figs. 5 and 6, traced, and Kliore et al. (1980)
Table 1. SPEC_15 Step 1.

**Lindal et al. (1985), AJ 90, 1136**, page 5 of the PDF (journal page 1140). Fig. 5: temperature
against pressure for 36.3 N (solid, the anchor, Voyager 2 ingress), 31.2 S (dashed, Voyager 2
egress) and 74.2 S (dash-dot, Voyager 1 ingress, 13 November 1980). Fig. 6: 3.0 S (solid, Voyager 1
egress, 13 November 1980) and 11.9 S (stippled, Pioneer 11 exit, 1 September 1979, redrawn from
Kliore et al.). Latitudes as Lindal labels them, planetographic (SPEC_15 section 2 ruling 1). Both
figures are computed for 94 percent hydrogen and 6 percent helium by number, 2.135 amu, so nothing
is rescaled.

**Tracing**, every step coded so that running this reproduces the CSVs:

1. Page 5 is rendered at 600 dpi by Ghostscript as a one-bit image (darker than mid-gray is ink).
2. Each figure's frame is found inside its `box` from the rows and columns inked along most of it.
   Pressure is logarithmic: a least-squares line of log10 p in page row through the ticks at 1, 10,
   100 and 1000 mbar and the frame lines at 0.1 and 10^4 mbar. Temperature is a least-squares line
   in page column through the ticks every 10 K and the frame lines at 80 and 160 K.
3. **Each curve is traced along a guide**, not scanned: Fig. 5's curves loop, curl, cross and
   touch, so no one scan direction separates them. The guide is a polyline of control points placed
   by hand on the scan (`CURVES`), close together in bends and sparse where the curve is smooth;
   these are the decisions made by hand, and they are coded. The guide is resampled every
   `STEP_PX` along its length, and each point is snapped across the guide's direction to the centre
   of the ink run it lies on, or of a run within `NEAR_PX` if it lies on none. Where there is none
   (a dash's gap: the neighbouring curve is not borrowed), or the run is wider than `MERGED_WIDTHS`
   line widths (two curves crossing or touching, a label's leader, a dot), the guide's point stands,
   so the guide is placed on the line's centre to within a few pixels.
4. The positions are converted to temperature and pressure, and reduced to the median of each
   `BIN_DEX` bin of log10 p: a profile is single-valued in pressure, and the bins remove the
   tracing's jitter. An empty bin is filled linearly in log p from its neighbours.
5. **Coincidence.** Where two curves of a figure lie within one line width (`LINE_WIDTH_PX`) of
   each other the scan cannot tell them apart; both traces then follow the shared ink, and
   `coincident_ranges` lists those pressure ranges.
6. **Error bars** (`BARS`): each horizontal bar is measured on its own row, from its end caps, and
   read as a temperature uncertainty of half its length at that row's pressure, assigned to the
   curve it sits on (SPEC_15 section 2 ruling 3). Elsewhere the uncertainty is NaN.

**Kliore et al. (1980), JGR 85, 5857, Table 1**, the Pioneer 11 exit (`read_pioneer11`), is
transcribed in `occul_data/kliore1980/pioneer11_exit.csv` as published, two misprints corrected
there with a comment each. It is computed for 90 percent hydrogen and 10 percent helium;
`rescale_composition` takes it to another mixture by Lindal's rule (Fig. 5 caption): temperature
in proportion to the mean molecular mass, pressure in proportion to the mean molecular mass over
the mean molecular refractivity.
"""

from __future__ import annotations

import argparse
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


PAGE = 5
DPI = 600
STEP_PX = 3.0
SNAP_PX = 9
#: A point on no ink snaps to a run only this close; farther, it is in a gap and stays.
NEAR_PX = 4
LINE_WIDTH_PX = 10.0
MERGED_WIDTHS = 1.8
BIN_DEX = 0.01
#: An error bar: the rows searched either side of its point, and the gaps joined along it.
BAR_ROWS_PX = 6
BAR_GAP_PX = 3

FIGURES = {
    "fig5": {"box": (400, 400, 2500, 2600), "pressure_ticks_mbar": (1.0, 10.0, 100.0, 1000.0),
             "frame_pressure_mbar": (0.1, 1.0e4), "frame_temperature_K": (80.0, 160.0),
             "temperature_tick_K": 10.0},
    "fig6": {"box": (2500, 400, 4700, 2500), "pressure_ticks_mbar": (1.0, 10.0, 100.0, 1000.0),
             "frame_pressure_mbar": (0.1, 1.0e4), "frame_temperature_K": (80.0, 160.0),
             "temperature_tick_K": 10.0},
}

#: Fig. 5, 36.3 N, solid: from the dot on its top bar, round the nose and down the hook; the
#: upper middle by rows; the upper branch by columns; the left side by rows; the lower branch
#: by columns, to the dot on its bottom bar.
_FIG5_SOLID = [
    # top: the bar's dot, the nose and the hook
    (1902, 743), (1950, 768), (2000, 782), (2045, 797), (2080, 818), (2098, 838), (2085, 852),
    (2055, 856), (2035, 862), (2030, 885), (2018, 905), (1997, 925), (1990, 950), (1990, 980),
    # upper middle, by rows
    (1997, 1000), (2008, 1020), (2008, 1040), (1981, 1060), (1954, 1080), (1930, 1100),
    (1922, 1120), (1907, 1140), (1886, 1160), (1834, 1180), (1796, 1200), (1752, 1220),
    (1691, 1240), (1638, 1260),
    # upper branch, by columns
    (1590, 1276), (1550, 1289), (1510, 1301), (1470, 1315), (1430, 1326), (1390, 1341),
    (1350, 1349), (1310, 1362), (1270, 1374), (1230, 1390), (1190, 1406), (1150, 1424),
    # left, by rows
    (1115, 1440), (1076, 1460), (1040, 1480), (1008, 1500), (980, 1520), (964, 1540), (957, 1560),
    (946, 1580), (936, 1600), (932, 1620), (934, 1640), (938, 1660), (947, 1680), (954, 1700),
    (965, 1720), (976, 1740), (994, 1760), (1012, 1780), (1033, 1800), (1060, 1820), (1088, 1840),
    (1122, 1860), (1160, 1880),
    # lower branch, by columns
    (1190, 1893), (1230, 1904), (1270, 1918), (1310, 1932), (1350, 1944), (1390, 1955),
    (1430, 1964), (1470, 1974), (1510, 1982), (1550, 1991), (1590, 1998), (1630, 2007),
    (1670, 2015), (1710, 2024), (1750, 2031), (1790, 2039), (1830, 2048), (1870, 2056),
    (1910, 2062), (1950, 2071), (1990, 2080), (2030, 2087), (2062, 2092),
]

#: Fig. 5, 31.2 S, dashed: from its dot through the curls, then as the solid. From x 1200 to
#: 1350 on the lower branch it lies just under the solid, half merged with it, and from 1380
#: it joins the dash-dot's band before parting from it near 1500.
_FIG5_DASHED = [
    (1792, 762), (1812, 768), (1826, 778), (1828, 790), (1818, 797), (1806, 796), (1793, 805),
    (1808, 818), (1824, 828), (1838, 843), (1832, 860), (1816, 867), (1799, 876), (1815, 881),
    (1832, 882), (1845, 877), (1856, 885), (1845, 894), (1828, 899), (1811, 906), (1828, 914),
    (1843, 921), (1862, 937), (1880, 952), (1890, 970),
    (1897, 1000), (1891, 1020), (1882, 1040), (1862, 1060), (1848, 1080), (1842, 1100),
    (1828, 1120), (1813, 1140), (1798, 1160), (1780, 1180), (1758, 1200), (1726, 1220),
    (1691, 1240), (1638, 1260),
    (1590, 1276), (1550, 1289), (1510, 1301), (1470, 1315), (1430, 1326), (1390, 1341),
    (1350, 1359), (1310, 1374), (1270, 1391), (1230, 1408), (1190, 1426), (1164, 1440),
    (1121, 1460), (1080, 1480), (1050, 1500), (1014, 1520), (996, 1540), (978, 1560), (964, 1580),
    (957, 1600), (951, 1620), (946, 1640), (947, 1660), (947, 1680), (954, 1700), (965, 1720),
    (976, 1740), (994, 1760), (1012, 1780), (1033, 1800), (1060, 1820), (1088, 1840),
    (1122, 1860), (1164, 1880),
    (1190, 1896), (1230, 1912), (1270, 1926), (1310, 1939), (1350, 1952), (1380, 1965),
    (1430, 1978), (1470, 1990), (1510, 1996), (1550, 2005), (1590, 2012), (1630, 2021),
    (1670, 2030), (1710, 2036), (1750, 2042), (1790, 2050), (1830, 2058), (1870, 2066),
    (1910, 2074), (1950, 2084), (1990, 2090), (2020, 2098),
]

#: Fig. 5, 74.2 S, dash-dot: from the dot on its top bar, under the solid's hook, down the
#: right of the other two and round the left, to the dot on its bottom bar at 1040 mbar.
_FIG5_DASHDOT = [
    (1846, 700), (1848, 740), (1866, 770), (1882, 795), (1900, 822), (1925, 838), (1960, 851),
    (1993, 862), (2015, 870), (2042, 888), (2046, 905), (2035, 915), (2008, 928), (1958, 945),
    (1945, 965), (1938, 980),
    (1932, 1000), (1930, 1020), (1926, 1040), (1925, 1060), (1924, 1080), (1928, 1100),
    (1924, 1120), (1926, 1140), (1916, 1160), (1904, 1180), (1890, 1200), (1869, 1220),
    (1838, 1240), (1800, 1260), (1764, 1280), (1704, 1300),
    (1670, 1313), (1630, 1326), (1590, 1335), (1550, 1347), (1510, 1358), (1470, 1370),
    (1430, 1384), (1390, 1396), (1350, 1411), (1310, 1426), (1270, 1442), (1230, 1458),
    (1190, 1470), (1162, 1480), (1132, 1500), (1101, 1520), (1076, 1540), (1056, 1560),
    (1038, 1580), (1026, 1600), (1015, 1620), (1014, 1640), (1016, 1660), (1022, 1680),
    (1034, 1700), (1046, 1720), (1060, 1740), (1076, 1760), (1088, 1780), (1112, 1800),
    (1123, 1820), (1136, 1840), (1150, 1860), (1170, 1880),
    (1200, 1905), (1215, 1914), (1240, 1922), (1260, 1932), (1271, 1940), (1300, 1947),
    (1340, 1958), (1350, 1962), (1390, 1971),
    (1430, 1980), (1470, 1992), (1510, 2006), (1550, 2017), (1590, 2026), (1630, 2034),
    (1670, 2042), (1710, 2050), (1750, 2059), (1765, 2066),
]

#: Fig. 6, 3.0 S, solid: from its top bar, the S near 1 mbar, the long run left above 30
#: mbar, round the left and down the lower branch to its bottom bar.
_FIG6_SOLID = [
    # top: the bar, then left along the curve
    (3935, 857), (3910, 860), (3880, 867), (3850, 872), (3820, 876), (3790, 881), (3760, 886),
    (3730, 890), (3700, 895), (3670, 900), (3645, 906), (3633, 920), (3632, 940), (3645, 955),
    (3662, 965), (3678, 980), (3688, 1000), (3682, 1015), (3665, 1028), (3645, 1042),
    (3625, 1052), (3605, 1065), (3597, 1085), (3600, 1100), (3610, 1120), (3618, 1140),
    (3622, 1160), (3620, 1180), (3610, 1200), (3590, 1220), (3570, 1240), (3552, 1260),
    (3526, 1280), (3512, 1300), (3519, 1320), (3540, 1340), (3540, 1360),
    (3500, 1376), (3460, 1383), (3420, 1389), (3380, 1394), (3340, 1400), (3300, 1406),
    (3260, 1410), (3220, 1416), (3180, 1420), (3140, 1425), (3100, 1432), (3026, 1440),
    (2966, 1460), (2938, 1480), (2923, 1500), (2910, 1520), (2906, 1540), (2909, 1560),
    (2918, 1580), (2914, 1600), (2918, 1620), (2937, 1640), (2962, 1660), (2992, 1680),
    (3016, 1700), (3017, 1720), (3016, 1740), (3020, 1760), (3027, 1780), (3040, 1805),
    (3060, 1825), (3100, 1842),
    (3140, 1862), (3180, 1880), (3220, 1894), (3260, 1907), (3300, 1917), (3340, 1928),
    (3380, 1936), (3420, 1947), (3460, 1956), (3500, 1964), (3540, 1974), (3580, 1982),
    (3620, 1990), (3660, 1998), (3700, 2007), (3740, 2015), (3780, 2022), (3820, 2030),
    (3860, 2038), (3880, 2041),
]

#: Fig. 6, 11.9 S, stippled: the check trace (SPEC_15 section 1 deliverable 1), its spike to
#: 130 K near 1.3 mbar included.
_FIG6_STIPPLED = [
    (3507, 903), (3505, 925), (3520, 932), (3540, 932), (3560, 942), (3580, 946), (3600, 952),
    (3620, 958), (3640, 963), (3645, 985), (3650, 1003), (3670, 1015), (3690, 1020),
    (3720, 1027), (3745, 1030), (3768, 1033), (3790, 1041), (3772, 1051), (3745, 1060),
    (3725, 1063), (3705, 1066), (3683, 1076), (3668, 1086), (3657, 1098), (3668, 1110),
    (3660, 1120), (3640, 1130), (3615, 1137), (3580, 1150), (3555, 1157), (3530, 1163), (3500, 1171),
    (3460, 1180), (3420, 1188), (3380, 1195), (3340, 1201), (3300, 1209), (3260, 1218),
    (3220, 1228), (3180, 1238), (3140, 1252), (3122, 1260), (3103, 1280), (3098, 1300),
    (3097, 1320), (3098, 1340), (3102, 1360), (3135, 1380), (3130, 1400), (3114, 1420),
    (3100, 1440), (3075, 1460), (3052, 1480), (3028, 1500), (3010, 1520), (2994, 1540),
    (2984, 1560), (2976, 1580), (2973, 1600), (2975, 1620), (2980, 1640), (2995, 1660),
    (3016, 1680), (3050, 1695), (3064, 1705), (3072, 1720), (3074, 1740), (3075, 1760),
    (3076, 1780), (3096, 1800),
]

#: Each curve: its figure and its guide, page px (x, y), placed by hand on the scan (module
#: docstring, step 3).
CURVES = {
    "36.3N": {"figure": "fig5", "guide": _FIG5_SOLID},
    "31.2S": {"figure": "fig5", "guide": _FIG5_DASHED},
    "74.2S": {"figure": "fig5", "guide": _FIG5_DASHDOT},
    "3.0S": {"figure": "fig6", "guide": _FIG6_SOLID},
    "11.9S": {"figure": "fig6", "guide": _FIG6_STIPPLED},
}

#: The error bars: the curve each sits on (where its dot is), and a page point on the bar clear of
#: the curves (module docstring, step 6).
BARS = {
    "74.2S top": {"curve": "74.2S", "point": (1720, 697)},
    "74.2S bottom": {"curve": "74.2S", "point": (1700, 2063)},
    "36.3N top": {"curve": "36.3N", "point": (2040, 742)},
    "36.3N bottom": {"curve": "36.3N", "point": (2120, 2093)},
    "31.2S top": {"curve": "31.2S", "point": (1680, 761)},
    "31.2S bottom": {"curve": "31.2S", "point": (1980, 2112)},
    "3.0S top": {"curve": "3.0S", "point": (4050, 855)},
    "3.0S bottom": {"curve": "3.0S", "point": (3930, 2041)},
}

#: The two files written, their curves and their top comments (SPEC_15 section 1 deliverable 3).
OUTPUTS = {
    "74.2S": ("occul_data/lindal/voyager1_ingress.csv", [
        "Voyager 1 ingress, 74.2 S as published (planetographic), 13 November 1980.",
        "Lindal et al. (1985), AJ 90, 1136, Fig. 5 (dash-dot), traced by",
        "src/casspian/tools/lindal/radio_profiles.py (SPEC_15 Step 1). 94 percent H2, 6 percent He by",
        "number, as the figure was computed. temperature_uncertainty_K is half a drawn error bar,",
        "at the pressure where it is drawn; NaN elsewhere.",
    ]),
    "3.0S": ("occul_data/lindal/voyager1_egress.csv", [
        "Voyager 1 egress, 3.0 S as published (planetographic), 13 November 1980.",
        "Lindal et al. (1985), AJ 90, 1136, Fig. 6 (solid), traced by",
        "src/casspian/tools/lindal/radio_profiles.py (SPEC_15 Step 1). 94 percent H2, 6 percent He by",
        "number, as the figure was computed. temperature_uncertainty_K is half a drawn error bar,",
        "at the pressure where it is drawn; NaN elsewhere.",
    ]),
}


def render(pdf) -> np.ndarray:
    """Page `PAGE` at `DPI`, ink True."""
    import subprocess
    import tempfile

    from PIL import Image

    from casspian.tools.lindal.egress_profile import find_ghostscript

    with tempfile.TemporaryDirectory() as directory:
        out = Path(directory) / "page.png"
        subprocess.run([find_ghostscript(None), "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                        "-sDEVICE=pnggray", f"-r{DPI}", f"-dFirstPage={PAGE}", f"-dLastPage={PAGE}",
                        f"-sOutputFile={out}", str(pdf)], check=True)
        with Image.open(out) as image:
            return np.asarray(image.convert("L")) < 128


def _runs(mask_1d, offset=0):
    edges = np.flatnonzero(np.diff(np.concatenate([[0], mask_1d.astype(np.int8), [0]])))
    return [(int(a) + offset, int(b) + offset) for a, b in zip(edges[::2], edges[1::2])]


@dataclass
class Calibration:
    frame: dict                 #: top, bottom, left, right: the inner edge of each frame line
    log_p_fit: tuple            #: (slope, intercept) of log10 p (mbar) in page row
    t_fit: tuple                #: (slope, intercept) of T (K) in page column
    p_residual_dex: np.ndarray  #: the ticks' departure from the fit
    t_residual_K: np.ndarray

    def temperature(self, x):
        return np.polyval(self.t_fit, np.asarray(x, dtype="float64"))

    def pressure_mbar(self, y):
        return 10.0 ** np.polyval(self.log_p_fit, np.asarray(y, dtype="float64"))

    def column(self, T):
        a, b = self.t_fit
        return (np.asarray(T, dtype="float64") - b) / a

    def row(self, p_mbar):
        a, b = self.log_p_fit
        return (np.log10(np.asarray(p_mbar, dtype="float64")) - b) / a


def calibrate(ink, name) -> Calibration:
    """The frame and the two scales of one figure (module docstring, step 2)."""
    spec = FIGURES[name]
    x0, y0, x1, y1 = spec["box"]
    box = ink[y0:y1, x0:x1]
    rows = _runs(box.sum(axis=1) > 0.5 * box.shape[1], y0)
    cols = _runs(box.sum(axis=0) > 0.4 * box.shape[0], x0)
    if len(rows) != 2 or len(cols) != 2:
        raise ValueError(f"{name}: found {len(rows)} frame rows and {len(cols)} frame columns; "
                         "the figure has 2 and 2")
    centre = lambda run: 0.5 * (run[0] + run[1] - 1)
    top, bottom, left, right = rows[0][1], rows[1][0], cols[0][1], cols[1][0]
    strip = 30
    tick_rows = [centre(r) for r in _runs(ink[top:bottom, left:left + strip].mean(axis=1) > 0.6, top)]
    tick_rows = [r for r in tick_rows if top + 20 < r < bottom - 20]
    tick_cols = [centre(c) for c in _runs(ink[bottom - 15:bottom, left:right].mean(axis=0) > 0.6, left)]
    tick_cols = [c for c in tick_cols if left + 20 < c < right - 20]
    p_ticks = spec["pressure_ticks_mbar"]
    t_lo, t_hi = spec["frame_temperature_K"]
    t_ticks = np.arange(t_lo + spec["temperature_tick_K"], t_hi - 1e-9, spec["temperature_tick_K"])
    if len(tick_rows) != len(p_ticks) or len(tick_cols) != t_ticks.size:
        raise ValueError(f"{name}: found {len(tick_rows)} pressure and {len(tick_cols)} temperature "
                         f"ticks; the figure has {len(p_ticks)} and {t_ticks.size}")
    p_rows = [centre(rows[0])] + tick_rows + [centre(rows[1])]
    p_values = np.log10([spec["frame_pressure_mbar"][0], *p_ticks, spec["frame_pressure_mbar"][1]])
    t_cols = [centre(cols[0])] + tick_cols + [centre(cols[1])]
    t_values = np.concatenate([[t_lo], t_ticks, [t_hi]])
    p_fit = np.polyfit(p_rows, p_values, 1)
    t_fit = np.polyfit(t_cols, t_values, 1)
    return Calibration(frame={"top": top, "bottom": bottom, "left": left, "right": right},
                       log_p_fit=tuple(p_fit), t_fit=tuple(t_fit),
                       p_residual_dex=np.polyval(p_fit, p_rows) - p_values,
                       t_residual_K=np.polyval(t_fit, t_cols) - t_values)


# ---------------------------------------------------------------------------
# The tracer
# ---------------------------------------------------------------------------

def resample(guide, step=STEP_PX):
    """Points every `step` px along the polyline `guide`, and the polyline's direction at each."""
    g = np.asarray(guide, dtype="float64")
    seg = np.diff(g, axis=0)
    length = np.hypot(seg[:, 0], seg[:, 1])
    s_nodes = np.concatenate([[0.0], np.cumsum(length)])
    s_out = np.arange(0.0, s_nodes[-1] + 1e-9, step)
    x = np.interp(s_out, s_nodes, g[:, 0])
    y = np.interp(s_out, s_nodes, g[:, 1])
    k = np.clip(np.searchsorted(s_nodes, s_out, side="right") - 1, 0, len(seg) - 1)
    angle = np.arctan2(seg[k, 1], seg[k, 0])
    return np.column_stack([x, y]), angle


def snap(ink, point, angle):
    """`point` moved across `angle` to the centre of its ink run (module docstring, step 3), and
    whether it moved."""
    x, y = point
    nx, ny = -np.sin(angle), np.cos(angle)
    offsets = np.arange(-SNAP_PX - 2 * LINE_WIDTH_PX, SNAP_PX + 2 * LINE_WIDTH_PX + 1)
    xs = np.clip(np.round(x + offsets * nx).astype(int), 0, ink.shape[1] - 1)
    ys = np.clip(np.round(y + offsets * ny).astype(int), 0, ink.shape[0] - 1)
    on = ink[ys, xs]
    runs = [(offsets[a], offsets[b - 1]) for a, b in _runs(on)]
    candidates = [(lo, hi) for lo, hi in runs if lo <= NEAR_PX and hi >= -NEAR_PX]
    if not candidates:
        return (x, y), False
    lo, hi = min(candidates, key=lambda r: 0.0 if r[0] <= 0 <= r[1] else min(abs(r[0]), abs(r[1])))
    if hi - lo + 1 > MERGED_WIDTHS * LINE_WIDTH_PX:
        return (x, y), False
    shift = 0.5 * (lo + hi)
    return (x + shift * nx, y + shift * ny), True


def trace_guide(ink, guide):
    """The traced path along a guide, page px, and which points the ink placed."""
    points, angles = resample(guide)
    out, snapped = [], []
    for point, angle in zip(points, angles):
        q, moved = snap(ink, point, angle)
        out.append(q)
        snapped.append(moved)
    return np.array(out), np.array(snapped)


@dataclass
class Trace:
    name: str
    figure: str
    path_px: np.ndarray                  #: (n, 2) page x, y in path order
    pressure_mbar: np.ndarray            #: on `BIN_DEX` bins, increasing
    temperature_K: np.ndarray
    filled: np.ndarray                   #: bins filled across a gap
    bars: list = field(default_factory=list)   #: (pressure mbar, half length K, centre K)


def binned(path_px, calibration: Calibration):
    """Median temperature in each `BIN_DEX` bin of log10 p (module docstring, step 4)."""
    T = calibration.temperature(path_px[:, 0])
    log_p = np.log10(calibration.pressure_mbar(path_px[:, 1]))
    lo, hi = np.ceil(log_p.min() / BIN_DEX), np.floor(log_p.max() / BIN_DEX)
    grid = np.arange(lo, hi + 1) * BIN_DEX
    index = np.round(log_p / BIN_DEX).astype(int) - int(lo)
    ok = (index >= 0) & (index < grid.size)
    value = np.full(grid.size, np.nan)
    for i in np.unique(index[ok]):
        value[i] = np.median(T[ok][index[ok] == i])
    filled = ~np.isfinite(value)
    value[filled] = np.interp(grid[filled], grid[~filled], value[~filled])
    return 10.0 ** grid, value, filled


def measure_bar(ink, point, calibration: Calibration):
    """A horizontal error bar through `point`, as (pressure mbar, half length K, centre K).

    The bar is a line two pixels thick between two end caps several rows tall, broken in places by
    a pixel or two, so on each row within `BAR_ROWS_PX` of `point` the runs `BAR_GAP_PX` or less
    apart are joined, and the longest joined run through `point` is the bar, end cap to end cap."""
    x, y = map(int, point)
    left, right = calibration.frame["left"], calibration.frame["right"]
    best = None
    for row in range(y - BAR_ROWS_PX, y + BAR_ROWS_PX + 1):
        joined = []
        for a, b in _runs(ink[row, left:right], left):
            if joined and a - joined[-1][1] <= BAR_GAP_PX:
                joined[-1] = (joined[-1][0], b)
            else:
                joined.append((a, b))
        for a, b in joined:
            if a <= x < b and (best is None or b - a > best[2] - best[1]):
                best = (row, a, b)
    if best is None or best[2] - best[1] < 60:
        raise ValueError(f"no horizontal bar through {point}")
    row, a, b = best
    T_lo, T_hi = calibration.temperature([a, b - 1])
    return (float(calibration.pressure_mbar(row)), 0.5 * float(T_hi - T_lo), 0.5 * float(T_hi + T_lo))


def trace_all(ink):
    calibrations = {name: calibrate(ink, name) for name in FIGURES}
    traces = {}
    for name, spec in CURVES.items():
        cal = calibrations[spec["figure"]]
        path, _ = trace_guide(ink, spec["guide"])
        p, T, filled = binned(path, cal)
        traces[name] = Trace(name=name, figure=spec["figure"], path_px=path, pressure_mbar=p,
                             temperature_K=T, filled=filled)
    for bar_name, bar in BARS.items():
        trace = traces[bar["curve"]]
        trace.bars.append(measure_bar(ink, bar["point"], calibrations[trace.figure]))
    return calibrations, traces


def coincident_ranges(a: Trace, b: Trace, calibration: Calibration):
    """Pressure ranges (mbar) where the two traces lie within `LINE_WIDTH_PX` of each other."""
    lo = max(a.pressure_mbar.min(), b.pressure_mbar.min())
    hi = min(a.pressure_mbar.max(), b.pressure_mbar.max())
    if lo >= hi:
        return []
    p = 10.0 ** np.arange(np.log10(lo), np.log10(hi), BIN_DEX)
    Ta = np.interp(np.log10(p), np.log10(a.pressure_mbar), a.temperature_K)
    Tb = np.interp(np.log10(p), np.log10(b.pressure_mbar), b.temperature_K)
    close = np.abs(calibration.column(Ta) - calibration.column(Tb)) <= LINE_WIDTH_PX
    return [(float(p[s]), float(p[e - 1])) for s, e in _runs(close) if e - s >= 3]


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

def read_csv(path):
    """A data CSV whose top comment lines (`#`) precede the header, as a structured array."""
    lines = [line for line in Path(path).read_text(encoding="utf-8").splitlines()
             if not line.lstrip().startswith("#") and line.strip()]
    return np.genfromtxt(lines, delimiter=",", names=True, comments="#")


def write_profile(trace: Trace, path, header_lines) -> Path:
    """`pressure_mbar, temperature_K, temperature_uncertainty_K`, top down, the last NaN except
    at the bars' pressures."""
    sigma = np.full(trace.pressure_mbar.shape, np.nan)
    for p_bar, half, _ in trace.bars:
        sigma[int(np.argmin(np.abs(np.log10(trace.pressure_mbar / p_bar))))] = half
    lines = [f"# {line}" for line in header_lines]
    lines.append("pressure_mbar,temperature_K,temperature_uncertainty_K")
    for p, T, s in zip(trace.pressure_mbar, trace.temperature_K, sigma):
        lines.append(f"{p:.5g},{T:.2f},{'' if not np.isfinite(s) else f'{s:.2f}'}")
    path = Path(path)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


# ---------------------------------------------------------------------------
# Kliore et al. (1980) Table 1, and the composition
# ---------------------------------------------------------------------------

#: Kliore et al.'s mixture and Lindal's, by number.
KLIORE_MIXTURE = {"H2": 0.90, "He": 0.10}
LINDAL_MIXTURE = {"H2": 0.94, "He": 0.06}
#: The refractivity set of `species_master.toml` the rescaling reads: Lindal's own.
REFRACTIVITY_SET = "lindal1985"


def mean_molar_mass(mixture, master="data_static/species_master.toml") -> float:
    with open(master, "rb") as handle:
        table = tomllib.load(handle)
    return float(sum(x * table["molar_mass"][s]["value"] for s, x in mixture.items()))


def mean_refractivity(mixture, master="data_static/species_master.toml") -> float:
    with open(master, "rb") as handle:
        table = tomllib.load(handle)
    sets = table["refractivity"][REFRACTIVITY_SET]
    return float(sum(x * sets[s]["refractivity_stp_e6"] for s, x in mixture.items()))


def rescale_composition(temperature_K, pressure_mbar, source=KLIORE_MIXTURE, target=LINDAL_MIXTURE,
                        master="data_static/species_master.toml"):
    """Lindal's rule (Fig. 5 caption): T in proportion to the mean molecular mass, p in proportion
    to the mean molecular mass over the mean molecular refractivity. Returns (T, p, factors)."""
    m = mean_molar_mass(target, master) / mean_molar_mass(source, master)
    r = mean_refractivity(target, master) / mean_refractivity(source, master)
    return (np.asarray(temperature_K, dtype="float64") * m, np.asarray(pressure_mbar, dtype="float64") * m / r,
            {"temperature": m, "pressure": m / r, "refractivity": r})


def read_pioneer11(path="occul_data/kliore1980/pioneer11_exit.csv"):
    return read_csv(path)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pdf", default="docs/Lindal_et_al_1985_AJ90_1136.pdf")
    args = parser.parse_args(argv)
    ink = render(args.pdf)
    _, traces = trace_all(ink)
    for name, (path, header) in OUTPUTS.items():
        trace = traces[name]
        write_profile(trace, path, header)
        print(f"wrote {path}: {trace.pressure_mbar.size} levels, {trace.pressure_mbar.min():.3g} to "
              f"{trace.pressure_mbar.max():.4g} mbar, {len(trace.bars)} error bars")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
