"""The Voyager 2 egress temperature profile, digitized, and its pressure. SPEC_12 Step 1.

Source: Lindal et al. (1985), AJ 90, 1136, Fig. 4 (page 4 of the PDF, journal page 1139): the
temperature against altitude above each profile's 1 bar level, for ingress (solid, labeled 36.3 N)
and egress (dashed, labeled 31.2 S), and the standard deviation sigma_T of the computed temperature
in a panel beside them.

**Digitizing**, every step coded so that running this reproduces the CSV:

1. Page 4 is rendered at 600 dpi by Ghostscript as a one-bit image (darker than mid-gray is ink).
2. The two panels' frames are found inside `FIGURE_BOX` from the rows and columns inked along most
   of it. Temperature is a least-squares line through its ticks, every 10 K from 90 to 150 K on the
   bottom edge, and sigma_T one through its tick at 10 K and the panel's two frame lines at 0 and
   20 K. Altitude, ticked every 50 km from 350 to -50 km on the left edge, is piecewise linear
   between its ticks and extended by the end intervals: the ticks depart from one straight line by
   up to 3.5 px (1.1 km), far more than a tick's position is uncertain, so the scan's scale varies
   along the axis and is taken from the ticks either side of each row.
3. Inside the temperature panel the two labels (`LABEL_BOXES`) and the two leader lines are
   blanked: the 36.3 N leader is the band `LEADER_36N`, the 31.2 S leader the connected piece
   holding `LEADER_31S_SEED`.
4. **The tracing.** Each profile is single-valued in altitude, but the curves are steep in some
   stretches and shallow in others, and where shallow they lie one above the other. So the panel is
   cut into `ZONES` by altitude, each scanned across the curves: by rows (temperature at each
   altitude) where the curves are steep, by columns (altitude at each temperature) where they are
   shallow. On each scan line the ink runs are found (runs `GAP_CLOSE_PX` or less apart joined),
   and the two curves are followed from the zone's `start`, kept in their order on each line: the
   zones end where the egress and ingress cross, near 175 km. Which of the two is the dashed curve
   is read from the figure for each zone (`ZONES`, `dashed`); as a check, the dashed curve should be
   found in more separate pieces along the scan, each dash starting one, and the record says whether
   it is. A single run is read as follows:
   - at least `SPLIT_WIDTH` line widths wide, with both curves predicted within it: the two lines
     side by side, each centre half a line in from its edge (a single line crossed obliquely is
     wide too, but has only one curve near it);
   - narrow, with the two predictions within a line width of each other: both curves;
   - narrow, one curve unseen for more than `HOLD_AFTER_LINES` lines: the other curve's;
   - otherwise: the curve on its side of the two predictions, the other absent.
   A track's prediction is linear from its last two positions, and held at its last once it has
   gone unseen for `HOLD_AFTER_LINES` lines. Each scan begins `ZONE_LEAD_KM` outside its zone.
5. Each curve's samples are reduced to the median of each `GRID_KM` bin of altitude; a bin with no
   sample, or a spike more than `SPIKE_K` from the median of its neighbors, is filled linearly in
   altitude from its neighbors. The sigma_T curve is traced by rows, kept within the altitudes the
   temperature curves span, and binned the same way.

**Altitude to pressure** (SPEC_12 section 1 deliverable 2, section 2 item 3). Fig. 4's altitude is a
distance along the local vertical above the profile's 1 bar level. From 1 bar,

    d ln p / dh = - m_bar g / (R T),

with `T` the profile's own, `m_bar` Lindal's composition at the level (`lindal_composition.nc`),
and `g` the magnitude of the effective gravity of `lib.gravity` (Lindal's gravity, rotation and
wind) at the profile's planetocentric latitude and at the radius `r = r_1bar + h cos(psi)`, the
reduction's own placement of an altitude (`lib.reduction.absolute_radius`). `r_1bar` is Lindal's
1 bar surface of `lindal_geodesy.nc` at that latitude, the ellipse of its equatorial and polar
radii. The integration is the trapezoid rule in `h` on the profile's own altitudes, repeated once
so that `m_bar` is read at the pressures it gives.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as ndi

PAGE = 4
DPI = 600
#: The figure on the 600 dpi page (left, top, right, bottom), generous around both panels.
FIGURE_BOX = (2563, 627, 4527, 2454)
#: The temperature panel's right limit: the axis break, left of the sigma_T panel.
TEMPERATURE_PANEL_RIGHT = 4040
ALTITUDE_TICKS_KM = tuple(float(v) for v in range(350, -51, -50))
TEMPERATURE_TICKS_K = tuple(float(v) for v in range(90, 151, 10))
SIGMA_TICK_K = 10.0
SIGMA_FRAME_K = (0.0, 20.0)
TICK_STRIP_PX = 25
#: Margins kept clear of the frame lines and their ticks, px: past the top and bottom ticks, and
#: narrow at the left, where the curves' minimum comes within about 30 px of the frame; the left
#: ticks themselves are blanked (`TICK_BLANK_PX` long, `TICK_HALF_PX` either side of each).
FRAME_MARGIN_PX = {"top": 45, "bottom": 45, "left": 10}
TICK_BLANK_PX = 45
TICK_HALF_PX = 6
#: The labels, blanked (left, top, right, bottom), page px at 600 dpi.
LABEL_BOXES = {"36.3 N": (2958, 1294, 3122, 1341), "31.2 S": (3028, 1589, 3182, 1634)}
#: The 36.3 N leader, a thin line joined to the solid curve: the band within `half_width` of the
#: segment from `start` to `stop` (x, y), page px, and its arrowhead, the box `arrowhead`, which
#: stands about 3 K left of the solid curve and touches it only at its tip.
LEADER_36N = {"start": (2909.0, 1326.0), "stop": (2933.0, 1460.0), "half_width": 7.0,
              "arrowhead": (2915, 1450, 2950, 1490)}
#: A pixel of the 31.2 S leader, which is not joined to either curve; its connected piece is blanked.
LEADER_31S_SEED = (1590, 2963)
#: Altitude zones (km, top first) and how each is scanned, and from which end the two curves are
#: followed (where both are present and apart).
#: `dashed` is the side of the two ordered tracks the egress curve lies on, read from the figure:
#: by rows, "first" is the cooler; by columns, "first" is the higher.
ZONES = (
    {"top_km": 400.0, "bottom_km": 172.0, "scan": "rows", "start": "bottom", "dashed": "first"},
    {"top_km": 172.0, "bottom_km": 140.0, "scan": "columns", "start": "low", "dashed": "second"},
    {"top_km": 140.0, "bottom_km": 70.0, "scan": "rows", "start": "top", "dashed": "second"},
    {"top_km": 70.0, "bottom_km": -30.0, "scan": "columns", "start": "low", "dashed": "second"},
)
#: A single run at least this many line widths wide holds both lines side by side; a line's width
#: is measured in each zone, along its scan, as this percentile of the zone's run widths.
SPLIT_WIDTH = 1.3
LINE_WIDTH_PERCENTILE = 30
#: Runs separated by at most this many px are one: a line broken by a hairline in the scan.
GAP_CLOSE_PX = 3
#: A binned value this far (K) from the median of its neighbours (this many either side) is a spike.
SPIKE_K = 2.0
SPIKE_NEIGHBOURS = 3
#: How far outside its zone each scan begins, km.
ZONE_LEAD_KM = 5.0
#: A track unseen for more than this many scan lines is predicted at its last position.
HOLD_AFTER_LINES = 30
GRID_KM = 1.0
COLUMNS = ("altitude_km", "pressure_mbar", "temperature_K", "temperature_uncertainty_K")
ONE_BAR_Pa = 1.0e5
#: Each profile's planetocentric latitude: Lindal's label converted as the anchor's is
#: (`tests/step11_2/target_latitude.py`); the ingress's is the anchor's own `phi_c`.
EGRESS_LATITUDE_PLANETOCENTRIC_DEG = -26.06095353762054
INGRESS_LATITUDE_PLANETOCENTRIC_DEG = 30.80556842739218


@dataclass(frozen=True)
class AltitudeScale:
    """Altitude (km) from page row, piecewise linear between the ticks, and its inverse."""

    rows: np.ndarray            #: tick rows, ascending
    km: np.ndarray              #: their altitudes, descending

    def __call__(self, row):
        row = np.asarray(row, dtype="float64")
        inside = np.interp(row, self.rows, self.km)
        low = self.km[0] + (row - self.rows[0]) * (self.km[1] - self.km[0]) / (self.rows[1] - self.rows[0])
        high = self.km[-1] + (row - self.rows[-1]) * (self.km[-1] - self.km[-2]) / (self.rows[-1] - self.rows[-2])
        return np.where(row < self.rows[0], low, np.where(row > self.rows[-1], high, inside))

    def row_of(self, km):
        km = np.asarray(km, dtype="float64")
        inside = np.interp(km, self.km[::-1], self.rows[::-1])
        low = self.rows[0] + (km - self.km[0]) * (self.rows[1] - self.rows[0]) / (self.km[1] - self.km[0])
        high = self.rows[-1] + (km - self.km[-1]) * (self.rows[-1] - self.rows[-2]) / (self.km[-1] - self.km[-2])
        return np.where(km > self.km[0], low, np.where(km < self.km[-1], high, inside))


@dataclass
class Trace:
    """One digitized curve: altitude (km) and temperature (K) on `GRID_KM` bins."""

    altitude_km: np.ndarray
    temperature_K: np.ndarray
    samples_h_km: np.ndarray
    samples_T_K: np.ndarray
    filled: np.ndarray          #: bins filled across a gap rather than measured


@dataclass
class Digitization:
    ink: np.ndarray
    frame: dict
    altitude_scale: AltitudeScale
    temperature_fit: tuple      #: (slope K/px, intercept K) in page columns
    sigma_fit: tuple
    altitude_residual_km: np.ndarray  #: the ticks' departure from one straight line, the reason
                                      #: the scale is piecewise
    temperature_residual_K: np.ndarray
    solid: Trace                #: ingress, 36.3 N
    dashed: Trace               #: egress, 31.2 S
    sigma: Trace                #: sigma_T, as a Trace of altitude and K
    zone_absence: list          #: per zone, the pieces each track was found in, the line width,
                                #: which track is dashed and whether the pieces agree


def find_ghostscript(given=None) -> str:
    for name in ([given] if given else ["gswin64c", "gswin32c", "gs"]):
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError("Ghostscript is needed to render the PDF and none was found; pass --gs.")


def render(pdf, gs=None) -> np.ndarray:
    with tempfile.TemporaryDirectory() as directory:
        out = Path(directory) / "page.png"
        subprocess.run([find_ghostscript(gs), "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                        "-sDEVICE=pnggray", f"-r{DPI}", f"-dFirstPage={PAGE}", f"-dLastPage={PAGE}",
                        f"-sOutputFile={out}", str(pdf)], check=True)
        with Image.open(out) as image:
            return np.asarray(image.convert("L")) < 128


def _runs(mask_1d, offset=0):
    """(start, stop) of each run of True, stop exclusive, shifted by `offset`."""
    edges = np.flatnonzero(np.diff(np.concatenate([[0], mask_1d.astype(np.int8), [0]])))
    return [(int(a) + offset, int(b) + offset) for a, b in zip(edges[::2], edges[1::2])]


def _line_runs(mask_1d):
    """(centre, width) of each run on a scan line, runs `GAP_CLOSE_PX` or less apart joined."""
    joined = []
    for a, b in _runs(mask_1d):
        if joined and a - joined[-1][1] <= GAP_CLOSE_PX:
            joined[-1] = (joined[-1][0], b)
        else:
            joined.append((a, b))
    return [(0.5 * (a + b - 1), b - a) for a, b in joined]


def find_frame(ink) -> dict:
    """The frame lines: top and bottom rows, the temperature panel's left edge, the sigma panel's
    left and right edges (each the line's center), from the lines inked along most of the box."""
    x0, y0, x1, y1 = FIGURE_BOX
    box = ink[y0:y1, x0:x1]
    rows = _runs(box.sum(axis=1) > 0.5 * box.shape[1], y0)
    cols = _runs(box.sum(axis=0) > 0.4 * box.shape[0], x0)
    centre = lambda run: 0.5 * (run[0] + run[1] - 1)
    if len(rows) != 2 or len(cols) != 3:
        raise ValueError(f"found {len(rows)} frame rows and {len(cols)} frame columns; the figure "
                         "has 2 and 3")
    return {"top": rows[0], "bottom": rows[1], "left": cols[0], "sigma_left": cols[1],
            "sigma_right": cols[2], "centres": {"sigma_left": centre(cols[1]),
                                               "sigma_right": centre(cols[2])}}


def _ticks(strip, axis, offset, low, high):
    centres = [0.5 * (a + b - 1) for a, b in _runs(strip.mean(axis=axis) > 0.5, offset)]
    return np.array([c for c in centres if low + 8 < c < high - 8])


def calibrate(ink, frame):
    top, bottom = frame["top"][1], frame["bottom"][0]
    left = frame["left"][1]
    rows = _ticks(ink[top:bottom, left:left + TICK_STRIP_PX], 1, top, top, bottom)
    cols = _ticks(ink[bottom - TICK_STRIP_PX:bottom, left:TEMPERATURE_PANEL_RIGHT], 0, left, left,
                  TEMPERATURE_PANEL_RIGHT)
    cols = cols[cols < TEMPERATURE_PANEL_RIGHT - 60]
    if rows.size != len(ALTITUDE_TICKS_KM) or cols.size != len(TEMPERATURE_TICKS_K):
        raise ValueError(f"found {rows.size} altitude and {cols.size} temperature ticks; the figure "
                         f"has {len(ALTITUDE_TICKS_KM)} and {len(TEMPERATURE_TICKS_K)}")
    frame["tick_rows"] = rows
    a_fit = np.polyfit(rows, ALTITUDE_TICKS_KM, 1)
    scale = AltitudeScale(rows=np.asarray(rows, dtype="float64"),
                          km=np.asarray(ALTITUDE_TICKS_KM, dtype="float64"))
    t_fit = np.polyfit(cols, TEMPERATURE_TICKS_K, 1)
    s_left, s_right = frame["sigma_left"][1], frame["sigma_right"][0]
    sigma_cols = _ticks(ink[bottom - TICK_STRIP_PX:bottom, s_left:s_right], 0, s_left, s_left, s_right)
    if sigma_cols.size != 1:
        raise ValueError(f"found {sigma_cols.size} sigma_T ticks inside the panel; the figure has 1")
    points = np.array([frame["centres"]["sigma_left"], sigma_cols[0], frame["centres"]["sigma_right"]])
    s_fit = np.polyfit(points, [SIGMA_FRAME_K[0], SIGMA_TICK_K, SIGMA_FRAME_K[1]], 1)
    return (scale, tuple(t_fit), tuple(s_fit),
            np.polyval(a_fit, rows) - np.array(ALTITUDE_TICKS_KM),
            np.polyval(t_fit, cols) - np.array(TEMPERATURE_TICKS_K))


def curve_mask(ink, frame, tick_rows) -> np.ndarray:
    """The temperature panel's ink, clear of the frame, its ticks, the labels and the leaders."""
    top, bottom, left = frame["top"][1], frame["bottom"][0], frame["left"][1]
    m = FRAME_MARGIN_PX
    mask = np.zeros_like(ink)
    window = (slice(top + m["top"], bottom - m["bottom"]),
              slice(left + m["left"], TEMPERATURE_PANEL_RIGHT))
    mask[window] = ink[window]
    for row in tick_rows:
        mask[int(row) - TICK_HALF_PX:int(row) + TICK_HALF_PX + 1, left:left + TICK_BLANK_PX] = False
    for x0, y0, x1, y1 in LABEL_BOXES.values():
        mask[y0:y1, x0:x1] = False
    x0, y0, x1, y1 = LEADER_36N["arrowhead"]
    mask[y0:y1, x0:x1] = False
    (xa, ya), (xb, yb) = LEADER_36N["start"], LEADER_36N["stop"]
    yy, xx = np.mgrid[int(ya) - 10:int(yb) + 10, int(min(xa, xb)) - 15:int(max(xa, xb)) + 15]
    t = np.clip(((xx - xa) * (xb - xa) + (yy - ya) * (yb - ya)) / ((xb - xa) ** 2 + (yb - ya) ** 2),
                0.0, 1.0)
    near = np.hypot(xx - (xa + t * (xb - xa)), yy - (ya + t * (yb - ya))) <= LEADER_36N["half_width"]
    mask[yy[near], xx[near]] = False
    labels, _ = ndi.label(mask, structure=np.ones((3, 3)))
    seed = labels[LEADER_31S_SEED]
    if seed == 0:
        raise ValueError("the 31.2 S leader's seed pixel is not ink")
    mask[labels == seed] = False
    return mask


def _follow(lines, start_reversed, single_width):
    """Two tracks through a sequence of scan lines, each a list of (centre, width) runs.

    Returns two arrays of positions (NaN where absent) in the input order."""
    order = range(len(lines) - 1, -1, -1) if start_reversed else range(len(lines))
    first = np.full(len(lines), np.nan)
    second = np.full(len(lines), np.nan)
    history = ([], [])

    def predict(track, at):
        # Linear from the last two positions, held at the last once the track has gone unseen for
        # more than `HOLD_AFTER_LINES` scan lines, so a curve that has ended is not extrapolated.
        h = history[track]
        if not h:
            return None
        (i1, v1) = h[-1]
        if len(h) == 1 or abs(at - i1) > HOLD_AFTER_LINES:
            return v1
        (i0, v0) = h[-2]
        return v1 + (v1 - v0) / (i1 - i0) * (at - i1)

    def stale(track, at):
        h = history[track]
        return bool(h) and abs(at - h[-1][0]) > HOLD_AFTER_LINES

    for i in order:
        runs = lines[i]
        if not runs:
            continue
        p = (predict(0, i), predict(1, i))
        if p[0] is None:
            if len(runs) >= 2:
                runs = sorted(runs)[:2]
                first[i], second[i] = runs[0][0], runs[1][0]
            elif runs[0][1] >= SPLIT_WIDTH * single_width:
                centre, width = runs[0]
                half = 0.5 * (width - single_width)
                first[i], second[i] = centre - half, centre + half
            else:
                continue
        elif len(runs) == 1:
            centre, width = runs[0]
            reach = 0.5 * width + single_width
            if (width >= SPLIT_WIDTH * single_width and abs(p[0] - centre) <= reach
                    and abs(p[1] - centre) <= reach):
                # Two lines side by side in one run: each centre half a line in from its edge.
                half = 0.5 * (width - single_width)
                first[i], second[i] = centre - half, centre + half
            elif abs(p[1] - p[0]) < single_width:
                # The two lines on top of each other: both curves are here.
                first[i] = second[i] = centre
            elif stale(0, i) != stale(1, i):
                # One curve has ended (unseen for more than `HOLD_AFTER_LINES` lines): a lone run
                # is the other's.
                if stale(0, i):
                    second[i] = centre
                else:
                    first[i] = centre
            elif centre <= 0.5 * (p[0] + p[1]):
                first[i] = centre
            else:
                second[i] = centre
        else:
            best = None
            for a in range(len(runs)):
                for b in range(len(runs)):
                    if a == b:
                        continue
                    cost = abs(runs[a][0] - p[0]) + abs(runs[b][0] - p[1])
                    if best is None or cost < best[0]:
                        best = (cost, a, b)
            # The curves keep their order within a zone (`ZONES` ends at their crossings).
            pair = sorted((runs[best[1]][0], runs[best[2]][0]))
            first[i], second[i] = pair
        for track, value in ((0, first[i]), (1, second[i])):
            if np.isfinite(value):
                history[track].append((i, value))
    return first, second


def _line_width(lines) -> float:
    """One line's width across this zone's scan direction: the 30th percentile of its run widths,
    where lone lines, not two lines touching, are most of the runs (`LINE_WIDTH_PERCENTILE`)."""
    widths = [w for runs in lines for _, w in runs if w >= 4]
    return float(np.percentile(widths, LINE_WIDTH_PERCENTILE))


def trace_curves(mask, frame, altitude_scale, temperature_fit):
    """The solid and dashed samples (h km, T K) and each zone's absent fractions."""
    to_h = altitude_scale
    to_T = lambda col: np.polyval(temperature_fit, col)
    row_of = altitude_scale.row_of
    samples = {"solid": ([], []), "dashed": ([], [])}
    absence = []
    left = frame["left"][1] + FRAME_MARGIN_PX["left"]
    for zone in ZONES:
        # The scan begins `ZONE_LEAD_KM` outside the zone on its starting side (both sides for a
        # column zone), so the two tracks are settled when they enter it; only samples inside count.
        lead_top = ZONE_LEAD_KM if zone["start"] in ("top", "low", "high") else 0.0
        lead_bottom = ZONE_LEAD_KM if zone["start"] in ("bottom", "low", "high") else 0.0
        r0 = int(np.ceil(max(row_of(zone["top_km"] + lead_top), frame["top"][1] + FRAME_MARGIN_PX["top"])))
        r1 = int(np.floor(min(row_of(zone["bottom_km"] - lead_bottom),
                              frame["bottom"][0] - FRAME_MARGIN_PX["bottom"])))
        part = mask[r0:r1]
        if zone["scan"] == "rows":
            lines = [_line_runs(part[i]) for i in range(part.shape[0])]
            single = _line_width(lines)
            a, b = _follow(lines, zone["start"] == "bottom", single)
            coord = np.arange(r0, r1)
            pairs = [(to_h(coord), to_T(a)), (to_h(coord), to_T(b))]
        else:
            sub = part[:, left:TEMPERATURE_PANEL_RIGHT]
            lines = [_line_runs(sub[:, j]) for j in range(sub.shape[1])]
            single = _line_width(lines)
            a, b = _follow(lines, zone["start"] == "high", single)
            coord = np.arange(left, TEMPERATURE_PANEL_RIGHT)
            pairs = [(to_h(r0 + a), to_T(coord)), (to_h(r0 + b), to_T(coord))]
        segments = []
        for track in (a, b):
            present = np.isfinite(track)
            segments.append(int(present[0]) + int(np.sum(present[1:] & ~present[:-1])))
        dashed_index = 1 if zone["dashed"] == "second" else 0
        absence.append({"zone": zone, "segments": segments, "line_width_px": single,
                        "dashed": zone["dashed"],
                        "pieces_agree": segments[dashed_index] >= segments[1 - dashed_index]})
        for k, (h, T) in enumerate(pairs):
            keep = np.isfinite(h) & np.isfinite(T)
            keep &= (h <= zone["top_km"]) & (h >= zone["bottom_km"])
            name = "dashed" if k == dashed_index else "solid"
            samples[name][0].append(h[keep])
            samples[name][1].append(T[keep])
    return ({name: (np.concatenate(hs), np.concatenate(Ts)) for name, (hs, Ts) in samples.items()},
            absence, single)


def binned(h, T) -> Trace:
    """`GRID_KM` bins from the lowest to the highest sample; empty bins filled linearly."""
    order = np.argsort(h)
    h, T = h[order], T[order]
    grid = np.arange(np.ceil(h.min() / GRID_KM) * GRID_KM, np.floor(h.max() / GRID_KM) * GRID_KM + 1e-9,
                     GRID_KM)
    index = np.round(h / GRID_KM).astype(int) - int(round(grid[0] / GRID_KM))
    ok = (index >= 0) & (index < grid.size)
    counts = np.bincount(index[ok], minlength=grid.size)
    value = np.full(grid.size, np.nan)
    for i in np.flatnonzero(counts):
        value[i] = np.median(T[ok][index[ok] == i])
    filled = counts == 0
    value[filled] = np.interp(grid[filled], grid[~filled], value[~filled])
    # A bin more than `SPIKE_K` from the median of its `SPIKE_NEIGHBOURS` neighbours either side is
    # a sample of the other line taken across a gap; it is refilled from its neighbours.
    n = SPIKE_NEIGHBOURS
    spike = np.zeros(grid.size, dtype=bool)
    for i in range(grid.size):
        around = np.concatenate([value[max(0, i - n):i], value[i + 1:i + 1 + n]])
        spike[i] = around.size >= n and abs(value[i] - np.median(around)) > SPIKE_K
    if spike.any():
        value[spike] = np.interp(grid[spike], grid[~spike], value[~spike])
    filled |= spike
    return Trace(altitude_km=grid, temperature_K=value, samples_h_km=h, samples_T_K=T, filled=filled)


def trace_sigma(ink, frame, altitude_scale, sigma_fit):
    """sigma_T against altitude, by rows inside the sigma panel."""
    top, bottom = frame["top"][1], frame["bottom"][0]
    left, right = frame["sigma_left"][1] + 10, frame["sigma_right"][0] - 10
    rows, values = [], []
    for row in range(top + FRAME_MARGIN_PX["top"], bottom - FRAME_MARGIN_PX["bottom"]):
        runs = _runs(ink[row, left:right], left)
        if len(runs) == 1:
            rows.append(row)
            values.append(0.5 * (runs[0][0] + runs[0][1] - 1))
    rows, values = np.array(rows), np.array(values)
    return altitude_scale(rows), np.polyval(sigma_fit, values)


def digitize(pdf, gs=None) -> Digitization:
    ink = render(pdf, gs)
    frame = find_frame(ink)
    a_fit, t_fit, s_fit, a_res, t_res = calibrate(ink, frame)
    mask = curve_mask(ink, frame, frame["tick_rows"])
    samples, absence, _ = trace_curves(mask, frame, a_fit, t_fit)
    solid, dashed = binned(*samples["solid"]), binned(*samples["dashed"])
    s_h, s_v = trace_sigma(ink, frame, a_fit, s_fit)
    low = min(solid.altitude_km.min(), dashed.altitude_km.min())
    high = max(solid.altitude_km.max(), dashed.altitude_km.max())
    keep = (s_h >= low) & (s_h <= high)
    sigma = binned(s_h[keep], s_v[keep])
    return Digitization(ink=ink, frame=frame, altitude_scale=a_fit, temperature_fit=t_fit,
                        sigma_fit=s_fit, altitude_residual_km=a_res, temperature_residual_K=t_res,
                        solid=solid, dashed=dashed, sigma=sigma, zone_absence=absence)


# ---------------------------------------------------------------------------
# Altitude to pressure
# ---------------------------------------------------------------------------

@dataclass
class Profile:
    """A profile on its altitude grid with its pressure."""

    altitude_km: np.ndarray
    pressure_Pa: np.ndarray
    temperature_K: np.ndarray
    temperature_uncertainty_K: np.ndarray
    latitude_planetocentric_deg: float
    radius_1bar_m: float
    mean_molar_mass_kg_mol: np.ndarray
    gravity_ms2: np.ndarray


def one_bar_radius(geodesy, latitude_planetocentric_deg) -> float:
    """Lindal's 1 bar surface at a planetocentric latitude: the ellipse of its two radii."""
    row = int(np.flatnonzero(np.asarray(geodesy["surface_pressure_Pa"].values) == ONE_BAR_Pa)[0])
    a = float(geodesy["radius_equatorial_m"].values[row])
    b = float(geodesy["radius_polar_m"].values[row])
    phi = np.radians(latitude_planetocentric_deg)
    return a * b / np.sqrt((b * np.cos(phi)) ** 2 + (a * np.sin(phi)) ** 2)


def to_pressure(altitude_km, temperature_K, latitude_planetocentric_deg, composition, gravity,
                rotation, wind, geodesy, sigma_K=None) -> Profile:
    """Hydrostatic pressure from 1 bar at altitude zero, as the module docstring states."""
    from casspian.lib.constants import MOLAR_GAS_CONSTANT
    from casspian.lib.gravity import g_eff_vector
    from casspian.refrac.anchor import wind_of_latitude
    from casspian.tools.lindal.lindal_wind import mean_molar_mass

    h = np.asarray(altitude_km, dtype="float64") * 1e3
    T = np.asarray(temperature_K, dtype="float64")
    phi = np.radians(float(latitude_planetocentric_deg))
    constants = (float(rotation["angular_rate_rad_s"]), float(gravity["GM_m3s2"]),
                 np.asarray(gravity["J"].values, dtype="float64"),
                 np.asarray(gravity["degree"].values), float(gravity["normalization_radius_m"]))
    u = float(wind_of_latitude(wind)(np.array([phi]))[0])
    r1 = one_bar_radius(geodesy, latitude_planetocentric_deg)
    psi = float(g_eff_vector(u, r1, phi, *constants)[3])
    radius = r1 + h * np.cos(psi)
    g = np.asarray(g_eff_vector(np.full(h.shape, u), radius, np.full(h.shape, phi), *constants)[2])
    zero = int(np.argmin(np.abs(h)))
    if h[zero] != 0.0:
        raise ValueError("the profile's altitudes do not include zero, the 1 bar level")
    p = np.full(h.shape, ONE_BAR_Pa)
    for _ in range(2):
        m = mean_molar_mass(composition, p, np.array([float(latitude_planetocentric_deg)]))[:, 0]
        rate = -m * g / (MOLAR_GAS_CONSTANT * T)
        step = 0.5 * (rate[1:] + rate[:-1]) * np.diff(h)
        ln_p = np.concatenate([[0.0], np.cumsum(step)])
        p = ONE_BAR_Pa * np.exp(ln_p - ln_p[zero])
    sigma = np.full(h.shape, np.nan) if sigma_K is None else np.asarray(sigma_K, dtype="float64")
    return Profile(altitude_km=h / 1e3, pressure_Pa=p, temperature_K=T,
                   temperature_uncertainty_K=sigma,
                   latitude_planetocentric_deg=float(latitude_planetocentric_deg),
                   radius_1bar_m=r1, mean_molar_mass_kg_mol=m, gravity_ms2=g)


def sigma_on(trace: Trace, sigma: Trace) -> np.ndarray:
    """sigma_T at a trace's altitudes, NaN beyond the sigma_T curve's span."""
    out = np.interp(trace.altitude_km, sigma.altitude_km, sigma.temperature_K)
    beyond = (trace.altitude_km < sigma.altitude_km.min()) | (trace.altitude_km > sigma.altitude_km.max())
    return np.where(beyond, np.nan, out)


def write_csv(profile: Profile, path) -> Path:
    path = Path(path)
    lines = [",".join(COLUMNS)]
    for h, p, T, s in zip(profile.altitude_km, profile.pressure_Pa, profile.temperature_K,
                          profile.temperature_uncertainty_K):
        lines.append(f"{h:.1f},{p / 100.0:.6g},{T:.2f},{'' if not np.isfinite(s) else f'{s:.2f}'}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def _load(path, kind):
    from casspian.lib import io as cio

    handle = cio.read(path, kind)
    try:
        return handle.load()
    finally:
        handle.close()


def lindal_inputs(directory="occul_data/lindal"):
    """Lindal's composition, gravity, rotation, wind and geodesy, the reduction's own."""
    d = Path(directory)
    return {"composition": _load(d / "lindal_composition.nc", "composition"),
            "gravity": _load(d / "lindal_gravity.nc", "gravity"),
            "rotation": _load(d / "lindal_rotation.nc", "rotation"),
            "wind": _load(d / "lindal_wind.nc", "wind"),
            "geodesy": _load(d / "lindal_geodesy.nc", "geodesy")}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pdf", default="docs/Lindal_et_al_1985_AJ90_1136.pdf")
    parser.add_argument("--output", default="occul_data/lindal/voyager2_egress.csv")
    parser.add_argument("--gs", default=None, help="the Ghostscript executable")
    args = parser.parse_args(argv)
    result = digitize(args.pdf, args.gs)
    inputs = lindal_inputs()
    egress = to_pressure(result.dashed.altitude_km, result.dashed.temperature_K,
                         EGRESS_LATITUDE_PLANETOCENTRIC_DEG, sigma_K=sigma_on(result.dashed, result.sigma),
                         **inputs)
    path = write_csv(egress, args.output)
    print(f"wrote {path}: {egress.altitude_km.size} levels, {egress.altitude_km.min():g} to "
          f"{egress.altitude_km.max():g} km, {egress.pressure_Pa.max() / 100:.4g} to "
          f"{egress.pressure_Pa.min() / 100:.4g} mbar")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
