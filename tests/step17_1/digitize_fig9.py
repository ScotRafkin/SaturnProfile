"""Lindal et al. (1985) Fig. 9, right panel, digitized: SPEC_17 v0.2 Step 1 item 1.

The panel gives the altitude of the 100 mbar surface above Lindal's reference geoid against
planetocentric latitude: five points, and the stippled curve of the geoid that best fits them.

1. Journal page 1142 (page 7 of the PDF) is rendered at 600 dpi by Ghostscript as a one-bit image
   (darker than mid-gray is ink), as SPEC_15 rendered Figs. 5 and 6.
2. The panel's frame is the top and bottom rows and the right column inked along most of the box;
   the zero-altitude line is the long column at its left. Altitude: a least-squares line in column
   through the top axis's ticks every 10 km, from the zero line (0 km) to the right frame (140 km).
   Latitude: a least-squares line in row through the right axis's ticks every 10 deg, from the top
   frame (90 N) to the bottom (90 S). The ticks are drafted unevenly, by up to 5 px (0.5 deg).
3. The points are the filled dots: what survives an opening of the ink by a disc of 11 px radius
   (the dashed curve is 8 to 10 px wide, a dot 30 px across), at the centroid of each.
4. The stippled curve is scanned row by row between the frames, the dots, the ticks and the zero line
   masked; in each row the run of ink nearest the row before is taken, and dash gaps are left empty.

Writes `occul_data/lindal/isobar_100mbar.csv` (the points) when asked (`--write`), and returns
the points, the curve and the calibration to the acceptance script. Run from the repository root.
"""

import argparse
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from PIL import Image
from scipy import ndimage as nd

from casspian.tools.lindal.egress_profile import find_ghostscript

PDF = "docs/Lindal_et_al_1985_AJ90_1136.pdf"
PAGE = 7
DPI = 600
#: The right panel's box on the page, in pixels at 600 dpi (column, row, column, row).
BOX = (2650, 500, 4350, 2500)
DOT_RADIUS_PX = 11
OUT = Path("occul_data/lindal/isobar_100mbar.csv")
#: The points from north to south, as Lindal labels them in the left panel, and the occultation
#: each is (Lindal et al. 1985, p. 1141): Pioneer 11's is from radio tracking, not an occultation.
OBSERVATIONS = ("voyager2", "voyager1", "pioneer11", "voyager2", "voyager1")
HEADER = """\
# Altitude of the 100 mbar isobaric surface above Lindal's reference geoid, against planetocentric
# latitude. Lindal, G. F., Sweetnam, D. N., and Eshleman, V. R. 1985, Astron. J. 90, 1136, Fig. 9,
# right panel (journal page 1142), digitized by tests/step17_1/digitize_fig9.py (SPEC_17 Step 1
# item 1): rendered at 600 dpi, calibrated by the ticks, each point the centroid of its dot. The
# ticks are drafted unevenly, by up to 0.5 deg in latitude; the calibration RMS is in the report.
# The reference geoid (Fig. 9 caption and Appendix Eqs. 12 to 14): the level surface through a polar
# radius of 54,438 km of the external gravity field and a uniform rotation at the System III period,
# 10 h 39 m 22.4 s, with no wind. Lindal gives no error bars on the points; the geoids he fits to
# the isobaric surfaces fit them to 3 to 4 km (p. 1141).
# The points and their epochs:
#   voyager2   Voyager 2 ingress, 36.3 N planetographic, the anchor, 26 August 1981
#   voyager1   Voyager 1 egress, 3.0 S planetographic, 13 November 1980
#   pioneer11  Pioneer 11, from radio tracking, not occultation (p. 1141), September 1979; for
#              reference only, two years from the anchor (SPEC_17 section 0)
#   voyager2   Voyager 2 egress, 31.2 S planetographic, 26 August 1981
#   voyager1   Voyager 1 ingress, 74.2 S planetographic, 13 November 1980
# Latitudes are planetocentric, as Lindal plots them.
"""


@dataclass
class Digitized:
    latitude_deg: np.ndarray        #: the points, north to south
    altitude_km: np.ndarray
    observation: tuple
    curve_latitude_deg: np.ndarray  #: the stippled curve, one value per row with ink
    curve_altitude_km: np.ndarray
    altitude_rms_km: float          #: the ticks about their fits
    latitude_rms_deg: float
    altitude_ticks: int
    latitude_ticks: int


def render(pdf=PDF) -> np.ndarray:
    """The page at `DPI`, ink True, cropped to `BOX`."""
    with tempfile.TemporaryDirectory() as directory:
        out = Path(directory) / "page.png"
        subprocess.run([find_ghostscript(None), "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                        "-sDEVICE=pnggray", f"-r{DPI}", f"-dFirstPage={PAGE}", f"-dLastPage={PAGE}",
                        f"-sOutputFile={out}", str(pdf)], check=True)
        Image.MAX_IMAGE_PIXELS = None
        page = np.asarray(Image.open(out)) < 128
    x0, y0, x1, y1 = BOX
    return page[y0:y1, x0:x1]


def _centres(mask, offset=0):
    labels, n = nd.label(mask)
    return np.array([offset + np.flatnonzero(labels == k).mean() for k in range(1, n + 1)])


def _fit(positions, values):
    a, b = np.polyfit(positions, values, 1)
    return a, b, float(np.sqrt(np.mean((a * positions + b - values) ** 2)))


def digitize(ink) -> Digitized:
    height, width = ink.shape
    rows = _centres(ink.sum(axis=1) > 0.4 * width)
    cols = _centres(ink.sum(axis=0) > 0.45 * height)
    if rows.size != 2 or cols.size != 3:
        raise ValueError(f"found {rows.size} frame rows and {cols.size} long columns; the panel has "
                         "2 rows and 3 columns (the left panel's frame, the zero line, the right frame)")
    top, bottom = rows
    zero, right = cols[1], cols[2]

    # Altitude: the ticks hang below the top frame; each is assigned its 10 km by the frame.
    below = ink[int(top) + 4:int(top) + 14, :].sum(axis=0) >= 7
    x = _centres(below)
    x = x[(x > zero - 5) & (x < right + 5)]
    rough = (x - zero) / (right - zero) * 140.0
    km = np.round(rough / 10.0) * 10.0
    keep = np.abs(rough - km) < 3.0
    a_x, b_x, rms_x = _fit(x[keep], km[keep])

    # Latitude: the ticks stand left of the right frame, the frame rows at 90 N and 90 S.
    left = ink[:, int(right) - 18:int(right) - 4].sum(axis=1) >= 13
    y = _centres(left)
    y = y[(y > top - 5) & (y < bottom + 5)]
    rough = 90.0 - (y - top) / (bottom - top) * 180.0
    deg = np.round(rough / 10.0) * 10.0
    keep_y = np.abs(rough - deg) < 3.0
    a_y, b_y, rms_y = _fit(y[keep_y], deg[keep_y])

    # The dots.
    disc = np.hypot(*np.mgrid[-DOT_RADIUS_PX:DOT_RADIUS_PX + 1, -DOT_RADIUS_PX:DOT_RADIUS_PX + 1]) <= DOT_RADIUS_PX
    opened = nd.binary_opening(ink, structure=disc)
    opened[:, int(right) - 2:] = False
    labels, n = nd.label(opened)
    dots = []
    for k in range(1, n + 1):
        ys, xs = np.nonzero(labels == k)
        if ys.size > 200:
            dots.append((ys.mean(), xs.mean()))
    if len(dots) != len(OBSERVATIONS):
        raise ValueError(f"found {len(dots)} dots; Fig. 9's right panel has {len(OBSERVATIONS)}")
    dots.sort()
    dot_y, dot_x = np.array(dots).T
    latitude = a_y * dot_y + b_y
    altitude = a_x * dot_x + b_x

    # The stippled curve.
    masked = ink.copy()
    masked[:int(top) + 30, :] = False
    masked[int(bottom) - 30:, :] = False
    masked[:, int(right) - 25:] = False
    masked[:, int(zero) - 5:int(zero) + 6] = False
    masked[:, :int(zero) - 120] = False
    masked &= ~nd.binary_dilation(opened, iterations=6)
    curve_y, curve_x = [], []
    previous = None
    for row in range(int(top) + 30, int(bottom) - 30):
        runs = _runs(masked[row])
        if not runs:
            continue
        centres = np.array([0.5 * (s + e - 1) for s, e in runs])
        pick = centres[np.argmin(np.abs(centres - previous))] if previous is not None else centres[0]
        if previous is not None and abs(pick - previous) > 40:
            continue
        previous = pick
        curve_y.append(row)
        curve_x.append(pick)
    curve_y, curve_x = np.array(curve_y, dtype="float64"), np.array(curve_x)
    return Digitized(latitude, altitude, OBSERVATIONS, a_y * curve_y + b_y, a_x * curve_x + b_x,
                     rms_x, rms_y, int(keep.sum()), int(keep_y.sum()))


def _runs(line):
    edges = np.flatnonzero(np.diff(np.concatenate([[0], line.astype(np.int8), [0]])))
    return list(zip(edges[::2], edges[1::2]))


def write(result: Digitized, path=OUT):
    lines = [HEADER, "latitude_planetocentric_deg,altitude_km,observation\n"]
    lines += [f"{lat:.2f},{alt:.1f},{obs}\n" for lat, alt, obs in
              zip(result.latitude_deg, result.altitude_km, result.observation)]
    Path(path).write_bytes("".join(lines).encode("utf-8"))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--write", action="store_true", help=f"write {OUT}")
    args = parser.parse_args(argv)
    result = digitize(render())
    print(f"altitude: {result.altitude_ticks} ticks, RMS {result.altitude_rms_km:.3f} km; "
          f"latitude: {result.latitude_ticks} ticks, RMS {result.latitude_rms_deg:.3f} deg")
    for lat, alt, obs in zip(result.latitude_deg, result.altitude_km, result.observation):
        print(f"  {obs:10s} {lat:+7.2f} deg  {alt:6.1f} km")
    print(f"curve: {result.curve_latitude_deg.size} rows, {result.curve_latitude_deg.min():.1f} to "
          f"{result.curve_latitude_deg.max():.1f} deg")
    if args.write:
        write(result)
        print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
