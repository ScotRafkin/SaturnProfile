"""Digitize the Voyager IRIS temperatures at three levels. SPEC_09 Step 1.

Source: Conrath and Pirraglia (1983), Icarus 53, 286, Fig. 1 (page 2 of the PDF, journal page
287): retrieved temperatures against latitude at three pressure levels, one dot per retrieval.
The same data are Lindal et al. (1985) Fig. 10. The top level is written as 110 mbar, the label
of Lindal and of Conrath et al. (1984); Conrath and Pirraglia label it 150 mbar (SPEC_09 section
2, ruling 1). Latitude is planetographic (ruling 2).

Method, every step coded so that running this reproduces the CSV:

1. The page is rendered at 600 dpi by Ghostscript, as a one-bit image (ink is any pixel darker
   than mid-gray).
2. The plot frame is found inside `FIGURE_BOX` from the rows and columns that are inked along
   most of the box.
3. Both axes are calibrated by their ticks. The ticks are the runs of rows (left edge) or columns
   (bottom edge) inked across most of a strip just inside the frame. The temperature ticks are
   every 5 K from 120 down to 80 K, the latitude ticks every 10 deg from 80 to -80 (the frame's own
   edges are 90 and -90 and carry no separate tick). Each axis is a least-squares line through
   its ticks; the residuals are reported.
4. Inside the frame the three level labels (`LABEL_BOXES`) are blanked, and the image is opened by
   a disk of radius `OPEN_RADIUS`. A dot is about 11 px across and a stroke of the broken curve,
   a tick or the frame about 4 to 7 px, so the opening keeps the dots and removes the rest.
5. Each connected component holds `k = round(area / A1)` dots, `A1` the median area of the
   components between half and one and a half times the median of all. A component with
   `k > 1` is split into `k` centers by k-means on its pixels (fixed seed). This is the method of
   `data_static/winds/digitize_smith1982_fig4.py`. A few dots are printed narrower than the
   opening's disk and are removed with the strokes; the ink the opening removed is taken back as
   a dot where it has the size, shape and thickness of one (`NARROW_DOT`).
6. A center in `CURVE_REGION` (south of -18 deg and colder than 85 K, where the broken curve runs
   and no dot lies) is a remnant of the curve and is dropped.
7. Each dot is assigned to its level by temperature: 730 mbar above `BOUNDARY_730_290_K`, 150 mbar
   (written 110) below the line `BOUNDARY_290_150`, 290 mbar between. The second boundary is a
   line in latitude because the 290 and 150 mbar dots come within about 4 K of each other near
   75 N and near -90.

The figure's printed values are not interpreted: no smoothing, no merging of dots, nothing
dropped except by the rules above.

`fit` (SPEC_10 Step 1) is the local polynomial fit of the digitized temperatures in latitude, with
its gradient and the gradient's standard error; `running_mean` (SPEC_09) only characterizes the
spread.
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
from scipy.cluster.vq import kmeans2

PAGE = 2
DPI = 600
#: The figure on the 600 dpi page (left, top, right, bottom), generous around the frame.
FIGURE_BOX = (333, 500, 3708, 2791)
#: The level labels on the 600 dpi page, measured from the scan (left, top, right, bottom).
LABEL_BOXES = {"730 mbar": (1488, 632, 1722, 678), "290 mbar": (1412, 1618, 1648, 1664),
               "150 mbar": (1284, 2214, 1514, 2262)}
TEMPERATURE_TICKS_K = (120.0, 115.0, 110.0, 105.0, 100.0, 95.0, 90.0, 85.0, 80.0)
LATITUDE_TICKS_DEG = tuple(float(v) for v in range(80, -81, -10))
TICK_STRIP_PX = 40
OPEN_RADIUS = 4.0
KMEANS_SEED = 1
#: Step 5b: a mark the opening removed is a dot if its area is within these fractions of a single
#: dot's, its bounding box no more elongated than `aspect`, and its ink reaches `thickness_px`
#: from its edge (measured: dots 4.0, strokes of the broken curve 3.6 at most outside its ends).
NARROW_DOT = {"area_low": 0.6, "area_high": 1.4, "aspect": 1.75, "thickness_px": 4.0}
CURVE_REGION = {"south_of_deg": -18.0, "colder_than_K": 85.0}
BOUNDARY_730_290_K = 105.0
#: (planetographic latitude deg, temperature K) nodes of the 290 / 150 mbar boundary.
BOUNDARY_290_150 = ((90.0, 86.0), (40.0, 88.0), (15.0, 92.0), (0.0, 90.0), (-15.0, 93.5),
                    (-90.0, 92.5))
#: The printed level and the level written for it (SPEC_09 section 2, ruling 1).
LEVELS_MBAR = {"730": 730, "290": 290, "150": 110}
COLUMNS = ("latitude_planetographic_deg", "pressure_mbar", "temperature_K")


@dataclass
class Digitization:
    """The rendered page, the frame, the calibration and every dot found."""

    ink: np.ndarray            #: the page, True where inked
    frame: tuple               #: inner edges (top, bottom, left, right) in page px
    temperature_fit: tuple     #: (slope K/px, intercept K) in page rows
    latitude_fit: tuple        #: (slope deg/px, intercept deg) in page columns
    temperature_residual_K: np.ndarray
    latitude_residual_deg: np.ndarray
    row: np.ndarray            #: dot centers, page px
    col: np.ndarray
    k: np.ndarray              #: dots in the component each center came from
    latitude_deg: np.ndarray
    temperature_K: np.ndarray
    pressure_mbar: np.ndarray
    single_dot_area_px: float
    dropped_curve: int         #: centers dropped by `CURVE_REGION`


def find_ghostscript(given=None) -> str:
    """The Ghostscript executable: `given`, or the first of gswin64c, gswin32c, gs on the path."""
    for name in ([given] if given else ["gswin64c", "gswin32c", "gs"]):
        found = shutil.which(name)
        if found:
            return found
    raise FileNotFoundError(
        "Ghostscript is needed to render the PDF and none was found on the path; pass --gs.")


def render(pdf, gs=None) -> np.ndarray:
    """The page `PAGE` of `pdf` at `DPI`, True where inked."""
    with tempfile.TemporaryDirectory() as directory:
        out = Path(directory) / "page.png"
        subprocess.run([find_ghostscript(gs), "-q", "-dSAFER", "-dBATCH", "-dNOPAUSE",
                        "-sDEVICE=pnggray", f"-r{DPI}", f"-dFirstPage={PAGE}",
                        f"-dLastPage={PAGE}", f"-sOutputFile={out}", str(pdf)], check=True)
        with Image.open(out) as image:
            return np.asarray(image.convert("L")) < 128


def find_frame(ink) -> tuple:
    """Inner edges (top, bottom, left, right) of the plot frame inside `FIGURE_BOX`."""
    x0, y0, x1, y1 = FIGURE_BOX
    box = ink[y0:y1, x0:x1]
    rows = np.flatnonzero(box.sum(axis=1) > 0.6 * box.shape[1])
    cols = np.flatnonzero(box.sum(axis=0) > 0.6 * box.shape[0])
    top_lines, bottom_lines = rows[rows < box.shape[0] // 2], rows[rows >= box.shape[0] // 2]
    left_lines, right_lines = cols[cols < box.shape[1] // 2], cols[cols >= box.shape[1] // 2]
    return (int(top_lines.max()) + 1 + y0, int(bottom_lines.min()) - 1 + y0,
            int(left_lines.max()) + 1 + x0, int(right_lines.min()) - 1 + x0)


def _tick_centers(strip, axis, offset, inner_low, inner_high) -> np.ndarray:
    """Centers of the runs inked across most of `strip`, away from the frame's corners."""
    inked = strip.mean(axis=axis) > 0.5
    labels, count = ndi.label(inked)
    centers = np.array([np.flatnonzero(labels == i).mean() for i in range(1, count + 1)]) + offset
    return centers[(centers > inner_low + 8) & (centers < inner_high - 8)]


def calibrate(ink, frame):
    """The two axis fits, from the left edge's and the bottom edge's ticks."""
    top, bottom, left, right = frame
    rows = _tick_centers(ink[top:bottom + 1, left:left + TICK_STRIP_PX], 1, top, top, bottom)
    cols = _tick_centers(ink[bottom - TICK_STRIP_PX + 1:bottom + 1, left:right + 1], 0, left,
                         left, right)
    if rows.size != len(TEMPERATURE_TICKS_K) or cols.size != len(LATITUDE_TICKS_DEG):
        raise ValueError(f"found {rows.size} temperature ticks and {cols.size} latitude ticks; "
                         f"the figure has {len(TEMPERATURE_TICKS_K)} and "
                         f"{len(LATITUDE_TICKS_DEG)}")
    t_fit = np.polyfit(rows, TEMPERATURE_TICKS_K, 1)
    l_fit = np.polyfit(cols, LATITUDE_TICKS_DEG, 1)
    return (tuple(t_fit), tuple(l_fit),
            np.polyval(t_fit, rows) - np.array(TEMPERATURE_TICKS_K),
            np.polyval(l_fit, cols) - np.array(LATITUDE_TICKS_DEG))


def find_dots(ink, frame):
    """Dot centers in page px, and the dot count of each center's component."""
    top, bottom, left, right = frame
    inner = np.zeros_like(ink)
    inner[top:bottom + 1, left:right + 1] = ink[top:bottom + 1, left:right + 1]
    for x0, y0, x1, y1 in LABEL_BOXES.values():
        inner[y0:y1, x0:x1] = False
    reach = int(np.ceil(OPEN_RADIUS))
    yy, xx = np.mgrid[-reach:reach + 1, -reach:reach + 1]
    opened = ndi.binary_opening(inner, structure=(xx ** 2 + yy ** 2) <= OPEN_RADIUS ** 2)
    labels, count = ndi.label(opened)
    areas = np.bincount(labels.ravel())[1:]
    typical = np.median(areas)
    single = float(np.median(areas[(areas > 0.5 * typical) & (areas < 1.5 * typical)]))
    rows, cols, ks = [], [], []
    # Step 5b: a dot printed narrower than the opening's disk is removed with the strokes. The
    # ink the opening removed, away from what it kept, is taken as one dot where it is the size
    # and shape of one and as thick as a dot, which no stroke of the curve is (`NARROW_DOT`).
    residual = inner & ~ndi.binary_dilation(opened, iterations=2)
    thickness = ndi.distance_transform_edt(inner)
    rest, rest_count = ndi.label(residual)
    for index, region in enumerate(ndi.find_objects(rest), start=1):
        mask = rest[region] == index
        height, width = mask.shape
        area = int(mask.sum())
        if (NARROW_DOT["area_low"] * single <= area <= NARROW_DOT["area_high"] * single
                and max(height, width) <= NARROW_DOT["aspect"] * min(height, width)
                and thickness[region][mask].max() >= NARROW_DOT["thickness_px"]):
            rr, cc = np.nonzero(mask)
            rows.append(rr.mean() + region[0].start)
            cols.append(cc.mean() + region[1].start)
            ks.append(1)
    for index, region in enumerate(ndi.find_objects(labels), start=1):
        rr, cc = np.nonzero(labels[region] == index)
        rr, cc = rr + region[0].start, cc + region[1].start
        k = max(1, int(round(rr.size / single)))
        if k == 1:
            centers = np.array([[rr.mean(), cc.mean()]])
        else:
            centers, _ = kmeans2(np.c_[rr, cc].astype(float), k, minit="++", seed=KMEANS_SEED)
        rows += list(centers[:, 0])
        cols += list(centers[:, 1])
        ks += [k] * k
    return np.array(rows), np.array(cols), np.array(ks), single


def assign_level(latitude, temperature) -> np.ndarray:
    """The level written for each dot, in mbar, by the boundaries of the module docstring."""
    nodes = np.array(BOUNDARY_290_150)
    order = np.argsort(nodes[:, 0])
    boundary = np.interp(latitude, nodes[order, 0], nodes[order, 1])
    return np.where(temperature > BOUNDARY_730_290_K, LEVELS_MBAR["730"],
                    np.where(temperature > boundary, LEVELS_MBAR["290"], LEVELS_MBAR["150"]))


def digitize(pdf, gs=None) -> Digitization:
    """Every step of the module docstring, from the PDF to the dots with their levels."""
    ink = render(pdf, gs)
    frame = find_frame(ink)
    t_fit, l_fit, t_res, l_res = calibrate(ink, frame)
    row, col, k, single = find_dots(ink, frame)
    latitude, temperature = np.polyval(l_fit, col), np.polyval(t_fit, row)
    curve = ((latitude < CURVE_REGION["south_of_deg"])
             & (temperature < CURVE_REGION["colder_than_K"]))
    keep = ~curve
    row, col, k, latitude, temperature = (a[keep] for a in (row, col, k, latitude, temperature))
    return Digitization(ink=ink, frame=frame, temperature_fit=t_fit, latitude_fit=l_fit,
                        temperature_residual_K=t_res, latitude_residual_deg=l_res, row=row,
                        col=col, k=k, latitude_deg=latitude, temperature_K=temperature,
                        pressure_mbar=assign_level(latitude, temperature),
                        single_dot_area_px=single, dropped_curve=int(curve.sum()))


def running_mean(latitude, temperature, at, fwhm_deg=4.0):
    """The Gaussian-weighted running mean of `temperature` in latitude, at the latitudes `at`.

    SPEC_09 acceptance 1. It characterizes the spread of the dots and is not the fit the shear is
    built from, which SPEC_10 defines (REVIEW_09_step1). The weight of a dot is
    `exp(-(at - latitude)^2 / (2 sigma^2))` with `sigma = fwhm_deg / sqrt(8 ln 2)`: the width is
    read as the full width at half maximum, the usual statement of a resolution. Where no dot lies
    within `fwhm_deg` of a latitude the mean is NaN: the fit does not bridge a gap in the data.
    """
    latitude = np.asarray(latitude, dtype="float64")
    temperature = np.asarray(temperature, dtype="float64")
    at = np.atleast_1d(np.asarray(at, dtype="float64"))
    sigma = fwhm_deg / np.sqrt(8.0 * np.log(2.0))
    distance = at[:, None] - latitude[None, :]
    weight = np.exp(-0.5 * (distance / sigma) ** 2)
    mean = (weight @ temperature) / weight.sum(axis=1)
    covered = (np.abs(distance) <= fwhm_deg).any(axis=1)
    return np.where(covered, mean, np.nan)


@dataclass
class Fit:
    """`fit`'s result on its latitude grid. NaN outside the data's own latitude range."""

    latitude_deg: np.ndarray
    value_K: np.ndarray
    gradient_K_per_deg: np.ndarray      #: dT/dphi, planetographic, per degree
    gradient_error_K_per_deg: np.ndarray
    fwhm_deg: np.ndarray                #: the window used: 4 deg, or wider where widened


def _sigma(fwhm_deg):
    return fwhm_deg / np.sqrt(8.0 * np.log(2.0))


def _effective_points(distance, fwhm_deg) -> float:
    weight = np.exp(-0.5 * (distance / _sigma(fwhm_deg)) ** 2)
    return float(weight.sum() ** 2 / (weight ** 2).sum())


def fit(latitude, temperature, degree=1, min_points=3, *, grid, fwhm_deg=4.0) -> Fit:
    """A local polynomial fit of `temperature` in latitude, its value, slope and slope's error.

    SPEC_10 Step 1. The defaults, local linear with at least 3 effective points, are the author's
    choice from the comparison of the three candidates (REVIEW_10_step1). At each latitude `phi0` of `grid` a polynomial of `degree` (1, local linear,
    or 2, local quadratic) in `x = phi - phi0` (planetographic degrees) is fitted by weighted least
    squares with Gaussian weights of full width at half maximum `fwhm_deg`. Where the effective
    number of points `(sum w)^2 / sum w^2` is below `min_points`, the width is widened to the
    smallest value that reaches it, solved continuously so that the fit stays continuous in
    latitude. The value is the polynomial's constant term and the gradient its linear term.

    The gradient's standard error is that of a linear smoother with independent errors of one
    variance: `Cov = s^2 A^-1 (X^T W^2 X) A^-1`, `A = X^T W X`, with `s^2` the weighted mean square
    residual `sum w r^2 / sum w` scaled by `n_eff / (n_eff - (degree + 1))`.

    The fit is defined from the southernmost to the northernmost point; elsewhere on `grid` it is
    NaN (SPEC_10 section 1: beyond the data is SPEC_11's).
    """
    from scipy.optimize import brentq

    latitude = np.asarray(latitude, dtype="float64")
    temperature = np.asarray(temperature, dtype="float64")
    grid = np.asarray(grid, dtype="float64")
    parameters = degree + 1
    if min_points <= parameters:
        raise ValueError(f"min_points = {min_points} leaves no residual degree of freedom for a "
                         f"polynomial of degree {degree}")
    shape = grid.shape
    value, gradient, error, width = (np.full(shape, np.nan) for _ in range(4))
    inside = (grid >= latitude.min()) & (grid <= latitude.max())
    for i in np.flatnonzero(inside):
        x = latitude - grid[i]
        used = fwhm_deg
        if _effective_points(x, used) < min_points:
            high = 2.0 * used
            while _effective_points(x, high) < min_points:
                high *= 2.0
            used = brentq(lambda f: _effective_points(x, f) - min_points, used, high,
                          xtol=1e-9)
        w = np.exp(-0.5 * (x / _sigma(used)) ** 2)
        design = np.vander(x, parameters, increasing=True)
        a = design.T @ (w[:, None] * design)
        a_inv = np.linalg.inv(a)
        beta = a_inv @ (design.T @ (w * temperature))
        residual = temperature - design @ beta
        n_eff = w.sum() ** 2 / (w ** 2).sum()
        s2 = (w @ residual ** 2) / w.sum() * n_eff / (n_eff - parameters)
        covariance = s2 * a_inv @ (design.T @ ((w ** 2)[:, None] * design)) @ a_inv
        value[i], gradient[i] = beta[0], beta[1]
        error[i], width[i] = np.sqrt(covariance[1, 1]), used
    return Fit(latitude_deg=grid, value_K=value, gradient_K_per_deg=gradient,
               gradient_error_K_per_deg=error, fwhm_deg=width)


def write_csv(result: Digitization, path) -> Path:
    """One row per dot, by level from the top down and then from north to south."""
    path = Path(path)
    order = np.lexsort((-result.latitude_deg, result.pressure_mbar))
    lines = [",".join(COLUMNS)]
    lines += [f"{result.latitude_deg[i]:.2f},{int(result.pressure_mbar[i])},"
              f"{result.temperature_K[i]:.2f}" for i in order]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pdf", default="docs/conrath_etal_IRIS_1983.pdf")
    parser.add_argument("--output", default="occul_data/lindal/iris_temperatures.csv")
    parser.add_argument("--gs", default=None, help="the Ghostscript executable")
    args = parser.parse_args(argv)
    result = digitize(args.pdf, args.gs)
    path = write_csv(result, args.output)
    counts = {int(p): int((result.pressure_mbar == p).sum()) for p in sorted(set(result.pressure_mbar))}
    print(f"wrote {path}: {result.latitude_deg.size} dots, by level {counts}; single dot "
          f"{result.single_dot_area_px:g} px; calibration residuals at most "
          f"{np.abs(result.temperature_residual_K).max():.3f} K and "
          f"{np.abs(result.latitude_residual_deg).max():.3f} deg; {result.dropped_curve} "
          "curve remnant(s) dropped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
