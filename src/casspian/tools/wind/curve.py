"""Assembly of a published zonal wind curve into a function of latitude on the whole globe.

SPEC_01 v0.15 Step 7, items 1 to 3. The wind is the solid curve of Ingersoll and Pollard
(1982) Fig. 5, digitized in two segments with a gap where the rings obscured the southern
equatorial latitudes and with ends short of both poles. This module turns those samples into
`u(phi_g)` on [-90, 90] under declared rules, and reports which rule produced every value.

Manuscript equations implemented: none. This is interpolation and declared extension.

Every rule here is named in the control file and flagged in `value_provenance`. Nothing is
smoothed, fitted or adjusted: the curve passes through its samples, and the regions the curve
does not cover carry the rule that filled them.

The three regions the samples do not cover:

* **The ring gap**, -10.9 to +1.3 degrees. Filled by reflecting the northern segment about
  the equator, which is the authors' own dashed curve (their Fig. 3 caption). The strip
  between the northern end and its mirror is bridged by a cubic Hermite matched to value and
  slope at both ends, which by the symmetry is even and reduces to `a + b phi^2`.
* **The join**, where the reflection meets the southern segment. They do not meet: at -10.9
  the reflection reads about 377 m/s and the southern solid curve 346. The two are blended
  linearly across a declared window so the assembled wind is continuous.
* **The polar caps**, poleward of +81.1 and -72.8. Brought to exactly zero at the poles,
  which SPEC_00 section 6.6 requires of every kind W file.

**A tabulated wind** (SPEC_13 Step 1 item 2, `curve_format = "table"`): the Sanchez-Lavega et al.
(2000) table, one row per latitude bin with the wind and the standard deviation of the
measurements in the bin (`u_rms`, their Eq. 6). `read_table` reads it as one segment, dropping the
rows the file marks as not data; `TableCurve` assembles it under `gap_rule = "pchip_bridge"`: one
PCHIP through every row, so the ring gap and the small gaps are bridged between their data edges
with no reflection and no join, and the poles by the polar rule. Its uncertainty is the table's
`u_rms`, linear in latitude, except in the declared ring gap, where the scatter was never measured
and it is NaN (SPEC_13 section 3 ruling 2), and poleward of the rows, as for the curve.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ["Segments", "read_curve", "AssembledCurve", "Table", "read_table", "TableCurve"]

#: `value_provenance` codes, SPEC_00 section 6.6.
OBSERVED, INTERPOLATED, PARAMETERIZED, EXTRAPOLATED, EXTENDED = 0, 1, 2, 3, 4

#: A tabulated wind: the marker of a row the file itself declares is not data.
NOT_DATA_MARKER = "fake data"
#: A tabulated wind: an interval between adjacent rows wider than this many times the median
#: interval is a gap, and the latitudes inside it are `interpolated` rather than `observed`. The
#: rows of the Sanchez-Lavega table are 0.5 deg apart planetocentric, 0.45 to 0.62 deg
#: planetographic; its gaps are 1.5 deg and wider.
GAP_FACTOR = 2.0


@dataclass(frozen=True)
class Segments:
    """The digitized curve, one entry per solid segment, each sorted by latitude."""

    north: tuple
    south: tuple

    @property
    def north_range(self):
        return float(self.north[0].min()), float(self.north[0].max())

    @property
    def south_range(self):
        return float(self.south[0].min()), float(self.south[0].max())


def read_curve(path) -> Segments:
    """Read the digitized curve CSV into its named segments."""
    with open(Path(path), newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    by_segment = {}
    for row in rows:
        by_segment.setdefault(row["segment"], []).append(
            (float(row["latitude_planetographic_deg"]), float(row["u_ms"]))
        )
    missing = {"solid_north", "solid_south"} - set(by_segment)
    if missing:
        raise ValueError(f"{path}: the curve is missing the segment(s) {sorted(missing)}")
    out = {}
    for name, pairs in by_segment.items():
        array = np.array(sorted(pairs), dtype="float64")
        out[name] = (array[:, 0], array[:, 1])
    return Segments(north=out["solid_north"], south=out["solid_south"])


def _polar_piece(x_end, y_end, x_prev, y_prev, pole, rule):
    """The cap from a segment end to exactly zero at `pole`, under the declared rule."""
    if rule == "linear_to_zero":
        xs = np.array([x_end, pole]) if pole > 0 else np.array([pole, x_end])
        ys = np.array([y_end, 0.0]) if pole > 0 else np.array([0.0, y_end])
        return lambda x: np.interp(x, xs, ys)
    if rule != "pchip_to_zero":
        raise ValueError(f"unknown polar_rule {rule!r}")
    # PCHIP through the last two digitized samples and the pole. Shape preserving, so it
    # cannot overshoot or turn back on its way to zero.
    pairs = sorted([(x_prev, y_prev), (x_end, y_end), (pole, 0.0)])
    xs = np.array([p[0] for p in pairs], dtype="float64")
    ys = np.array([p[1] for p in pairs], dtype="float64")
    interpolant = PchipInterpolator(xs, ys, extrapolate=False)
    return lambda x: np.nan_to_num(interpolant(np.clip(x, xs[0], xs[-1])), nan=0.0)


class AssembledCurve:
    """`u(phi_g)` on [-90, 90], with the rule that produced every value."""

    def __init__(self, segments: Segments, *, gap_rule: str, polar_rule: str,
                 join_window_deg, bin_centers=None, bin_values=None):
        self.segments = segments
        self.gap_rule = gap_rule
        self.polar_rule = polar_rule
        self.join_hi, self.join_lo = (float(join_window_deg[0]), float(join_window_deg[1]))
        if self.join_lo > self.join_hi:
            self.join_hi, self.join_lo = self.join_lo, self.join_hi
        self.bin_centers = bin_centers
        self.bin_values = bin_values

        north_lat, north_u = segments.north
        south_lat, south_u = segments.south
        self.north_min, self.north_max = float(north_lat[0]), float(north_lat[-1])
        self.south_min, self.south_max = float(south_lat[0]), float(south_lat[-1])

        self._north = PchipInterpolator(north_lat, north_u, extrapolate=False)
        self._south = PchipInterpolator(south_lat, south_u, extrapolate=False)

        # The even bridge across the strip: value and slope matched at +north_min, mirrored
        # at -north_min, so a + b phi^2 with b from the slope.
        edge_value = float(self._north(self.north_min))
        edge_slope = float(self._north(self.north_min, nu=1))
        self._bridge_b = edge_slope / (2.0 * self.north_min)
        self._bridge_a = edge_value - self._bridge_b * self.north_min**2

        self._north_cap = _polar_piece(self.north_max, float(north_u[-1]),
                                       float(north_lat[-2]), float(north_u[-2]),
                                       90.0, polar_rule)
        self._south_cap = _polar_piece(self.south_min, float(south_u[0]),
                                       float(south_lat[1]), float(south_u[1]),
                                       -90.0, polar_rule)

    # -- the pieces -------------------------------------------------------------------
    def _reflected(self, phi):
        return np.nan_to_num(self._north(-phi), nan=0.0)

    def _bridge(self, phi):
        return self._bridge_a + self._bridge_b * phi**2

    def _gap_value(self, phi):
        """The fill inside the ring gap, before the join blend."""
        inside_strip = np.abs(phi) <= self.north_min
        return np.where(inside_strip, self._bridge(phi), self._reflected(phi))

    def _bins_in_gap(self, phi):
        """The Smith bin means, for the `reflect_north_then_bins` alternative."""
        if self.bin_centers is None:
            raise ValueError("gap_rule 'reflect_north_then_bins' needs the Smith bin means")
        good = np.isfinite(self.bin_values)
        return np.interp(phi, self.bin_centers[good], self.bin_values[good],
                         left=np.nan, right=np.nan)

    # -- the whole curve --------------------------------------------------------------
    def __call__(self, phi):
        phi = np.asarray(phi, dtype="float64")
        out = np.empty(phi.shape, dtype="float64")

        north = (phi >= self.north_min) & (phi <= self.north_max)
        south = (phi >= self.south_min) & (phi < self.join_lo)
        blend = (phi >= self.join_lo) & (phi < self.join_hi)
        gap = (phi >= self.join_hi) & (phi < self.north_min)
        north_cap = phi > self.north_max
        south_cap = phi < self.south_min

        out[north] = np.nan_to_num(self._north(phi[north]), nan=0.0)
        out[south] = np.nan_to_num(self._south(phi[south]), nan=0.0)

        if blend.any():
            x = phi[blend]
            weight = np.clip((x - self.join_lo) / (self.join_hi - self.join_lo), 0.0, 1.0)
            reflection = self._gap_value(x)
            if self.gap_rule == "reflect_north_then_bins":
                from_bins = self._bins_in_gap(x)
                reflection = np.where(np.isfinite(from_bins), from_bins, reflection)
            out[blend] = weight * reflection + (1.0 - weight) * np.nan_to_num(
                self._south(x), nan=0.0
            )

        if gap.any():
            x = phi[gap]
            value = self._gap_value(x)
            if self.gap_rule == "reflect_north_then_bins":
                from_bins = self._bins_in_gap(x)
                value = np.where(np.isfinite(from_bins), from_bins, value)
            out[gap] = value

        out[north_cap] = self._north_cap(phi[north_cap])
        out[south_cap] = self._south_cap(phi[south_cap])
        # The poles are exactly zero, not nearly: SPEC_00 section 6.6 is checked with ==.
        out[np.abs(phi) >= 90.0] = 0.0
        return out

    def provenance(self, phi):
        """The `value_provenance` code for every latitude, SPEC_01 v0.15 Step 7 item 5."""
        phi = np.asarray(phi, dtype="float64")
        code = np.full(phi.shape, OBSERVED, dtype="int8")
        code[(phi >= self.join_lo) & (phi < self.north_min)] = PARAMETERIZED
        code[phi > self.north_max] = EXTRAPOLATED
        code[phi < self.south_min] = EXTRAPOLATED
        code[np.abs(phi) >= 90.0] = EXTRAPOLATED
        return code

    def uncertainty(self, phi, centers, rms, counts, min_count):
        """1 sigma at each latitude, from the Smith bins, per SPEC_01 v0.15 Step 7 item 6."""
        phi = np.asarray(phi, dtype="float64")
        usable = (counts >= min_count) & np.isfinite(rms)
        value = np.interp(phi, centers[usable], rms[usable], left=np.nan, right=np.nan)

        # In the gap, where the reflection supplies the wind, the uncertainty is mirrored
        # from the north; where Smith bins exist (the southern flank), theirs is used.
        gap = (phi >= self.join_hi) & (phi < self.north_min)
        if gap.any():
            here = phi[gap]
            direct = np.interp(here, centers[usable], rms[usable], left=np.nan, right=np.nan)
            mirrored = np.interp(-here, centers[usable], rms[usable],
                                 left=np.nan, right=np.nan)
            value[gap] = np.where(np.isfinite(direct), direct, mirrored)

        # No information poleward of the data, and none at the poles.
        value[phi > self.north_max] = np.nan
        value[phi < self.south_min] = np.nan
        value[np.abs(phi) >= 90.0] = np.nan
        return value


@dataclass(frozen=True)
class Table:
    """A tabulated wind, sorted by planetographic latitude, with the rows dropped as not data."""

    latitude_deg: np.ndarray        #: planetographic
    u_ms: np.ndarray
    u_rms_ms: np.ndarray
    dropped: int
    header: str                     #: the file's first line, its own statement of its source


def read_table(path) -> Table:
    """The Sanchez-Lavega table: two header lines, then planetocentric latitude, planetographic
    latitude, `u` and `u_rms` in whitespace-separated columns. A row carrying `NOT_DATA_MARKER`
    in its trailing comment is dropped."""
    rows, dropped = [], 0
    with open(Path(path), encoding="utf-8") as handle:
        lines = handle.read().splitlines()
    for line in lines[2:]:
        if not line.strip():
            continue
        values, _, comment = line.partition("!")
        if NOT_DATA_MARKER in comment:
            dropped += 1
            continue
        fields = values.split()
        if len(fields) != 4:
            raise ValueError(f"{path}: a row has {len(fields)} columns, not 4: {line!r}")
        rows.append((float(fields[1]), float(fields[2]), float(fields[3])))
    array = np.array(sorted(rows), dtype="float64")
    if np.any(np.diff(array[:, 0]) <= 0.0):
        raise ValueError(f"{path}: the planetographic latitudes are not distinct")
    return Table(latitude_deg=array[:, 0], u_ms=array[:, 1], u_rms_ms=array[:, 2], dropped=dropped,
                 header=lines[0].strip())


class TableCurve:
    """`u(phi_g)` on [-90, 90] from a tabulated wind, with the rule that produced every value.

    `ring_gap_deg` is the declared ring gap, whose two ends must be rows of the table."""

    def __init__(self, table: Table, *, gap_rule: str, polar_rule: str, ring_gap_deg):
        if gap_rule != "pchip_bridge":
            raise ValueError(f"gap_rule {gap_rule!r} is not a rule for a tabulated wind; "
                             "the rule is 'pchip_bridge'")
        self.table = table
        self.gap_rule = gap_rule
        self.polar_rule = polar_rule
        lat, u = table.latitude_deg, table.u_ms
        self.ring_lo, self.ring_hi = sorted(float(v) for v in ring_gap_deg)
        for edge in (self.ring_lo, self.ring_hi):
            if not np.any(lat == edge):
                raise ValueError(f"the ring gap's edge {edge} deg is not a row of the table")
        inside = (lat > self.ring_lo) & (lat < self.ring_hi)
        if inside.any():
            raise ValueError(f"the table has rows inside the declared ring gap: {lat[inside]}")
        self.south_min, self.north_max = float(lat[0]), float(lat[-1])
        self._u = PchipInterpolator(lat, u, extrapolate=False)
        spacing = np.diff(lat)
        self.gaps = [(float(lat[i]), float(lat[i + 1]))
                     for i in np.flatnonzero(spacing > GAP_FACTOR * np.median(spacing))]
        self._north_cap = _polar_piece(self.north_max, float(u[-1]), float(lat[-2]), float(u[-2]),
                                       90.0, polar_rule)
        self._south_cap = _polar_piece(self.south_min, float(u[0]), float(lat[1]), float(u[1]),
                                       -90.0, polar_rule)

    def __call__(self, phi):
        phi = np.asarray(phi, dtype="float64")
        out = np.empty(phi.shape, dtype="float64")
        rows = (phi >= self.south_min) & (phi <= self.north_max)
        out[rows] = self._u(phi[rows])
        north_cap, south_cap = phi > self.north_max, phi < self.south_min
        out[north_cap] = self._north_cap(phi[north_cap])
        out[south_cap] = self._south_cap(phi[south_cap])
        out[np.abs(phi) >= 90.0] = 0.0
        return out

    def _in_gap(self, phi):
        inside = np.zeros(phi.shape, dtype=bool)
        for lo, hi in self.gaps:
            inside |= (phi > lo) & (phi < hi)
        return inside

    def provenance(self, phi):
        """`observed` among the rows, `interpolated` inside a gap, `extrapolated` beyond the rows
        (SPEC_13 section 3 ruling 2)."""
        phi = np.asarray(phi, dtype="float64")
        code = np.full(phi.shape, OBSERVED, dtype="int8")
        code[self._in_gap(phi)] = INTERPOLATED
        code[(phi > self.north_max) | (phi < self.south_min) | (np.abs(phi) >= 90.0)] = EXTRAPOLATED
        return code

    def uncertainty(self, phi):
        """`u_rms` linear in latitude; NaN in the ring gap and beyond the rows."""
        phi = np.asarray(phi, dtype="float64")
        value = np.interp(phi, self.table.latitude_deg, self.table.u_rms_ms, left=np.nan,
                          right=np.nan)
        value[(phi > self.ring_lo) & (phi < self.ring_hi)] = np.nan
        value[np.abs(phi) >= 90.0] = np.nan
        return value
