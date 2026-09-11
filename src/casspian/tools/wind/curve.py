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
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.interpolate import PchipInterpolator

__all__ = ["Segments", "read_curve", "AssembledCurve"]

#: `value_provenance` codes, SPEC_00 section 6.6.
OBSERVED, INTERPOLATED, PARAMETERIZED, EXTRAPOLATED, EXTENDED = 0, 1, 2, 3, 4


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
