"""What is Lindal's in kind C: the ammonia column of Table I, and the words that go with it.

SPEC_18 v0.3 section 5, separated from `tools/composition/build_composition.py`, which calls this when
a raw bundle's table carries an ammonia column. No other source has one.

Manuscript equations implemented: none.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d

from casspian.lib.control import ControlFileError

#: The uncertainty method of `x_H2_uncertainty`, Lindal's stated range.
H2_UNCERTAINTY_METHOD = "the source's stated 94 +/- 3 percent (p. 1137)"

NH3_ALOFT_ASSUMPTION = (
    "Above the highest tabulated level Table I is blank, which is 'not measured', not "
    "'zero'. Zero is assigned there as an assumption, justified by physics rather than "
    "by the table: Lindal states the troposphere at these levels is saturated with "
    "ammonia, and the saturation mixing ratio over NH3 ice on his own temperatures "
    "falls from a few ppm near 0.8 bar (consistent with the tabulated 2.6 ppm at "
    "831.76 mbar) to below 1e-7 by 0.5 bar and below 1e-12 at the tropopause. The "
    "assumption costs under 1e-6 of the mean molar mass and nothing in the "
    "refractivity, since the lindal1985 set carries NH3 at zero. Levels so filled "
    "carry nh3_provenance = assumed."
)

NH3_NOTE = (
    "The Lindal set carries NH3 refractivity zero. Lindal et al. (1985) Fig. 3 caption "
    "states the density was computed assuming 94 percent hydrogen with the remainder "
    "mostly helium, that is from H2 and He only, so the ammonia contributes no "
    "refractivity in their reduction. Its mole fraction still enters the mean molar "
    "mass."
)


def fill_ammonia(pressure, nh3):
    """Fill the ammonia column under the declared rule, returning values and provenance.

    SPEC_01 v0.17 Step 8, in three parts, because the three regions are known in three
    different ways:

    * **Between the tabulated levels**, linear in pressure: `measured` where Table I gives a
      value, `interpolated` at the one interior gap.
    * **Below the deepest tabulated level**, linear extrapolation from the last two:
      `extrapolated`.
    * **Above the highest tabulated level**, exactly zero: `assumed`. Table I is blank there,
      and blank means not measured, not zero. What justifies zero is physics rather than the
      table, and the reason is written into the file as `nh3_aloft_assumption`.

    The v0.5 wording of the rule extrapolated upward and clamped at zero, which put 0.2 ppm at
    794.33 mbar and made the clamp level a puzzle (REPORT_01_step8, finding 1). There is no
    upward extrapolation and no clamp here: a value nobody measured is an assumption, and it is
    labelled as one.
    """
    pressure = np.asarray(pressure, dtype="float64")
    nh3 = np.asarray(nh3, dtype="float64")
    tabulated = np.isfinite(nh3)
    if tabulated.sum() < 2:
        raise ControlFileError("the ammonia column needs at least two tabulated values")

    lo, hi = float(pressure[tabulated].min()), float(pressure[tabulated].max())
    rule = interp1d(pressure[tabulated], nh3[tabulated], kind="linear",
                    bounds_error=False, fill_value="extrapolate")

    filled = np.asarray(rule(pressure), dtype="float64")
    aloft = pressure < lo
    filled[aloft] = 0.0

    provenance = np.where(
        aloft, "assumed",
        np.where(tabulated, "measured",
                 np.where(pressure <= hi, "interpolated", "extrapolated")),
    )
    return filled, provenance, lo, hi
