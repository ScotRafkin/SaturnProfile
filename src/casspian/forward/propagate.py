"""The seasonal propagation hook. SPEC_04 Step 0 deliverable 7, decision J.

Manuscript equations implemented: none. Eq. A36 states that the model relates latitudes
analytically and epochs not at all, so a season is crossed by moving the anchor itself, in `N`
at fixed pressure at its own latitude, before it is transferred (SPEC_04 decision I, the
seasonal design note (S8), (S9), (S9b) and section 5.4).

**This is the one place an anchor is moved in season, and here it is the identity.** It is
called for every anchor at every M and every season, between the loader and Step 1's placement,
so the code path is the same whether or not a propagator exists. The propagator specification
replaces the body of `propagate` and nothing that calls it: the object it returns and the two
uncertainty columns on it are the interface, and they do not change.

What it does today: returns the anchor with `N` unchanged and its season term zero, and records
the anchor's season, the run's, and that no propagation was applied. The absence is recorded in
the product rather than assumed, because a season term that is silently zero is a term that is
trusted once forgotten (decision J).
"""

from __future__ import annotations

import dataclasses
from types import MappingProxyType

__all__ = ["PROPAGATION_ABSENT", "propagate"]

#: What the record says while the propagator does not exist. The product carries it verbatim.
PROPAGATION_ABSENT = "none: propagator not implemented"


def propagate(anchor, season: float):
    """Move `anchor` to `season`. SPEC_04 Step 0 deliverable 7: today, the identity.

    `anchor` is a `casspian.lib.control.LoadedAnchor` and `season` the run's declared solar
    longitude in degrees. Returns a `LoadedAnchor` with its `propagation` record filled: the
    anchor's own season, the run's, the seasonal distance where both are known, and what was
    applied. `N`, the levels and the season term are untouched, which is what makes this the
    identity; a run whose anchors are all at the run's season passes through it the same way as
    one whose anchors are not.
    """
    own = anchor.anchor_season_deg
    record = {
        "anchor_solar_longitude_deg": "uniform" if own is None else float(own),
        "run_solar_longitude_deg": float(season),
        "propagation": PROPAGATION_ABSENT,
        "season_term": anchor.season_term,
    }
    if own is not None:
        # The seasonal circle: the shorter way round, in [-180, 180) degrees.
        distance = (float(season) - float(own) + 180.0) % 360.0 - 180.0
        record["seasonal_distance_deg"] = distance
        record["season_matches_run"] = bool(distance == 0.0)
    else:
        record["seasonal_distance_deg"] = "not defined: the anchor declares itself uniform"
        record["season_matches_run"] = False
    return dataclasses.replace(anchor, propagation=MappingProxyType(record))
