"""Kind C read as a field. SPEC_04 Step 0 deliverable 2, decisions L and O.

Manuscript equations implemented: none. This is the one interpolation rule for a composition,
stated once so that every reader uses it.

**Kind C is one structure** (decision O, author, 21 September 2026): a field on
`(level, latitude_planetocentric)` with pressure as the vertical coordinate, uniform in
latitude when that is the hypothesis. The one-latitude form is retired, so a reader that wants
a column asks for it by latitude rather than by the file's shape.

**The rule** (decision L): linear between the file's latitude nodes, the rule the reduction
used to register the anchor, with no extrapolation. A latitude outside the grid is refused,
which is the one refusal decision N leaves here: an interpolant does not invent data it does
not have.

**Flags are not interpolated.** `value_provenance` and the like say how a value was obtained,
and the average of two such codes means nothing, so an integer variable takes the value of the
nearer node. For the uniform field in hand every node carries the same code, so this is a
choice about files that do not exist yet, recorded rather than deferred.
"""

from __future__ import annotations

import numpy as np
import xarray as xr

__all__ = ["LATITUDE_DIM", "column_at", "latitudes_of"]

#: The latitude dimension of kind C (SPEC_00 section 6.2).
LATITUDE_DIM = "latitude_planetocentric"

#: Its coordinate variable.
LATITUDE_NAME = "latitude_planetocentric_deg"


def latitudes_of(composition) -> np.ndarray:
    """The file's planetocentric latitude nodes, in degrees, ascending."""
    root = composition.to_dataset(inherit=False) if hasattr(composition, "to_dataset") \
        else composition
    return np.asarray(root[LATITUDE_NAME].values, dtype="float64")


def column_at(composition, latitude_deg: float):
    """The composition at one planetocentric latitude, as a tree of the same shape.

    Returns a `DataTree` whose root carries every variable on `(level,)`, the latitude
    dimension gone, and whose `species` group is the file's, untouched. `latitude_deg` is in
    degrees, the convention of the file's own coordinate.

    Refuses a latitude outside the file's grid, naming the span: the model does not extrapolate
    a composition (decision N).
    """
    root = composition.to_dataset(inherit=False)
    latitude = latitudes_of(composition)
    order = np.argsort(latitude)
    latitude = latitude[order]
    value = float(latitude_deg)
    if not latitude[0] <= value <= latitude[-1]:
        raise ValueError(
            f"the composition covers latitudes {latitude[0]!r} to {latitude[-1]!r} deg, which "
            f"do not include {value!r} deg; the model does not extrapolate a composition "
            "(SPEC_04 Step 0 deliverable 2)."
        )
    nearest = order[int(np.argmin(np.abs(latitude - value)))]

    out = {}
    for name, variable in root.variables.items():
        if LATITUDE_DIM not in variable.dims:
            out[name] = variable
            continue
        axis = variable.dims.index(LATITUDE_DIM)
        values = np.asarray(variable.values)
        dims = tuple(d for d in variable.dims if d != LATITUDE_DIM)
        if values.dtype.kind in "fc":
            moved = np.moveaxis(values, axis, -1)[..., order]
            picked = np.apply_along_axis(
                lambda row: np.interp(value, latitude, row), -1, moved)
        else:
            picked = np.take(values, nearest, axis=axis)
        out[name] = xr.Variable(dims, np.asarray(picked), variable.attrs)
    if LATITUDE_NAME in out:
        del out[LATITUDE_NAME]

    dataset = xr.Dataset(
        {n: v for n, v in out.items() if n not in root.coords},
        coords={n: v for n, v in out.items() if n in root.coords},
        attrs=dict(root.attrs),
    )
    tree = xr.DataTree(dataset=dataset)
    for name, child in composition.children.items():
        tree[name] = xr.DataTree(dataset=child.to_dataset(inherit=False))
    return tree
