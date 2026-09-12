"""The standard diagnostics. SPEC_02 Step 5; SPEC_00 section 3.5.

A generic tool, tied to no profile, that renders a fixed set of figures from a CASSPIAN file
alone, dispatching on `casspian_kind`. Figures are made from the file, never from a run's
memory, so any product ever written can be drawn again. Nothing here computes a number that
reaches a product.
"""

from casspian.tools.plots.render import (
    SUPPORTED_KINDS,
    RenderResult,
    UnsupportedKindError,
    describe,
    main,
    render,
)

__all__ = ["SUPPORTED_KINDS", "RenderResult", "UnsupportedKindError", "describe", "main", "render"]
