"""`render` and the `casspian-plots` entry point. SPEC_02 Step 5.

`render` reads `casspian_kind` and dispatches: kind N renders the figure set F1 to F6 (those its
fields allow); a kind W, T or C file renders its single figure view; any other kind is refused
by name. Figures go to `<out_dir>/<prefix>_diag_<key>_<name>.<format>`, with a combined PDF of
the kind N set at `<out_dir>/<prefix>_diag.pdf`. The default `out_dir` is `figures/` beside the
file.
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import netCDF4
from matplotlib.backends.backend_pdf import PdfPages

from casspian.lib import io as cio
from casspian.lib.schema import CasspianSchemaError
from casspian.tools.plots import figures_inputs, figures_product, style

TOOL = "casspian-plots"
SUPPORTED_KINDS = ("refractivity", "wind", "thermo", "composition")
FORMATS = ("png", "pdf")


class UnsupportedKindError(ValueError):
    """The file's kind has no diagnostic figures."""


@dataclass
class RenderResult:
    """What was written, what was skipped and why, and the arrays each figure drew."""

    path: Path
    kind: str
    written: list = field(default_factory=list)
    combined_pdf: Path | None = None
    skipped: dict = field(default_factory=dict)
    footers: dict = field(default_factory=dict)
    data: dict = field(default_factory=dict)
    generated_at: str = ""


def _declared_kind(path: Path) -> str:
    with netCDF4.Dataset(path, "r") as handle:
        if "casspian_kind" not in handle.ncattrs():
            raise CasspianSchemaError(f"{path.name}: no casspian_kind attribute, so this is not a "
                                      "CASSPIAN file (SPEC_00 section 5).")
        return str(handle.getncattr("casspian_kind"))


def render(path, out_dir=None, format: str = "png", dpi: int = 150,
           generated_at: str | None = None) -> RenderResult:
    """Render the standard figures of one CASSPIAN file. See the module docstring."""
    path = Path(path)
    kind = _declared_kind(path)
    if kind not in SUPPORTED_KINDS:
        raise UnsupportedKindError(
            f"{path.name}: casspian_kind {kind!r} has no diagnostic figures. casspian-plots "
            "renders kind N (refractivity) and the single figure views of kinds W, T and C "
            "(SPEC_02 Step 5)."
        )
    if format not in FORMATS:
        raise ValueError(f"format {format!r}; SPEC_02 Step 5 allows {list(FORMATS)}")
    out_dir = Path(out_dir) if out_dir is not None else path.parent / "figures"
    stamp = generated_at or datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    result = RenderResult(path=path, kind=kind, generated_at=stamp)

    handle = cio.read(path, kind)
    try:
        prefix = str(handle.attrs["profile_or_run"])
        commit = str(handle.attrs["casspian_git_commit"])
        sha12 = cio.sha256(path)[:12]
        with style.styled():
            if kind == "refractivity":
                figures, skipped, data = figures_product.build_all(handle)
            else:
                figures, skipped, data = figures_inputs.build(kind, handle)
            result.skipped, result.data = skipped, data
            out_dir.mkdir(parents=True, exist_ok=True)
            for key, name, fig in figures:
                result.footers[key] = style.footer(fig, path.name, commit, sha12, stamp)
                stem = f"{prefix}_diag_{key}_{name}" if key != name else f"{prefix}_diag_{kind}"
                target = out_dir / f"{stem}.{format}"
                fig.savefig(target, dpi=dpi, format=format)
                result.written.append(target)
            if kind == "refractivity":
                combined = out_dir / f"{prefix}_diag.pdf"
                with PdfPages(combined) as pdf:
                    for _, _, fig in figures:
                        pdf.savefig(fig, dpi=dpi)
                result.combined_pdf = combined
    finally:
        handle.close()
    return result


def describe(result: RenderResult) -> list[str]:
    """The lines a caller prints: each figure written, the combined PDF, each figure skipped."""
    lines = [f"figure {p}" for p in result.written]
    if result.combined_pdf is not None:
        lines.append(f"combined {result.combined_pdf}")
    lines += [f"skipped {key}: {reason}" for key, reason in result.skipped.items()]
    return lines


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog=TOOL, description="Render the standard diagnostic figures of a CASSPIAN file.")
    parser.add_argument("file", help="a CASSPIAN netCDF file")
    parser.add_argument("--out", default=None, help="output directory (default: figures/ beside "
                                                    "the file)")
    parser.add_argument("--format", default="png", choices=FORMATS)
    parser.add_argument("--dpi", type=int, default=150)
    args = parser.parse_args(argv)
    try:
        result = render(args.file, args.out, args.format, args.dpi)
    except (UnsupportedKindError, CasspianSchemaError) as exc:
        print(f"{TOOL}: refused: {exc}", file=sys.stderr)
        return 2
    for line in describe(result):
        print(line)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
