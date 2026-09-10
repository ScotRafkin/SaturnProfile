"""Figure 1 of the CASSPIAN architecture and data-file specification.

Top panel: repository layout as labeled boxes. Bottom panel: data flow.
Bracketed identifiers are the ones the specification text refers to.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

W, H = 16.0, 12.6
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(0, H)
ax.axis("off")

COL = {
    "code": ("#dbe9f6", "#2b5d8c"),
    "static": ("#e6e6e6", "#555555"),
    "profile": ("#dff0d8", "#3c763d"),
    "run": ("#fbe3c8", "#9a5b13"),
    "docs": ("#e8ddf2", "#5e3a87"),
    "root": ("#ffffff", "#000000"),
}
MONO = "DejaVu Sans Mono"
LH = 0.265


def box(x, y, w, h, title, kind, ident=None, lines=(), title_size=11,
        line_size=8.6, lw=1.4, ls="-"):
    face, edge = COL[kind]
    p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08",
                       facecolor=face, edgecolor=edge, linewidth=lw, linestyle=ls)
    ax.add_patch(p)
    head = title if ident is None else f"[{ident}]  {title}"
    ax.text(x + 0.12, y + h - 0.16, head, ha="left", va="top", fontsize=title_size,
            fontweight="bold", family=MONO, color=edge)
    yy = y + h - 0.16 - 0.34
    for ln in lines:
        ax.text(x + 0.18, yy, ln, ha="left", va="top", fontsize=line_size,
                family=MONO, color="#222222")
        yy -= LH


def sub(x, ytop, w, title, kind, ident, lines, ls="-"):
    """Sub box sized to its content; returns its bottom y."""
    h = 0.55 + LH * len(lines)
    box(x, ytop - h, w, h, title, kind, ident, lines, ls=ls)
    return ytop - h


# ---------------------------------------------------------------- header
ax.text(0.3, H - 0.25, "Figure 1. CASSPIAN repository layout (top) and data flow (bottom)",
        fontsize=14, fontweight="bold", va="top")
ax.text(0.3, H - 0.62,
        "Boxes are directories. Bracketed identifiers are the ones referenced in the specification. "
        "Letters T, C, D, G, R, W, N are the netCDF file kinds of specification section 6.",
        fontsize=9.5, va="top", color="#333333")

# ---------------------------------------------------------------- root
box(0.3, H - 1.55, W - 0.6, 0.62, "SaturnProfile/   repository root (git, private GitHub)", "root", "R",
    title_size=12)
ax.text(W - 0.5, H - 1.25, "pyproject.toml   README.md   .gitignore", ha="right", va="center",
        fontsize=8.6, family=MONO, color="#222222")

top = H - 1.75
tree_bottom = 2.55

# ---------------------------------------------------------------- column 1: src
sx, sw = 0.3, 4.55
box(sx, tree_bottom, sw, top - tree_bottom, "src/casspian/   Python package", "code", "S", title_size=11.5)
y = top - 0.4
y = sub(sx + 0.15, y, sw - 0.3, "lib/   pure library", "code", "S1", (
    "constants.py     CODATA physical constants",
    "thermo.py        ideal gas relations, Eq. A9",
    "gravity.py       V, g_N, G_phi: Eqs. A4, A5",
    "geoid.py         Eq. B3; Lindal App. Eqs. 11 to 14",
    "latitude.py      planetographic fixed point",
    "schema.py, io.py netCDF kinds, validators, hashes",
    "control.py       TOML manifest and namelist",
)) - 0.15
y = sub(sx + 0.15, y, sw - 0.3, "tools/   preprocessing", "code", "S2", (
    "lindal/       ascii > raw bundle > T C D G R W",
    "wind/         convert, extend, decompose",
    "              deep parameterization, Eq. A39",
    "composition/  abundances plus species group",
    "gravity/      harmonic set and rotation files",
    "<profile>/    one directory per new source",
)) - 0.15
y = sub(sx + 0.15, y, sw - 0.3, "refrac/   reduction, Appendix B3", "code", "S3", (
    "source agnostic: reads kinds T C D G R W",
    "writes <profile>_refractivity.nc (kind N)",
)) - 0.15
y = sub(sx + 0.15, y, sw - 0.3, "forward/   transfer, B4 to B7", "code", "S4", (
    "reads kind N plus run inputs C G R W (no D)",
    "writes output/ under the run directory",
))

# ---------------------------------------------------------------- column 2: static + docs
dx, dw = 5.15, 3.3
y = top
y = sub(dx, y, dw, "data_static/", "static", "D0", (
    "curated transcriptions; edited",
    "by hand, rarely; never a run input",
    "",
    "species_master      table",
    "harmonics/   Null1981, Iess2019",
    "rotation/    SystemIII, ...",
    "winds/       Smith1982,",
    "             IngersollPollard1982,",
    "             GarciaMelendo2011, ...",
    "thermal/     Fletcher retrievals",
)) - 0.3
y = sub(dx, y, dw, "docs/", "docs", "X", (
    "specs/   SPEC_*.md, this figure,",
    "         STATE.md, REPORT_*.md",
    "manuscript copy (reference only)",
    "numerical methods document",
    "Lindal_et_al_1985_AJ90_1136.pdf",
))

# ---------------------------------------------------------------- column 3: occul_data
ox, ow = 8.75, 3.5
box(ox, tree_bottom, ow, top - tree_bottom, "occul_data/", "profile", "D1", title_size=11.5)
y = top - 0.4
y = sub(ox + 0.15, y, ow - 0.3, "lindal/   one per profile", "profile", "D1a", (
    "raw/  table1.csv, scalars.toml,",
    "      lindal_raw.nc (bundle)",
    "",
    "lindal_thermo.nc          T",
    "lindal_composition.nc     C",
    "lindal_geodesy.nc         D",
    "lindal_gravity.nc         G",
    "lindal_rotation.nc        R",
    "lindal_wind.nc            W",
    "lindal_build.toml         tool control",
    "lindal_reduction.toml     manifest",
    "",
    "lindal_refractivity.nc    N",
    "  phi_c, r0, N(r); embedded",
    "  copies and hashes of T C D G R W",
    "",
    "every file carries the prefix",
)) - 0.25
y = sub(ox + 0.15, y, ow - 0.3, "<profile2>/", "profile", "D1b", (
    "same file set, own prefix",
    "made by tools/<profile2>/",
    "reduced by the same refrac code",
), ls="--")

# ---------------------------------------------------------------- column 4: forward
fx, fw = 12.55, 3.15
box(fx, tree_bottom, fw, top - tree_bottom, "forward/", "run", "D2", title_size=11.5)
y = top - 0.4
y = sub(fx + 0.15, y, fw - 0.3, "<run>/   user named", "run", "D2a", (
    "<run>.toml        namelist",
    "",
    "inputs/",
    "  <run>_composition.nc  C",
    "  <run>_gravity.nc      G",
    "  <run>_rotation.nc     R",
    "  <run>_wind.nc         W",
    "",
    "output/",
    "  fields, profiles,",
    "  diagnostics, resolved",
    "  namelist, input hashes",
    "",
    "points outside itself only to",
    "occul_data/*/*_refractivity.nc",
)) - 0.25
y = sub(fx + 0.15, y, fw - 0.3, "<ensemble>/   later", "run", "D2b", (
    "directory of <run>/ dirs",
    "written by the wrapper, which",
    "is not part of the model",
), ls="--")

# ---------------------------------------------------------------- bottom panel: flow
fy0, fh = 0.3, 1.95
ax.add_patch(FancyBboxPatch((0.3, fy0), W - 0.6, fh, boxstyle="round,pad=0.02,rounding_size=0.08",
                            facecolor="#fafafa", edgecolor="#999999", linewidth=1.0))
ax.text(0.45, fy0 + fh - 0.12, "Data flow", fontsize=11, fontweight="bold", va="top")
ax.text(0.45, fy0 + fh - 0.42,
        "Reduction phase (left of the dashed line) runs once per profile and is then frozen. "
        "Forward phase (right) runs per scenario. The one scalar crossing the line is p_b, declared in the namelist.",
        fontsize=8.6, va="top", color="#333333")

fy = fy0 + 0.25
fbh = 0.78


def fbox(x, w, label, kind, ident=None):
    face, edge = COL[kind]
    ax.add_patch(FancyBboxPatch((x, fy), w, fbh, boxstyle="round,pad=0.02,rounding_size=0.06",
                                facecolor=face, edgecolor=edge, linewidth=1.3))
    head = label if ident is None else f"[{ident}] {label}"
    ax.text(x + w / 2, fy + fbh / 2, head, ha="center", va="center", fontsize=8.6,
            family=MONO, color="#111111", linespacing=1.25)
    return x + w


def farrow(x0, x1, color="#333333"):
    a = FancyArrowPatch((x0 + 0.03, fy + fbh / 2), (x1 - 0.03, fy + fbh / 2), arrowstyle="-|>",
                        mutation_scale=13, linewidth=1.3, color=color)
    ax.add_patch(a)


x = 0.5
e = fbox(x, 1.75, "[D0] data_static\n[D1a] raw/", "static"); farrow(e, e + 0.3); x = e + 0.3
e = fbox(x, 1.55, "[S2] tools/\nlindal, wind, ...", "code"); farrow(e, e + 0.3); x = e + 0.3
e = fbox(x, 1.75, "[D1a] inputs\nT C D G R W + manifest", "profile"); farrow(e, e + 0.3); x = e + 0.3
e = fbox(x, 1.25, "[S3] refrac/", "code"); farrow(e, e + 0.3); x = e + 0.3
e = fbox(x, 1.75, "[D1a] product\n*_refractivity.nc  N", "profile"); farrow(e, e + 0.42); x = e + 0.42
# dashed boundary
bx = x - 0.21
ax.plot([bx, bx], [fy - 0.12, fy + fbh + 0.35], linestyle="--", color="#b00000", linewidth=1.4)
ax.text(bx, fy - 0.16, "reverse-engineering boundary (B1)", ha="center", va="top",
        fontsize=7.8, color="#b00000")
e = fbox(x, 1.4, "[S4] forward/", "code"); farrow(e, e + 0.3); x = e + 0.3
e = fbox(x, 1.75, "[D2a] output/\nT, p, z fields", "run")
# run inputs feeding forward from above
rx = 10.35
ax.add_patch(FancyBboxPatch((rx, fy + fbh + 0.02), 2.0, 0.28, boxstyle="round,pad=0.01",
                            facecolor=COL["run"][0], edgecolor=COL["run"][1], linewidth=1.0))
ax.text(rx + 1.0, fy + fbh + 0.16, "[D2a] namelist + inputs/ C G R W", ha="center", va="center",
        fontsize=7.6, family=MONO)
ax.add_patch(FancyArrowPatch((rx + 1.0, fy + fbh + 0.02), (rx + 1.0, fy + fbh - 0.02) , arrowstyle="-|>",
                             mutation_scale=11, linewidth=1.1, color=COL["run"][1]))
# note: tools also write run inputs
ax.text(W - 0.5, fy0 + 0.06, "tools/ also write the forward inputs under [D2a] inputs/ from [D0].",
        ha="right", va="bottom", fontsize=7.8, color="#333333")

# legend (right of docs box, inside column 2 area below docs)
lx, ly = dx, tree_bottom + 0.05
ax.text(lx, ly + 1.05, "Legend", fontsize=9.5, fontweight="bold", va="bottom")
for i, (k, lab) in enumerate([("code", "source code"), ("static", "static data"),
                              ("profile", "reduction phase data (frozen)"),
                              ("run", "forward phase data (chosen)"), ("docs", "documents")]):
    face, edge = COL[k]
    yy = ly + 0.85 - i * 0.19
    ax.add_patch(FancyBboxPatch((lx, yy - 0.06), 0.28, 0.13, boxstyle="round,pad=0.01",
                                facecolor=face, edgecolor=edge, linewidth=1))
    ax.text(lx + 0.38, yy, lab, fontsize=8.4, va="center")

out = "CASSPIAN_Fig1_RepositoryLayout"
fig.savefig(out + ".svg")
fig.savefig(out + ".png", dpi=170)
print("ok")
