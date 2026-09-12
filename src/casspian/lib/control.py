"""Control files: the reduction manifest, its inputs, and the sections the tools read.

SPEC_00 sections 7 and 8; SPEC_02 Step 1. A control file holds pointers and declared choices
only. A physical value in a control file is a defect, and an unknown key is an error rather than
a warning, because the parser knows which keys exist. A relative path is resolved against the
directory holding the control file.

Manuscript equations implemented: none. This module reads, checks and refuses; it computes
nothing that reaches a product.

Two layers. `load_section`, `reject_physical_values` and `ControlFileError` are the narrow
section reader the tools used from `tools/gravity/control.py`, absorbed here unchanged in
behavior (SPEC_02 Step 1), with its named-key path resolution kept. `read_reduction_manifest`
and `load_reduction_inputs` are the reduction's: the manifest parsed against the SPEC_00
section 7.1 vocabulary, then the six inputs read under their kinds and checked against each
other before any arithmetic.

**Units at this boundary.** The manifest carries `fixed_point_tolerance_deg` in degrees,
because it is a file a person edits. `lib` is radians without exception (Step 6 review), and
SPEC_02 section 0 has the conversion happen once, at the `refrac` boundary, and nowhere else.
So this module returns the tolerance in degrees as written and deliberately offers no radian
form: a second place that converts would be a second place to get it wrong.
"""

from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

import numpy as np

from casspian.lib import io as cio
from casspian.lib.schema import UNIT_SUFFIXES

__all__ = [
    "ControlFileError",
    "load_section",
    "reject_physical_values",
    "ReductionManifest",
    "ReductionInputs",
    "read_reduction_manifest",
    "load_reduction_inputs",
    "INPUT_KINDS",
]


class ControlFileError(Exception):
    """A control file, or the inputs it names, fails a rule of SPEC_00 section 7 or 8.

    The message always names the file and the rule.
    """


# ---------------------------------------------------------------------------
# The section reader the tools use (absorbed from tools/gravity/control.py)
# ---------------------------------------------------------------------------


def load_section(path, section: str, allowed: dict[str, bool], path_keys=()) -> dict:
    """Load one section of a TOML control file and check its keys against `allowed`.

    `allowed` maps key name to whether it is required. Any key outside `allowed` is an error
    (SPEC_00 section 7).

    `path_keys` names the keys whose values are paths; each is resolved against the control
    file's own directory and returned as an absolute `Path`. The keys are named explicitly
    rather than guessed from the file extension, because guessing was wrong the first time a
    tool pointed at a `.csv`.
    """
    path = Path(path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: control file does not exist")
    with open(path, "rb") as handle:
        document = tomllib.load(handle)

    if section not in document:
        raise ControlFileError(
            f"{path}: no [{section}] section, which the tool reading this file requires."
        )
    table = document[section]

    unknown = sorted(set(table) - set(allowed))
    if unknown:
        raise ControlFileError(
            f"{path}: [{section}] carries unknown key(s) {unknown}. SPEC_00 section 7 makes an "
            f"unknown key an error. The keys this section accepts are {sorted(allowed)}."
        )
    missing = sorted(name for name, required in allowed.items() if required and name not in table)
    if missing:
        raise ControlFileError(f"{path}: [{section}] is missing required key(s) {missing}.")

    resolved = {}
    for key, value in table.items():
        if key in set(path_keys) and isinstance(value, str):
            resolved[key] = (path.parent / value).resolve()
        else:
            resolved[key] = value
    return resolved


def reject_physical_values(table: dict, path, section: str) -> None:
    """Refuse a control file section that carries a bare number.

    SPEC_00 principle 2: data lives in netCDF, choices live in TOML. For sections whose every
    legitimate value is a path or a name, any number is a physical value that has strayed.
    """
    numeric = sorted(
        key for key, value in table.items()
        if isinstance(value, (int, float)) and not isinstance(value, bool)
    )
    if numeric:
        raise ControlFileError(
            f"{path}: [{section}] carries the numeric value(s) {numeric}. SPEC_00 principle 2 "
            "keeps physical values in data files; a control file holds pointers and choices."
        )


# ---------------------------------------------------------------------------
# The reduction manifest, SPEC_00 section 7.1
# ---------------------------------------------------------------------------

#: The six inputs a reduction reads, by manifest key and the kind each must validate as.
INPUT_KINDS = MappingProxyType({
    "thermo": "thermo",
    "composition": "composition",
    "geodesy": "geodesy",
    "gravity": "gravity",
    "rotation": "rotation",
    "wind": "wind",
})

_NUMBER = "number"
_INTEGER = "integer"
_TEXT = "text"
_FLAG = "flag"

#: SPEC_00 section 7.1 (v0.11). Section name -> (section required, {key: (type, required)}).
#: Inside the two optional sections the keys a switched-on feature needs are enforced in
#: `_check_optional_sections`, not here, because their requirement depends on the switch.
_VOCABULARY = {
    "profile": (True, {"slug": (_TEXT, True), "description": (_TEXT, True)}),
    "inputs": (True, {key: (_TEXT, True) for key in INPUT_KINDS}),
    "anchor_isobar": (True, {"pressure_Pa": (_NUMBER, True)}),
    "anchor_latitude": (True, {
        "fixed_point_tolerance_deg": (_NUMBER, True),
        "max_iterations": (_INTEGER, True),
    }),
    "geoid": (True, {
        "anchor_surface_Pa": (_NUMBER, True),
        "anchor_quantity": (_TEXT, True),
        "anchor_rule": (_TEXT, True),
        "convergence_m": (_NUMBER, True),
    }),
    "output": (True, {"product": (_TEXT, True)}),
    "diagnostics": (False, {
        "figures": (_FLAG, False),
        "format": (_TEXT, False),
        "dpi": (_INTEGER, False),
    }),
    "sensitivity": (False, {
        "enabled": (_FLAG, False),
        "wind_draws": (_INTEGER, False),
        "wind_draw_seed": (_INTEGER, False),
    }),
}

#: The anchor rules a manifest may name. `lib.geoid.wind_geoid` also knows "latitude", but that
#: rule needs an anchor latitude and the section 7.1 vocabulary has no key to carry one.
_MANIFEST_ANCHOR_RULES = frozenset({"mean_polar_radius", "north_pole", "south_pole"})
_DIAGNOSTIC_FORMATS = frozenset({"png", "pdf"})


@dataclass(frozen=True)
class ReductionManifest:
    """A parsed, checked reduction manifest. Paths are absolute; the text is kept verbatim."""

    path: Path
    sha256: str
    text: str
    slug: str
    description: str
    inputs: MappingProxyType
    anchor_isobar_Pa: float
    fixed_point_tolerance_deg: float
    max_iterations: int
    anchor_surface_Pa: float
    anchor_quantity: str
    anchor_rule: str
    convergence_m: float
    product: Path
    diagnostics: MappingProxyType
    sensitivity: MappingProxyType


def _is_type(value, kind: str) -> bool:
    if kind == _NUMBER:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if kind == _INTEGER:
        return isinstance(value, int) and not isinstance(value, bool)
    if kind == _FLAG:
        return isinstance(value, bool)
    return isinstance(value, str)


def _looks_physical(key: str, value) -> bool:
    """An unknown key naming a quantity with a unit, holding a number, is data in a control file."""
    numeric = isinstance(value, (int, float)) and not isinstance(value, bool)
    return numeric and any(key.endswith(suffix) for suffix in UNIT_SUFFIXES)


def _check_positive(path, section, key, value, integer=False):
    if (integer and value < 1) or (not integer and not value > 0):
        raise ControlFileError(
            f"{path}: [{section}] {key} = {value!r} must be "
            f"{'a positive integer' if integer else 'positive'}."
        )


def _check_optional_sections(path, document) -> tuple[dict, dict]:
    diagnostics = dict(document.get("diagnostics", {}))
    diagnostics.setdefault("figures", False)       # SPEC_00 section 7.1: default false
    if diagnostics["figures"]:
        for key in ("format", "dpi"):
            if key not in diagnostics:
                raise ControlFileError(
                    f"{path}: [diagnostics] figures = true requires {key!r}; no default is "
                    "declared for it."
                )
        if diagnostics["format"] not in _DIAGNOSTIC_FORMATS:
            raise ControlFileError(
                f"{path}: [diagnostics] format = {diagnostics['format']!r}; SPEC_02 Step 5 "
                f"allows {sorted(_DIAGNOSTIC_FORMATS)}."
            )
        _check_positive(path, "diagnostics", "dpi", diagnostics["dpi"], integer=True)

    sensitivity = dict(document.get("sensitivity", {}))
    sensitivity.setdefault("enabled", False)        # SPEC_00 section 7.1: default false
    if sensitivity["enabled"]:
        for key in ("wind_draws", "wind_draw_seed"):
            if key not in sensitivity:
                raise ControlFileError(
                    f"{path}: [sensitivity] enabled = true requires {key!r}; no default is "
                    "declared for it."
                )
        _check_positive(path, "sensitivity", "wind_draws", sensitivity["wind_draws"],
                        integer=True)
    return diagnostics, sensitivity


def read_reduction_manifest(path) -> ReductionManifest:
    """Parse and check `<profile>_reduction.toml` against SPEC_00 section 7.1.

    Refuses: a missing required section or key; an unknown section; an unknown key, named as a
    physical value when it names a quantity with a unit and holds a number (principle 2); a
    value of the wrong type or out of range; an anchor rule the manifest cannot express; a path
    that resolves outside the manifest's own directory (SPEC_00 section 8); a file name that
    does not carry the profile's prefix (SPEC_00 section 2.2).
    """
    path = Path(path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: manifest does not exist")
    raw_bytes = path.read_bytes()
    try:
        document = tomllib.loads(raw_bytes.decode("utf-8"))
    except tomllib.TOMLDecodeError as exc:
        raise ControlFileError(f"{path}: not valid TOML: {exc}") from None

    unknown_sections = sorted(set(document) - set(_VOCABULARY))
    if unknown_sections:
        raise ControlFileError(
            f"{path}: unknown section(s) {unknown_sections}. SPEC_00 section 7.1 defines "
            f"{sorted(_VOCABULARY)}."
        )

    for section, (section_required, keys) in _VOCABULARY.items():
        if section not in document:
            if section_required:
                raise ControlFileError(f"{path}: required section [{section}] is missing.")
            continue
        table = document[section]
        for key, value in table.items():
            if key in keys:
                continue
            if _looks_physical(key, value):
                raise ControlFileError(
                    f"{path}: [{section}] {key} = {value!r} is a physical value in a control "
                    "file. SPEC_00 principle 2 keeps physical values in data files; the "
                    "manifest holds pointers and declared choices, and the parser refuses any "
                    "key it does not know."
                )
            raise ControlFileError(
                f"{path}: [{section}] carries the unknown key {key!r}. SPEC_00 section 7 makes "
                f"an unknown key an error. [{section}] accepts {sorted(keys)}."
            )
        for key, (kind, key_required) in keys.items():
            if key not in table:
                if key_required:
                    raise ControlFileError(f"{path}: [{section}] is missing required key {key!r}.")
                continue
            if not _is_type(table[key], kind):
                raise ControlFileError(
                    f"{path}: [{section}] {key} = {table[key]!r} is not a {kind} value."
                )

    slug = document["profile"]["slug"]
    anchor = document["anchor_isobar"]
    latitude = document["anchor_latitude"]
    geoid = document["geoid"]
    _check_positive(path, "anchor_isobar", "pressure_Pa", anchor["pressure_Pa"])
    _check_positive(path, "anchor_latitude", "fixed_point_tolerance_deg",
                    latitude["fixed_point_tolerance_deg"])
    _check_positive(path, "anchor_latitude", "max_iterations", latitude["max_iterations"],
                    integer=True)
    _check_positive(path, "geoid", "anchor_surface_Pa", geoid["anchor_surface_Pa"])
    _check_positive(path, "geoid", "convergence_m", geoid["convergence_m"])
    if geoid["anchor_rule"] not in _MANIFEST_ANCHOR_RULES:
        raise ControlFileError(
            f"{path}: [geoid] anchor_rule = {geoid['anchor_rule']!r}; a manifest may name "
            f"{sorted(_MANIFEST_ANCHOR_RULES)}. The 'latitude' rule of lib.geoid needs an anchor "
            "latitude, and SPEC_00 section 7.1 has no key to carry one."
        )
    diagnostics, sensitivity = _check_optional_sections(path, document)

    directory = path.parent

    def inside(name: str, what: str) -> Path:
        resolved = (directory / name).resolve()
        if resolved.parent != directory:
            raise ControlFileError(
                f"{path}: {what} = {name!r} resolves to {resolved}, outside the manifest's own "
                "directory. SPEC_00 section 8 refuses a manifest that points outside it."
            )
        if not resolved.name.startswith(f"{slug}_"):
            raise ControlFileError(
                f"{path}: {what} = {name!r} does not carry the profile prefix '{slug}_' "
                "(SPEC_00 section 2.2)."
            )
        return resolved

    inputs = MappingProxyType({key: inside(document["inputs"][key], f"[inputs] {key}")
                               for key in INPUT_KINDS})
    product = inside(document["output"]["product"], "[output] product")

    return ReductionManifest(
        path=path,
        sha256=cio.sha256(path),
        text=raw_bytes.decode("utf-8"),
        slug=slug,
        description=document["profile"]["description"],
        inputs=inputs,
        anchor_isobar_Pa=float(anchor["pressure_Pa"]),
        fixed_point_tolerance_deg=float(latitude["fixed_point_tolerance_deg"]),
        max_iterations=int(latitude["max_iterations"]),
        anchor_surface_Pa=float(geoid["anchor_surface_Pa"]),
        anchor_quantity=geoid["anchor_quantity"],
        anchor_rule=geoid["anchor_rule"],
        convergence_m=float(geoid["convergence_m"]),
        product=product,
        diagnostics=MappingProxyType(diagnostics),
        sensitivity=MappingProxyType(sensitivity),
    )


# ---------------------------------------------------------------------------
# The six inputs, SPEC_02 Step 1
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ReductionInputs:
    """The six inputs, loaded into memory, with the SHA-256 and commit each carried on disk."""

    manifest: ReductionManifest
    thermo: object
    composition: object
    geodesy: object
    gravity: object
    rotation: object
    wind: object
    sha256: MappingProxyType
    commits: MappingProxyType


def _load_into_memory(path: Path, kind: str):
    """Read under the kind, pull the values into memory, and release the file."""
    handle = cio.read(path, kind)
    try:
        loaded = handle.load()
    finally:
        handle.close()
    return loaded


def load_reduction_inputs(manifest: ReductionManifest) -> ReductionInputs:
    """Read the six inputs under their kinds and check them against each other.

    Every check runs before any arithmetic, so a reduction either starts on inputs that agree
    or does not start. Refuses when (SPEC_00 section 8; SPEC_02 Step 1):

    * any input carries a `casspian_git_commit` ending in `-dirty`;
    * the wind file's rotation system name or rate differs from the rotation file's;
    * the wind file's `reference_level_pressure_Pa` is not a node of its pressure grid;
    * the manifest's anchor isobar is not a surface of the geodesy file, or not a tabulated
      level of the thermo file (the reduction needs `h_ref` at that exact level, not an
      interpolation); both facts are reported together;
    * the geoid anchor surface is not a geodesy surface, or the geodesy file lacks the named
      anchor quantity;
    * the composition file's levels are not the thermo file's levels, bit for bit.
    """
    loaded, hashes, commits = {}, {}, {}
    for key, kind in INPUT_KINDS.items():
        path = manifest.inputs[key]
        if not path.exists():
            raise ControlFileError(f"{manifest.path}: [inputs] {key} names {path.name}, which "
                                   "does not exist.")
        dataset = _load_into_memory(path, kind)
        commit = str(dataset.attrs.get("casspian_git_commit", ""))
        if commit.endswith("-dirty"):
            raise ControlFileError(
                f"{path.name} carries casspian_git_commit = {commit!r}. SPEC_00 section 8: a "
                "file offered as an input to refrac may not carry -dirty. Rebuild it from a "
                "committed tree."
            )
        loaded[key], hashes[key], commits[key] = dataset, cio.sha256(path), commit

    wind, rotation = loaded["wind"], loaded["rotation"]
    wind_name = str(wind.attrs["rotation_system_name"])
    rotation_name = str(rotation.attrs["system_name"])
    wind_rate = float(wind.attrs["rotation_rate_rad_s"])
    rotation_rate = float(rotation["angular_rate_rad_s"])
    if wind_name != rotation_name or wind_rate != rotation_rate:
        raise ControlFileError(
            f"the wind file is in {wind_name!r} at {wind_rate!r} rad/s but the rotation file is "
            f"{rotation_name!r} at {rotation_rate!r} rad/s. SPEC_00 section 6.6 requires them "
            "to be the same system."
        )

    reference = float(wind["reference_level_pressure_Pa"])
    if not np.any(np.asarray(wind["pressure_Pa"].values) == reference):
        raise ControlFileError(
            f"the wind file's reference_level_pressure_Pa = {reference!r} is not a node of its "
            "pressure grid."
        )

    geodesy, thermo = loaded["geodesy"], loaded["thermo"]
    surfaces = np.asarray(geodesy["surface_pressure_Pa"].values, dtype="float64")
    levels = np.asarray(thermo["pressure_Pa"].values, dtype="float64")
    anchor = manifest.anchor_isobar_Pa
    is_surface = bool(np.any(surfaces == anchor))
    is_level = bool(np.any(levels == anchor))
    if not (is_surface and is_level):
        failures = []
        if not is_surface:
            failures.append(f"not a surface of the geodesy file (surfaces {surfaces.tolist()} Pa)")
        if not is_level:
            failures.append("not a tabulated level of the thermo file, and the reduction needs "
                            "h_ref at that exact level rather than an interpolation")
        raise ControlFileError(
            f"{manifest.path.name}: [anchor_isobar] pressure_Pa = {anchor:g} Pa is "
            + " and ".join(failures) + "."
        )

    if not bool(np.any(surfaces == manifest.anchor_surface_Pa)):
        raise ControlFileError(
            f"{manifest.path.name}: [geoid] anchor_surface_Pa = {manifest.anchor_surface_Pa:g} "
            f"Pa is not a surface of the geodesy file (surfaces {surfaces.tolist()} Pa)."
        )
    if manifest.anchor_quantity not in geodesy.variables:
        raise ControlFileError(
            f"{manifest.path.name}: [geoid] anchor_quantity = {manifest.anchor_quantity!r} is "
            "not a variable of the geodesy file."
        )

    composition_levels = np.asarray(
        loaded["composition"].dataset["pressure_Pa"].values, dtype="float64"
    )
    if not np.array_equal(composition_levels, levels):
        raise ControlFileError(
            "the composition file's levels are not the thermo file's levels, bit for bit. The "
            "reduction forms refractivity level by level and does not interpolate between them."
        )

    return ReductionInputs(
        manifest=manifest,
        sha256=MappingProxyType(hashes),
        commits=MappingProxyType(commits),
        **loaded,
    )
