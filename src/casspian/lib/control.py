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
import warnings
from dataclasses import dataclass
from datetime import date as _date
from pathlib import Path
from types import MappingProxyType

import numpy as np
import xarray as xr

from casspian.lib import io as cio
from casspian.lib.schema import UNIT_SUFFIXES, check_wind_components, check_wind_poles

__all__ = [
    "ControlFileError",
    "load_section",
    "reject_physical_values",
    "build_role",
    "build_epoch",
    "BUILD_ROLES",
    "ReductionManifest",
    "ReductionInputs",
    "read_reduction_manifest",
    "load_reduction_inputs",
    "INPUT_KINDS",
    "RUN_INPUT_KINDS",
    "CLOSURE_DROPPED_ATTRIBUTES",
    "RunAnchor",
    "RunNamelist",
    "RunInputs",
    "LoadedAnchor",
    "ClosureComparison",
    "read_run_namelist",
    "load_run_inputs",
    "check_closure_inputs",
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


#: SPEC_01 v0.26: the values of the `role` key every build-file section carries.
BUILD_ROLES = ("reduction", "forward")


def build_role(control: dict, path, section: str) -> str:
    """The required `role` of a build-file section, checked (SPEC_01 v0.26).

    The key is required by every tool's section table, so a section without it is refused by
    `load_section` before this is reached; there is no default in code. The value becomes the
    product's `role` global.
    """
    role = control.get("role")
    if role not in BUILD_ROLES:
        raise ControlFileError(
            f"{path}: [{section}] role = {role!r}; a build control file section declares "
            f"role as one of {list(BUILD_ROLES)} (SPEC_01 v0.26)."
        )
    return role


def build_epoch(control: dict, path, section: str) -> str:
    """The `epoch` of a build-file section that dates a file with no date of its own.

    SPEC_01 v0.28, SPEC_00 v0.19 section 5: a harmonic set or a rotation system carries the date of
    the observation it serves, declared in its build-file section. Refused unless it is an ISO 8601
    date.
    """
    epoch = control.get("epoch")
    try:
        _date.fromisoformat(str(epoch))
    except ValueError:
        raise ControlFileError(
            f"{path}: [{section}] epoch = {epoch!r} is not an ISO 8601 date (SPEC_00 v0.19 "
            "section 5, SPEC_01 v0.28)."
        ) from None
    return str(epoch)


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

#: SPEC_00 section 7.1 (v0.12). Section name -> (section required, {key: (type, required)}).
#: Inside the optional section the keys a switched-on feature needs are enforced in
#: `_check_optional_sections`, not here, because their requirement depends on the switch.
#: v0.12 removed `[sensitivity]`: no sensitivity study lives in the pipeline, so a manifest
#: carrying one is refused as an unknown section.
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
}

#: The anchor rules a manifest may name. `lib.geoid.wind_geoid` also knows "latitude", but that
#: rule needs an anchor latitude and the section 7.1 vocabulary has no key to carry one.
_MANIFEST_ANCHOR_RULES = frozenset({"mean_polar_radius", "north_pole", "south_pole",
                                    "equatorial_radius"})

#: SPEC_00 section 7.1 v0.15: the kind D quantity each manifest anchor rule anchors on. A rule
#: and a quantity that do not belong together are refused, here and by the build tool.
ANCHOR_QUANTITY_FOR_RULE = MappingProxyType({
    "mean_polar_radius": "radius_polar_m",
    "north_pole": "radius_polar_m",
    "south_pole": "radius_polar_m",
    "equatorial_radius": "radius_equatorial_m",
})
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


def _check_optional_sections(path, document) -> dict:
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
    return diagnostics


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
    expected_quantity = ANCHOR_QUANTITY_FOR_RULE[geoid["anchor_rule"]]
    if geoid["anchor_quantity"] != expected_quantity:
        raise ControlFileError(
            f"{path}: [geoid] anchor_rule = {geoid['anchor_rule']!r} anchors on "
            f"{expected_quantity!r}, but anchor_quantity = {geoid['anchor_quantity']!r}. "
            "SPEC_00 section 7.1 v0.15 refuses a rule and a quantity that do not belong together."
        )
    diagnostics = _check_optional_sections(path, document)

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
        commit = _refuse_dirty_commit(path, dataset.attrs, "refrac")
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


# ---------------------------------------------------------------------------
# The run namelist, SPEC_00 section 7.2, the subset of SPEC_03 Step 3
# ---------------------------------------------------------------------------

#: The four inputs a forward run reads, by namelist key and kind. No geodesy (SPEC_00 section 8).
RUN_INPUT_KINDS = MappingProxyType({
    "composition": "composition",
    "gravity": "gravity",
    "rotation": "rotation",
    "wind": "wind",
})

#: The modes of SPEC_00 section 7.2 and the ones implemented: closure by SPEC_03 Step 3,
#: transfer by SPEC_04 Step 0.
_RUN_MODES = ("closure", "transfer")
_RUN_MODES_IMPLEMENTED = ("closure", "transfer")
_P_B_RULES = ("anchor_profile_top",)
_P_B_LOCATIONS = ("top_of_anchor_profile",)

#: SPEC_00 section 7.2 sections and keys that belong to a later specification in closure mode,
#: refused as not implemented rather than as unknown, so that a namelist written for one mode
#: fails in the other for the right reason. Transfer mode requires all of them.
_NOT_IMPLEMENTED_SECTIONS = ("grid", "numerics", "estimation")
_NOT_IMPLEMENTED_KEYS = {"isobars": ("datum_isobar_Pa",)}

#: SPEC_04 Step 0 deliverable 4: the `[numerics]` tables and the schemes implemented for each.
#: The geopotential and the hydrostatic integration are SPEC_03's closed rules and are not
#: namelist keys; a namelist that names either is refused as not implemented, as in closure
#: mode, so that a namelist written against SPEC_00 section 7.2 fails for the right reason.
_NUMERICS_SCHEMES = {
    "reference_surface": ("rk4",),
    "shear_integral": ("trapezoid",),
    "isobar_tracing": ("rk4",),
    "transfer": ("trapezoid",),
    "altitude": ("trapezoid",),
}
_NUMERICS_CLOSED_IN_SPEC_03 = ("geopotential", "hydrostatic")
_OUTER_LOOP_KEYS = {
    "relative_tolerance_ln_p": (_NUMBER, True),
    "max_iterations": (_INTEGER, True),
}

#: SPEC_04 Step 0 deliverable 4: the sections transfer mode adds, and their keys.
_TRANSFER_VOCABULARY = {
    "target": (True, {"latitude_planetocentric_deg": (_NUMBER, True)}),
    # SPEC_04 v0.9: the mesh's extent is not declared. Step 2 builds it from the anchors'
    # levels and grows it when a traced curve needs more, so the only numbers here are the two
    # spacings (decision N: a number the code can compute is never asked of the user).
    "grid": (True, {
        "geopotential_spacing_m2s2": (_NUMBER, True),
        "latitude_spacing_deg": (_NUMBER, True),
    }),
    "estimation": (True, {
        "gauge_latitude_rule": (_TEXT, True),
        "kernel_uncertainty_per_rad": (_NUMBER, True),
        "model_error_correlation_length_deg": (_NUMBER, True),
    }),
}
_GAUGE_LATITUDE_RULES = ("weighted_centroid",)

#: SPEC_04 Step 0 deliverable 4 and SPEC_00 section 7.2: the anchor's role. 1 construction,
#: 0 validation; anything else is refused, since the calibration weights between them are the
#: combination specification's and nothing here would honor them.
_ANCHOR_WEIGHTS = (0.0, 1.0)

_RUN_VOCABULARY = {
    "run": (True, {
        "name": (_TEXT, True),
        "description": (_TEXT, True),
        "mode": (_TEXT, True),
        "solar_longitude_deg": (_NUMBER, True),
        "date": (_TEXT, False),
    }),
    "inputs": (True, {key: (_TEXT, True) for key in RUN_INPUT_KINDS}),
    "hydrostatic_boundary": (True, {
        "p_b_rule": (_TEXT, False),
        "p_b_Pa": (_NUMBER, False),
        "p_b_location": (_TEXT, True),
    }),
    "isobars": (True, {"gauge_isobar_Pa": (_NUMBER, True)}),
    "output": (True, {"directory": (_TEXT, True), "product": (_TEXT, True)}),
    "diagnostics": (False, {
        "figures": (_FLAG, False),
        "format": (_TEXT, False),
        "dpi": (_INTEGER, False),
    }),
}
_ANCHOR_KEYS = {
    "slug": (_TEXT, True),
    "path": (_TEXT, True),
    "weight": (_NUMBER, True),
    "measurement_uncertainty_scale": (_NUMBER, True),
}

#: SPEC_03 v0.10 Step 3 deliverable 2 (section 8 ruling 2): what the closure comparison drops
#: before comparing. The writer globals, the naming and role attributes, and every provenance
#: attribute that carries a path or a hash of another file; any other attribute whose value
#: contains `sha256:` is dropped as well.
CLOSURE_DROPPED_ATTRIBUTES = frozenset({
    "created_by", "created_at", "casspian_git_commit", "casspian_version", "history",
    "title", "profile_or_run", "role", "composition_role",
    "input_hashes", "control_file", "raw_bundle", "raw_sources", "master_table",
    "master_table_hash", "latitude_conversion_inputs", "decomposition_geometry",
})


@dataclass(frozen=True)
class RunAnchor:
    """One `[[anchors]]` entry, its path resolved."""

    slug: str
    path: Path
    weight: float
    measurement_uncertainty_scale: float

    @property
    def is_construction(self) -> bool:
        """Whether this anchor builds the estimate (SPEC_00 section 7.2, SPEC_04 Step 0).

        `weight = 1` is a construction anchor and enters the gauge latitude and the anchor
        constant; `weight = 0` is a validation anchor, which is propagated, placed and traced
        like any other and reports its `C_i` and every `D_ij`, but enters neither.
        """
        return self.weight == 1.0


@dataclass(frozen=True)
class RunNamelist:
    """A parsed, checked run namelist. Paths are absolute; the text is kept verbatim."""

    path: Path
    sha256: str
    text: str
    run_directory: Path
    name: str
    description: str
    mode: str
    solar_longitude_deg: float
    date: object
    anchors: tuple
    inputs: MappingProxyType
    p_b_rule: object
    p_b_Pa: object
    p_b_location: str
    gauge_isobar_Pa: float
    output_directory: Path
    product: Path
    diagnostics: MappingProxyType
    # SPEC_04 Step 0 deliverable 4. All five are None or empty in closure mode, which refuses
    # every one of them as not implemented.
    target_latitude_deg: object = None
    datum_isobar_Pa: object = None
    grid: MappingProxyType = MappingProxyType({})
    numerics: MappingProxyType = MappingProxyType({})
    estimation: MappingProxyType = MappingProxyType({})


def _check_run_table(path, section, table, keys):
    """Keys of one namelist table: not implemented, geodesy, physical, unknown, missing, type."""
    for key, value in table.items():
        if key in keys:
            continue
        if key in _NOT_IMPLEMENTED_KEYS.get(section, ()):
            raise ControlFileError(
                f"{path}: [{section}] {key} is not implemented in this specification (SPEC_03 "
                "Step 3 implements the closure subset of SPEC_00 section 7.2; SPEC_04 adds it)."
            )
        if section == "inputs" and key == "geodesy":
            raise ControlFileError(
                f"{path}: [inputs] names a geodesy file. forward refuses one: the reference "
                "surface comes from the run's gravity, rotation and wind and the anchor's frozen "
                "r0(phi_c) (SPEC_00 section 8, decision 3)."
            )
        if _looks_physical(key, value):
            raise ControlFileError(
                f"{path}: [{section}] {key} = {value!r} is a physical value in a control file. "
                "SPEC_00 principle 2 keeps physical values in data files, and the namelist "
                "parser refuses any key it does not know."
            )
        raise ControlFileError(
            f"{path}: [{section}] carries the unknown key {key!r}. SPEC_00 section 7 makes an "
            f"unknown key an error. [{section}] accepts {sorted(keys)}."
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


def _check_numerics(path, document) -> MappingProxyType:
    """The `[numerics]` tables of transfer mode (SPEC_04 Step 0 deliverable 4).

    One table per integration and one for the outer loop, each required, each with the scheme
    implemented for it. The geopotential and hydrostatic rules are SPEC_03's closed ones and
    are not namelist keys, so a table naming either is refused as not implemented, distinctly
    from the unknown-table refusal.
    """
    numerics = document["numerics"]
    for name in _NUMERICS_CLOSED_IN_SPEC_03:
        if name in numerics:
            raise ControlFileError(
                f"{path}: [numerics.{name}] is not implemented as a namelist key: the "
                f"{name} rule is closed by SPEC_03 and is not chosen per run (SPEC_04 Step 0 "
                "deliverable 4)."
            )
    known = set(_NUMERICS_SCHEMES) | {"outer_loop"}
    unknown = sorted(set(numerics) - known)
    if unknown:
        raise ControlFileError(
            f"{path}: unknown [numerics] table(s) {unknown}. SPEC_04 Step 0 deliverable 4 "
            f"defines {sorted(known)}."
        )
    missing = sorted(known - set(numerics))
    if missing:
        raise ControlFileError(
            f"{path}: [numerics] is missing the required table(s) {missing}; every integration "
            "declares its scheme (SPEC_04 Step 0 deliverable 4)."
        )
    out = {}
    for name, schemes in _NUMERICS_SCHEMES.items():
        _check_run_table(path, f"numerics.{name}", numerics[name], {"scheme": (_TEXT, True)})
        scheme = numerics[name]["scheme"]
        if scheme not in schemes:
            raise ControlFileError(
                f"{path}: [numerics.{name}] scheme = {scheme!r} is not implemented; SPEC_04 "
                f"implements {list(schemes)} for it."
            )
        out[name] = scheme
    _check_run_table(path, "numerics.outer_loop", numerics["outer_loop"], _OUTER_LOOP_KEYS)
    loop = numerics["outer_loop"]
    _check_positive(path, "numerics.outer_loop", "relative_tolerance_ln_p",
                    loop["relative_tolerance_ln_p"])
    _check_positive(path, "numerics.outer_loop", "max_iterations", loop["max_iterations"],
                    integer=True)
    out["outer_loop"] = MappingProxyType({
        "relative_tolerance_ln_p": float(loop["relative_tolerance_ln_p"]),
        "max_iterations": int(loop["max_iterations"]),
    })
    return MappingProxyType(out)


def read_run_namelist(path) -> RunNamelist:
    """Parse and check `<run>.toml` against SPEC_00 section 7.2 in the mode it declares.

    SPEC_03 Step 3 deliverable 2 and SPEC_04 Step 0 deliverable 4. Refuses: invalid TOML; a
    mode that is not one of SPEC_00 section 7.2's or is not implemented; in closure mode
    `[grid]`, `[numerics]`, `[estimation]`, `datum_isobar_Pa` or `[target]`; in transfer mode a
    `[numerics]` table that is unknown, missing, or names a rule SPEC_03 closed, and a scheme
    that is not implemented; an unknown section or key, as a physical value where it names a
    quantity with a unit; a geodesy key; a missing required section or key; other than exactly
    one `[[anchors]]` entry in closure mode; an anchor `weight` that is neither 1 nor 0; not
    exactly one of `p_b_rule` and `p_b_Pa`, `p_b_Pa` in closure mode, or a rule or location
    other than the accepted ones; a gauge latitude rule that is not implemented; a grid spacing
    that is not positive or a range that is not two increasing values; an input path outside
    `inputs/` of the run directory or without the run prefix; an output directory outside the
    run directory or a product without the prefix; a namelist file not named for its run; a
    season outside [0, 360) degrees; a date that is not an ISO 8601 date.
    """
    path = Path(path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: namelist does not exist")
    raw_bytes = path.read_bytes()
    try:
        document = tomllib.loads(raw_bytes.decode("utf-8"))
    except tomllib.TOMLDecodeError as exc:
        raise ControlFileError(f"{path}: not valid TOML: {exc}") from None

    if "run" not in document:
        raise ControlFileError(f"{path}: required section [run] is missing.")
    _check_run_table(path, "run", document["run"], _RUN_VOCABULARY["run"][1])
    run = document["run"]
    mode = run["mode"]
    if mode not in _RUN_MODES:
        raise ControlFileError(
            f"{path}: [run] mode = {mode!r}; SPEC_00 section 7.2 defines {list(_RUN_MODES)}."
        )
    if mode not in _RUN_MODES_IMPLEMENTED:
        raise ControlFileError(
            f"{path}: [run] mode = {mode!r} is not implemented in this specification (SPEC_03 "
            "implements closure; SPEC_04 adds transfer)."
        )
    transfer = mode == "transfer"
    if not transfer:
        for section in document:
            if section in _NOT_IMPLEMENTED_SECTIONS:
                raise ControlFileError(
                    f"{path}: [{section}] is not implemented in this specification (SPEC_03 "
                    "Step 3 implements the closure subset of SPEC_00 section 7.2; SPEC_04 adds "
                    "it)."
                )
        if "target" in document:
            raise ControlFileError(
                f"{path}: [target] is not accepted in closure mode: the target of a closure is "
                "the anchor's own latitude phi_c (SPEC_03 Step 3 deliverable 2)."
            )
    vocabulary = dict(_RUN_VOCABULARY)
    if transfer:
        vocabulary.update(_TRANSFER_VOCABULARY)
        vocabulary["isobars"] = (True, dict(_RUN_VOCABULARY["isobars"][1],
                                            datum_isobar_Pa=(_NUMBER, True)))
        vocabulary["numerics"] = (True, {})       # its tables are checked by _check_numerics
    allowed = set(vocabulary) | {"anchors"}
    unknown = sorted(set(document) - allowed)
    if unknown:
        raise ControlFileError(
            f"{path}: unknown section(s) {unknown}. SPEC_00 section 7.2, as implemented for "
            f"{mode} mode, defines {sorted(allowed)}."
        )
    for section, (section_required, keys) in vocabulary.items():
        # [numerics] holds tables and not keys; `_check_numerics` checks it below.
        if section in ("run", "numerics"):
            continue
        if section not in document:
            if section_required:
                raise ControlFileError(f"{path}: required section [{section}] is missing.")
            continue
        _check_run_table(path, section, document[section], keys)
    if transfer and "numerics" not in document:
        raise ControlFileError(f"{path}: required section [numerics] is missing.")

    anchors_table = document.get("anchors")
    if not isinstance(anchors_table, list) or not anchors_table:
        raise ControlFileError(
            f"{path}: at least one [[anchors]] entry is required (SPEC_00 section 7.2)."
        )
    if mode == "closure" and len(anchors_table) != 1:
        raise ControlFileError(
            f"{path}: closure mode takes exactly one [[anchors]] entry; this namelist has "
            f"{len(anchors_table)} (SPEC_03 Step 3 deliverable 2)."
        )
    for entry in anchors_table:
        if not isinstance(entry, dict):
            raise ControlFileError(f"{path}: [[anchors]] entries must be tables.")
        _check_run_table(path, "anchors", entry, _ANCHOR_KEYS)
        if float(entry["weight"]) not in _ANCHOR_WEIGHTS:
            raise ControlFileError(
                f"{path}: [[anchors]] slug = {entry['slug']!r} carries weight = "
                f"{entry['weight']!r}. SPEC_00 section 7.2 makes the weight the anchor's role, "
                "1 for construction and 0 for validation; the calibration weights between them "
                "are the combination specification's and are not implemented here (SPEC_04 "
                "Step 0 deliverable 4)."
            )
        if float(entry["measurement_uncertainty_scale"]) <= 0.0:
            raise ControlFileError(
                f"{path}: [[anchors]] slug = {entry['slug']!r} carries "
                f"measurement_uncertainty_scale = {entry['measurement_uncertainty_scale']!r}; "
                "it multiplies the anchor's own uncertainty and must be positive."
            )

    name = run["name"]
    if path.stem != name:
        raise ControlFileError(
            f"{path}: the namelist of run {name!r} is {name}.toml (SPEC_00 section 7.2); this "
            f"file is {path.name}."
        )
    season = float(run["solar_longitude_deg"])
    if not 0.0 <= season < 360.0:
        raise ControlFileError(
            f"{path}: [run] solar_longitude_deg = {season!r} is not a season; SPEC_00 v0.17 "
            "section 5 defines it in [0, 360) degrees."
        )
    date = run.get("date")
    if date is not None:
        try:
            _date.fromisoformat(date)
        except ValueError:
            raise ControlFileError(
                f"{path}: [run] date = {date!r} is not an ISO 8601 date."
            ) from None

    boundary = document["hydrostatic_boundary"]
    has_rule, has_value = "p_b_rule" in boundary, "p_b_Pa" in boundary
    if has_rule == has_value:
        raise ControlFileError(
            f"{path}: [hydrostatic_boundary] declares {'both' if has_rule else 'neither of'} "
            "p_b_rule and p_b_Pa; SPEC_00 section 7.2 v0.16 requires exactly one."
        )
    if mode == "closure" and not has_rule:
        raise ControlFileError(
            f"{path}: closure mode takes the boundary pressure by p_b_rule = "
            "'anchor_profile_top', the only form accepted in closure mode (SPEC_00 section 7.2 "
            "v0.16)."
        )
    if has_rule and boundary["p_b_rule"] not in _P_B_RULES:
        raise ControlFileError(
            f"{path}: [hydrostatic_boundary] p_b_rule = {boundary['p_b_rule']!r}; accepted "
            f"{list(_P_B_RULES)}."
        )
    if has_value:
        _check_positive(path, "hydrostatic_boundary", "p_b_Pa", boundary["p_b_Pa"])
    if boundary["p_b_location"] not in _P_B_LOCATIONS:
        raise ControlFileError(
            f"{path}: [hydrostatic_boundary] p_b_location = {boundary['p_b_location']!r}; the "
            f"only value accepted is {_P_B_LOCATIONS[0]!r} (SPEC_00 section 7.2)."
        )
    gauge = document["isobars"]["gauge_isobar_Pa"]
    _check_positive(path, "isobars", "gauge_isobar_Pa", gauge)

    target_latitude = datum = None
    grid = numerics = estimation = MappingProxyType({})
    if transfer:
        target_latitude = float(document["target"]["latitude_planetocentric_deg"])
        if not -90.0 <= target_latitude <= 90.0:
            raise ControlFileError(
                f"{path}: [target] latitude_planetocentric_deg = {target_latitude!r} is not a "
                "latitude."
            )
        datum = float(document["isobars"]["datum_isobar_Pa"])
        _check_positive(path, "isobars", "datum_isobar_Pa", datum)

        table = document["grid"]
        _check_positive(path, "grid", "geopotential_spacing_m2s2",
                        table["geopotential_spacing_m2s2"])
        _check_positive(path, "grid", "latitude_spacing_deg", table["latitude_spacing_deg"])
        grid = MappingProxyType({
            "geopotential_spacing_m2s2": float(table["geopotential_spacing_m2s2"]),
            "latitude_spacing_deg": float(table["latitude_spacing_deg"]),
        })

        numerics = _check_numerics(path, document)

        table = document["estimation"]
        if table["gauge_latitude_rule"] not in _GAUGE_LATITUDE_RULES:
            raise ControlFileError(
                f"{path}: [estimation] gauge_latitude_rule = "
                f"{table['gauge_latitude_rule']!r} is not implemented; SPEC_04 implements "
                f"{list(_GAUGE_LATITUDE_RULES)}."
            )
        _check_positive(path, "estimation", "kernel_uncertainty_per_rad",
                        table["kernel_uncertainty_per_rad"])
        if float(table["model_error_correlation_length_deg"]) != 0.0:
            raise ControlFileError(
                f"{path}: [estimation] model_error_correlation_length_deg = "
                f"{table['model_error_correlation_length_deg']!r}; only 0.0 is implemented, "
                "which recovers Eq. A28. A nonzero length is the kernel correction of Eq. A35 "
                "and belongs to the combination specification (SPEC_04 decision D)."
            )
        estimation = MappingProxyType({
            "gauge_latitude_rule": str(table["gauge_latitude_rule"]),
            "kernel_uncertainty_per_rad": float(table["kernel_uncertainty_per_rad"]),
            "model_error_correlation_length_deg": 0.0,
        })

    diagnostics = _check_optional_sections(path, document)

    run_directory = path.parent
    inputs_directory = run_directory / "inputs"

    def run_input(key: str) -> Path:
        value = document["inputs"][key]
        resolved = (run_directory / value).resolve()
        if resolved.parent != inputs_directory:
            raise ControlFileError(
                f"{path}: [inputs] {key} = {value!r} resolves to {resolved}, not a file in the "
                "run's inputs/ directory. SPEC_00 section 8 refuses a namelist pointing outside "
                "its own directory to anything but a kind N file."
            )
        if not resolved.name.startswith(f"{name}_"):
            raise ControlFileError(
                f"{path}: [inputs] {key} = {value!r} does not carry the run prefix '{name}_' "
                "(SPEC_00 section 8)."
            )
        return resolved

    inputs = MappingProxyType({key: run_input(key) for key in RUN_INPUT_KINDS})
    anchors = tuple(
        RunAnchor(slug=e["slug"], path=(run_directory / e["path"]).resolve(),
                  weight=float(e["weight"]),
                  measurement_uncertainty_scale=float(e["measurement_uncertainty_scale"]))
        for e in anchors_table
    )
    output_directory = (run_directory / document["output"]["directory"]).resolve()
    if output_directory.parent != run_directory:
        raise ControlFileError(
            f"{path}: [output] directory = {document['output']['directory']!r} is not a "
            "directory of the run directory; forward writes nothing outside it (SPEC_00 "
            "section 2.3)."
        )
    product_name = document["output"]["product"]
    product = (output_directory / product_name).resolve()
    if product.parent != output_directory or not product.name.startswith(f"{name}_"):
        raise ControlFileError(
            f"{path}: [output] product = {product_name!r} must be a file name in the output "
            f"directory carrying the run prefix '{name}_' (SPEC_00 section 8)."
        )

    return RunNamelist(
        path=path,
        sha256=cio.sha256(path),
        text=raw_bytes.decode("utf-8"),
        run_directory=run_directory,
        name=name,
        description=run["description"],
        mode=mode,
        solar_longitude_deg=season,
        date=date,
        anchors=anchors,
        inputs=inputs,
        p_b_rule=boundary.get("p_b_rule"),
        p_b_Pa=float(boundary["p_b_Pa"]) if has_value else None,
        p_b_location=boundary["p_b_location"],
        gauge_isobar_Pa=float(gauge),
        output_directory=output_directory,
        product=product,
        diagnostics=MappingProxyType(diagnostics),
        target_latitude_deg=target_latitude,
        datum_isobar_Pa=datum,
        grid=grid,
        numerics=numerics,
        estimation=estimation,
    )


# ---------------------------------------------------------------------------
# The run's inputs and the closure comparison, SPEC_03 Step 3 deliverable 2
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class ClosureComparison:
    """One run input against the anchor's embedded copy of the same kind."""

    kind: str
    identical: bool
    differences: tuple
    dropped: MappingProxyType


@dataclass(frozen=True)
class LoadedAnchor:
    """One anchor as it arrives at the model. SPEC_04 Step 0 deliverable 6.

    Everything on the anchor's own levels, top down as kind N stores them. Every later step
    reads the anchor's uncertainty from here and never from the file.

    The two uncertainty columns are the measurement term and the season term of decision J.
    `sigma_ln_N_measurement` is the uncertainty of `ln N`, that is the file's
    `refractivity_uncertainty` divided by its `refractivity` and multiplied by the entry's
    `measurement_uncertainty_scale`: the weights of Eq. A30 act on `ln N`, so the column that
    carries them is an uncertainty of `ln N` (REPORT_04_step0 finding 1). `sigma_ln_N_season`
    is zero here and `season_term` says so; the propagator specification fills it.

    `label_pressure_Pa` is the label of each level's isobar. It arrives as the anchor's own
    tabulated pressure and is replaced by the anchor's produced pressure once the anchor has
    been produced (SPEC_04 Step 1 deliverable 3).
    """

    slug: str
    path: Path
    tree: object
    weight: float
    measurement_uncertainty_scale: float
    latitude_planetocentric_deg: float
    radius_m: object
    ln_N: object
    label_pressure_Pa: object
    sigma_ln_N_measurement: object
    sigma_ln_N_season: object
    season_term: str
    anchor_season_deg: object
    run_season_deg: float
    season_matches_run: bool
    gauge_level_index: int
    propagation: MappingProxyType = MappingProxyType({})

    @property
    def is_construction(self) -> bool:
        """Whether this anchor enters the gauge latitude and the anchor constant."""
        return self.weight == 1.0


@dataclass(frozen=True)
class RunInputs:
    """The anchors and the four run inputs in memory, with hashes, commits and check results."""

    namelist: RunNamelist
    anchor: object
    composition: object
    gravity: object
    rotation: object
    wind: object
    sha256: MappingProxyType
    commits: MappingProxyType
    gauge_level_index: int
    closure: tuple
    #: SPEC_04 Step 0 deliverable 6: one `LoadedAnchor` per `[[anchors]]` entry, in the
    #: namelist's order. Empty in closure mode, which reads the anchor through `anchor`.
    anchors: tuple = ()
    #: SPEC_04 Step 0 deliverable 5, decision N: what each file declares as its season, beside
    #: the run's under the key `run`. A difference is recorded here and warned, never refused.
    seasons: MappingProxyType = MappingProxyType({})


def _nodes(obj) -> dict:
    """`{group path relative to obj: Dataset}` for a Dataset, a DataTree, or a subtree."""
    if isinstance(obj, xr.DataTree):
        base = obj.path.rstrip("/")
        out = {}
        for node in obj.subtree:
            relative = node.path[len(base):] if base else node.path
            out[relative or "/"] = node.to_dataset(inherit=False)
        return out
    return {"/": obj}


def _dropped_names(attrs) -> list:
    return sorted(k for k, v in attrs.items()
                  if k in CLOSURE_DROPPED_ATTRIBUTES or "sha256:" in str(v))


def _strip(dataset):
    """A copy without the dropped attributes, on the group and its variables, and the names."""
    dropped = _dropped_names(dataset.attrs)
    out = dataset.copy()
    out.attrs = {k: v for k, v in dataset.attrs.items() if k not in dropped}
    for name in list(out.variables):
        removed = _dropped_names(out[name].attrs)
        if removed:
            dropped.extend(f"{name}.{k}" for k in removed)
            out[name].attrs = {k: v for k, v in out[name].attrs.items() if k not in removed}
    return out, tuple(dropped)


def _same_value(a, b) -> bool:
    a_arr, b_arr = np.asarray(a), np.asarray(b)
    if a_arr.shape != b_arr.shape:
        return False
    if a_arr.dtype.kind in "fc" and b_arr.dtype.kind in "fc":
        return bool(np.array_equal(a_arr, b_arr, equal_nan=True))
    return bool(np.all(a_arr == b_arr))


def _differences(run, anchor, where: str) -> list:
    out = []
    run_vars, anchor_vars = set(run.variables), set(anchor.variables)
    if run_vars != anchor_vars:
        out.append(f"{where}: variables only in the run's file {sorted(run_vars - anchor_vars)}, "
                   f"only in the anchor's copy {sorted(anchor_vars - run_vars)}")
    for name in sorted(run_vars & anchor_vars):
        x, y = run[name], anchor[name]
        if x.dims != y.dims or not _same_value(x.values, y.values):
            out.append(f"{where}: variable {name} differs in its values")
            continue
        keys = set(x.attrs) | set(y.attrs)
        bad = sorted(k for k in keys if k not in x.attrs or k not in y.attrs
                     or not _same_value(x.attrs[k], y.attrs[k]))
        if bad:
            out.append(f"{where}: variable {name} differs in the attributes {bad}")
    keys = set(run.attrs) | set(anchor.attrs)
    bad = sorted(k for k in keys if k not in run.attrs or k not in anchor.attrs
                 or not _same_value(run.attrs[k], anchor.attrs[k]))
    if bad:
        out.append(f"{where}: the attributes {bad} differ")
    return out


def check_closure_inputs(inputs: dict, anchor_tree) -> tuple:
    """Compare each run input with the anchor's embedded copy of its kind, root and every group.

    SPEC_03 v0.10 Step 3 deliverable 2. Both sides lose `CLOSURE_DROPPED_ATTRIBUTES` and any
    attribute whose value contains `sha256:`, on the group and on its variables; what is left must
    be identical (`xarray.Dataset.identical`). Returns one `ClosureComparison` per kind, carrying
    every difference found and the attributes dropped from the run's file, per group.
    """
    results = []
    for key in RUN_INPUT_KINDS:
        run_nodes = _nodes(inputs[key])
        anchor_nodes = _nodes(anchor_tree[f"inputs/{key}"])
        differences, dropped = [], {}
        if set(run_nodes) != set(anchor_nodes):
            differences.append(
                f"{key}: groups only in the run's file {sorted(set(run_nodes) - set(anchor_nodes))}, "
                f"only in the anchor's copy {sorted(set(anchor_nodes) - set(run_nodes))}")
        for group in sorted(set(run_nodes) & set(anchor_nodes)):
            run_ds, run_dropped = _strip(run_nodes[group])
            anchor_ds, _ = _strip(anchor_nodes[group])
            dropped[group] = run_dropped
            if not run_ds.identical(anchor_ds):
                found = _differences(run_ds, anchor_ds, f"{key} {group}")
                differences.extend(found or [f"{key} {group}: not identical"])
        results.append(ClosureComparison(kind=key, identical=not differences,
                                         differences=tuple(differences),
                                         dropped=MappingProxyType(dropped)))
    return tuple(results)


def _refuse_dirty_commit(path: Path, attrs, consumer: str, what: str | None = None) -> str:
    """Refuse a file whose `casspian_git_commit` ends in `-dirty`, and return the commit.

    SPEC_00 section 8. The one place `refrac` and `forward` make this refusal, so that a step
    that must rebuild an input file on a working tree has one point to relax, inside its
    acceptance script and named in its output (SPEC_03 section 0, the author's ruling of 14
    September 2026). The two message forms are unchanged: `forward` names the role the file
    plays in the run, `refrac` names the file alone.
    """
    commit = str(attrs.get("casspian_git_commit", ""))
    if commit.endswith("-dirty"):
        where = path.name if what is None else f"{path.name} ({what})"
        raise ControlFileError(
            f"{where} carries casspian_git_commit = {commit!r}. SPEC_00 section 8: a file "
            f"offered as an input to {consumer} may not carry -dirty. Rebuild it from a "
            "committed tree."
        )
    return commit


def _refuse_dirty(path: Path, attrs, what: str) -> str:
    return _refuse_dirty_commit(path, attrs, "forward", what)


def load_run_inputs(namelist: RunNamelist) -> RunInputs:
    """Read the anchor and the four inputs of a run under their kinds, and check them.

    SPEC_03 Step 3 deliverable 2 and SPEC_04 Step 0 deliverable 5. Every check runs before any
    arithmetic. Refuses when: a file does not exist or is not of its kind; an anchor or an input
    carries `-dirty`; an anchor's `profile_or_run` is not its `[[anchors]]` slug; an input does
    not carry the run name in `profile_or_run`, or `role = "forward"`, or (kind C)
    `composition_role = "forward"`; the wind's rotation system or rate differs from the run's
    kind R; the wind's parts do not sum to its total or its poles are not zero; the wind's
    coverage does not span the anchors' levels and latitudes. In closure mode, also when: the
    run's season differs from the anchor's; the composition's levels are not the anchor's; the
    gauge isobar is not a tabulated level of the anchor (the nearest level named); an input
    differs from the anchor's embedded copy. In transfer mode, also when: the wind's latitude
    coverage excludes the target; the composition does not cover every anchor and the target;
    the wind or the composition carries neither the run's season nor `uniform`; an anchor does
    not carry `radius_m` as its coordinate; the gauge isobar is not a tabulated level of every
    anchor. A season mismatch between an anchor and the run is recorded, not refused
    (SPEC_04 decision I).
    """
    if namelist.mode == "transfer":
        return _load_transfer_inputs(namelist)
    if namelist.mode != "closure":
        raise ControlFileError(f"{namelist.path}: mode {namelist.mode!r} is not implemented.")
    hashes, commits = {}, {}
    entry = namelist.anchors[0]
    if not entry.path.exists():
        raise ControlFileError(
            f"{namelist.path}: [[anchors]] path names {entry.path}, which does not exist."
        )
    anchor = _load_into_memory(entry.path, "refractivity")
    commits["anchor"] = _refuse_dirty(entry.path, anchor.attrs, "anchor")
    hashes["anchor"] = cio.sha256(entry.path)
    slug = str(anchor.attrs.get("profile_or_run", ""))
    if slug != entry.slug:
        raise ControlFileError(
            f"{namelist.path}: [[anchors]] slug = {entry.slug!r} but {entry.path.name} carries "
            f"profile_or_run = {slug!r} (SPEC_00 section 7.2)."
        )

    loaded = {}
    for key, kind in RUN_INPUT_KINDS.items():
        path = namelist.inputs[key]
        if not path.exists():
            raise ControlFileError(
                f"{namelist.path}: [inputs] {key} names {path.name}, which does not exist."
            )
        dataset = _load_into_memory(path, kind)
        commits[key] = _refuse_dirty(path, dataset.attrs, f"run input {key}")
        hashes[key] = cio.sha256(path)
        attrs = dataset.attrs
        if str(attrs.get("profile_or_run")) != namelist.name:
            raise ControlFileError(
                f"{path.name} carries profile_or_run = {attrs.get('profile_or_run')!r}; a run's "
                f"input carries the run prefix {namelist.name!r} (SPEC_00 section 8)."
            )
        if attrs.get("role") != "forward":
            raise ControlFileError(
                f"{path.name} carries role = {attrs.get('role')!r}; a run's input is role "
                "'forward' (SPEC_03 Step 3 deliverable 2)."
            )
        if kind == "composition" and attrs.get("composition_role") != "forward":
            raise ControlFileError(
                f"{path.name} carries composition_role = {attrs.get('composition_role')!r}; a "
                "run's composition is 'forward' (SPEC_00 section 6.2)."
            )
        loaded[key] = dataset

    wind, rotation = loaded["wind"], loaded["rotation"]
    wind_name = str(wind.attrs["rotation_system_name"])
    rotation_name = str(rotation.attrs["system_name"])
    wind_rate = float(wind.attrs["rotation_rate_rad_s"])
    rotation_rate = float(rotation["angular_rate_rad_s"])
    if wind_name != rotation_name or wind_rate != rotation_rate:
        raise ControlFileError(
            f"the run's wind is in {wind_name!r} at {wind_rate!r} rad/s but its rotation file is "
            f"{rotation_name!r} at {rotation_rate!r} rad/s (SPEC_00 section 6.6)."
        )
    try:
        check_wind_components(wind, namelist.inputs["wind"].name)
        check_wind_poles(wind, namelist.inputs["wind"].name)
    except Exception as exc:
        raise ControlFileError(str(exc)) from None

    root = anchor.to_dataset(inherit=False)
    levels = np.asarray(anchor["inputs/thermo"].to_dataset(inherit=False)["pressure_Pa"].values,
                        dtype="float64")
    phi_c = float(root["latitude_planetocentric_deg"].values)
    coverage_p = np.asarray(wind.attrs["coverage_pressure_Pa"], dtype="float64")
    coverage_lat = np.asarray(wind.attrs["coverage_latitude_planetocentric_deg"], dtype="float64")
    if not (coverage_p.min() <= levels.min() and coverage_p.max() >= levels.max()):
        raise ControlFileError(
            f"the run's wind covers {coverage_p.tolist()} Pa, which does not span the anchor's "
            f"levels {levels.min()!r} to {levels.max()!r} Pa; the model does not extrapolate "
            "the wind (SPEC_00 section 6.6)."
        )
    if not (coverage_lat.min() <= phi_c <= coverage_lat.max()):
        raise ControlFileError(
            f"the run's wind covers latitudes {coverage_lat.tolist()} deg, which do not include "
            f"the anchor's phi_c = {phi_c!r} deg (SPEC_00 section 6.6)."
        )

    anchor_season = anchor.attrs.get("solar_longitude_deg")
    if anchor_season is None or float(anchor_season) != namelist.solar_longitude_deg:
        raise ControlFileError(
            f"{namelist.path.name}: [run] solar_longitude_deg = {namelist.solar_longitude_deg!r} "
            f"but the anchor {entry.path.name} is at {anchor_season!r}; not a closure. In closure "
            "mode the run's season equals the anchor's exactly: a closure at another season is "
            "not a closure (SPEC_00 section 7.2 v0.17)."
        )
    composition_levels = np.asarray(loaded["composition"].dataset["pressure_Pa"].values,
                                    dtype="float64")
    if not np.array_equal(composition_levels, levels):
        raise ControlFileError(
            "the run's composition levels are not the anchor's levels, bit for bit; closure mode "
            "runs on the anchor's tabulated levels (SPEC_03 Step 3 deliverable 2)."
        )
    match = np.flatnonzero(levels == namelist.gauge_isobar_Pa)
    if match.size != 1:
        nearest = int(np.argmin(np.abs(levels - namelist.gauge_isobar_Pa)))
        raise ControlFileError(
            f"{namelist.path.name}: [isobars] gauge_isobar_Pa = {namelist.gauge_isobar_Pa!r} is "
            f"not a tabulated level of the anchor; the nearest is level {nearest} at "
            f"{levels[nearest]!r} Pa. In closure mode the gauge is a level, matched exactly "
            "(SPEC_03 Step 1)."
        )

    closure = check_closure_inputs(loaded, anchor)
    failing = [c for c in closure if not c.identical]
    if failing:
        first = failing[0]
        raise ControlFileError(
            f"the run's {first.kind} input is not content-identical to the anchor's embedded "
            f"copy, so this is not a closure (SPEC_03 Step 3 deliverable 2). First difference: "
            f"{first.differences[0]}. All differences: {'; '.join(first.differences)}."
        )

    return RunInputs(
        namelist=namelist,
        anchor=anchor,
        sha256=MappingProxyType(hashes),
        commits=MappingProxyType(commits),
        gauge_level_index=int(match[0]),
        closure=closure,
        **loaded,
    )


# ---------------------------------------------------------------------------
# Transfer mode, SPEC_04 Step 0 deliverables 5 and 6
# ---------------------------------------------------------------------------


def _season_of(attrs) -> object:
    """A file's declared season in degrees, or `None` when it declares itself uniform."""
    if "solar_longitude_deg" in attrs:
        return float(attrs["solar_longitude_deg"])
    return None


def _record_season(path: Path, attrs, run_season: float, what: str):
    """Record a file's season beside the run's, warning when they differ. Never refuses.

    SPEC_04 decision N, which amends decision I's "else refused": a season that differs from
    the run's is a recorded, unmodeled term until the propagator exists, and the warning is
    there because the record alone might be missed. Returns what the file declares: the solar
    longitude in degrees, or `"uniform"`.
    """
    season = _season_of(attrs)
    if season is None:
        return str(attrs.get("season_absent_meaning", "unstated"))
    if season != run_season:
        warnings.warn(
            f"{path.name}: the run's {what} is at solar_longitude_deg = {season!r} but the run "
            f"declares {run_season!r}. The difference is recorded in the product and is an "
            "unmodeled term until the propagator exists (SPEC_04 decisions I and N).",
            stacklevel=2,
        )
    return season


def _load_anchor_object(entry: RunAnchor, namelist: RunNamelist, tree) -> LoadedAnchor:
    """The anchor as it arrives (SPEC_04 Step 0 deliverable 6), from its kind N file.

    The retrieval instance of kind N, which carries no `radius_m`, does not pass the schema
    today, so there is no check of its own here; the retrieval leg adds what it needs
    (SPEC_04 Step 0 deliverable 5, decision N).
    """
    root = tree.to_dataset(inherit=False)
    levels = np.asarray(tree["inputs/thermo"].to_dataset(inherit=False)["pressure_Pa"].values,
                        dtype="float64")
    match = np.flatnonzero(levels == namelist.gauge_isobar_Pa)
    if match.size != 1:
        nearest = int(np.argmin(np.abs(levels - namelist.gauge_isobar_Pa)))
        raise ControlFileError(
            f"{namelist.path.name}: [isobars] gauge_isobar_Pa = {namelist.gauge_isobar_Pa!r} is "
            f"not a tabulated level of the anchor {entry.slug!r}; the nearest is level "
            f"{nearest} at {levels[nearest]!r} Pa. The gauge is a tabulated level of every "
            "occultation anchor, matched exactly; an anchor whose levels do not include it is "
            "not implemented in this specification (SPEC_04 Step 0 deliverable 4, decision M)."
        )
    N = np.asarray(root["refractivity"].values, dtype="float64")
    sigma_N = np.asarray(root["refractivity_uncertainty"].values, dtype="float64")
    anchor_season = _season_of(tree.attrs)
    return LoadedAnchor(
        slug=entry.slug,
        path=entry.path,
        tree=tree,
        weight=entry.weight,
        measurement_uncertainty_scale=entry.measurement_uncertainty_scale,
        latitude_planetocentric_deg=float(root["latitude_planetocentric_deg"].values),
        radius_m=np.asarray(root["radius_m"].values, dtype="float64"),
        ln_N=np.log(N),
        label_pressure_Pa=levels,
        sigma_ln_N_measurement=entry.measurement_uncertainty_scale * sigma_N / N,
        sigma_ln_N_season=np.zeros_like(N),
        season_term="absent",
        anchor_season_deg=anchor_season,
        run_season_deg=namelist.solar_longitude_deg,
        season_matches_run=(anchor_season is not None
                            and anchor_season == namelist.solar_longitude_deg),
        gauge_level_index=int(match[0]),
    )


def _load_transfer_inputs(namelist: RunNamelist) -> RunInputs:
    """The anchors and the four inputs of a transfer run. SPEC_04 Step 0 deliverable 5."""
    hashes, commits, anchors = {}, {}, []
    for entry in namelist.anchors:
        if not entry.path.exists():
            raise ControlFileError(
                f"{namelist.path}: [[anchors]] path names {entry.path}, which does not exist."
            )
        tree = _load_into_memory(entry.path, "refractivity")
        commits[f"anchor:{entry.slug}"] = _refuse_dirty(entry.path, tree.attrs,
                                                        f"anchor {entry.slug}")
        hashes[f"anchor:{entry.slug}"] = cio.sha256(entry.path)
        slug = str(tree.attrs.get("profile_or_run", ""))
        if slug != entry.slug:
            raise ControlFileError(
                f"{namelist.path}: [[anchors]] slug = {entry.slug!r} but {entry.path.name} "
                f"carries profile_or_run = {slug!r} (SPEC_00 section 7.2)."
            )
        anchors.append(_load_anchor_object(entry, namelist, tree))

    loaded = {}
    for key, kind in RUN_INPUT_KINDS.items():
        path = namelist.inputs[key]
        if not path.exists():
            raise ControlFileError(
                f"{namelist.path}: [inputs] {key} names {path.name}, which does not exist."
            )
        dataset = _load_into_memory(path, kind)
        commits[key] = _refuse_dirty(path, dataset.attrs, f"run input {key}")
        hashes[key] = cio.sha256(path)
        attrs = dataset.attrs
        if str(attrs.get("profile_or_run")) != namelist.name:
            raise ControlFileError(
                f"{path.name} carries profile_or_run = {attrs.get('profile_or_run')!r}; a run's "
                f"input carries the run prefix {namelist.name!r} (SPEC_00 section 8)."
            )
        if attrs.get("role") != "forward":
            raise ControlFileError(
                f"{path.name} carries role = {attrs.get('role')!r}; a run's input is role "
                "'forward' (SPEC_03 Step 3 deliverable 2)."
            )
        if kind == "composition" and attrs.get("composition_role") != "forward":
            raise ControlFileError(
                f"{path.name} carries composition_role = {attrs.get('composition_role')!r}; a "
                "run's composition is 'forward' (SPEC_00 section 6.2)."
            )
        loaded[key] = dataset

    wind, rotation = loaded["wind"], loaded["rotation"]
    wind_name = str(wind.attrs["rotation_system_name"])
    rotation_name = str(rotation.attrs["system_name"])
    wind_rate = float(wind.attrs["rotation_rate_rad_s"])
    rotation_rate = float(rotation["angular_rate_rad_s"])
    if wind_name != rotation_name or wind_rate != rotation_rate:
        raise ControlFileError(
            f"the run's wind is in {wind_name!r} at {wind_rate!r} rad/s but its rotation file is "
            f"{rotation_name!r} at {rotation_rate!r} rad/s (SPEC_00 section 6.6)."
        )
    try:
        check_wind_components(wind, namelist.inputs["wind"].name)
        check_wind_poles(wind, namelist.inputs["wind"].name)
    except Exception as exc:
        raise ControlFileError(str(exc)) from None

    # The seasons of W, C and every anchor are recorded beside the run's and a difference is
    # warned, never refused (SPEC_04 decision N, which amends decision I's "else refused").
    # Coverage is the interpolants' business, not the loader's: `lib.windfield` and
    # `lib.composition.column_at` refuse a point outside their data, which is where a coverage
    # failure belongs and where it names the point that failed.
    seasons = {"run": namelist.solar_longitude_deg}
    for key in ("wind", "composition"):
        seasons[key] = _record_season(namelist.inputs[key], loaded[key].attrs,
                                      namelist.solar_longitude_deg, key)
    for anchor in anchors:
        seasons[f"anchor:{anchor.slug}"] = _record_season(
            anchor.path, anchor.tree.attrs, namelist.solar_longitude_deg,
            f"anchor {anchor.slug}")

    return RunInputs(
        namelist=namelist,
        anchor=anchors[0].tree,
        sha256=MappingProxyType(hashes),
        commits=MappingProxyType(commits),
        gauge_level_index=anchors[0].gauge_level_index,
        closure=(),
        anchors=tuple(anchors),
        seasons=MappingProxyType(seasons),
        **loaded,
    )
