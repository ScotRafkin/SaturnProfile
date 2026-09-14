"""The file conventions of SPEC_00 section 5 and the file kinds of section 6, as code.

This module is the single validation layer called for in SPEC_00 section 8. It declares what
every kind requires and it refuses anything that does not conform. It performs no file access;
`lib.io` opens and closes files and calls `validate` here.

Manuscript equations implemented: none. This module encodes file structure only.

Two rules from SPEC_00 shape everything below.

Section 5, naming: a variable or attribute named `latitude`, `r`, `z`, or `p`, with no
convention or unit in its name, is a defect and is refused for every kind, including `raw`.

Section 5, dimensions: absent, assumed, and missing are three different things. A kind that
declares a dimension may omit it only with a global attribute `<dimension>_absent_meaning` of
`uniform` (the quantity does not vary along it, the reader broadcasts) or `point` (the quantity
exists only at the single stated coordinate, and the reader refuses to broadcast).
"""

from __future__ import annotations

from dataclasses import dataclass, field

# ---------------------------------------------------------------------------
# Fixed strings and vocabularies, SPEC_00 section 5
# ---------------------------------------------------------------------------

CONVENTIONS = "CF-1.10, CASSPIAN-0.1"

#: SPEC_00 section 6.4. `harmonic_convention` is a short code, not a description: a code is
#: compared exactly and its meaning lives in one place. CASSPIAN-J1 means even zonal harmonics
#: J_l of degree l, unscaled (the printed value times 1e-6), in the potential
#: V = -(GM/r) [1 - sum_l J_l (R/r)^l P_l(sin phi_c)], with R the normalization radius stored
#: in the file, phi_c planetocentric, and g_N = dV/dr positive inward (handoff section 9A.5).
#: A different convention (tesseral terms, a different sign, normalized coefficients) gets a
#: new code, never a reinterpretation of this one.
HARMONIC_CONVENTION = "CASSPIAN-J1"

#: The codes this reader knows. A kind G file carrying any other code is refused.
KNOWN_HARMONIC_CONVENTIONS = frozenset({HARMONIC_CONVENTION})

#: SPEC_00 section 5. Two of these are easily confused, so section 5 defines them: `derived`
#: is obtained from measurements by a fit or a computation with no additional physical
#: assumption (a harmonic coefficient fitted to tracking data, a pressure from a hydrostatic
#: integration); `inferred` is deduced through a model assumption the source declares (an
#: ammonia abundance from an opacity with an assumed line shape).
PROVENANCE_VALUES = frozenset(
    {
        "measured",
        "derived",
        "assumed",
        "inferred",
        "interpolated",
        "parameterized",
        "extrapolated",
        # Produced by `casspian.forward`. Appears only in files the model writes, so that a
        # delivered field is never mistaken for a retrieval when the two lie side by side.
        "modeled",
        # A coordinate that is a label rather than a measurement: `degree`, `species_name`,
        # `surface`.
        "index",
    }
)

#: SPEC_00 section 5. How an uncertainty was obtained (a sample standard deviation of binned
#: points, a formal error, a quoted range) goes in a free-text `uncertainty_method` attribute,
#: which this registry allows on any uncertainty variable and never validates. A bin standard
#: deviation is a 1sigma estimate and is spelled that way here.
UNCERTAINTY_KINDS = frozenset({"1sigma", "range", "stated"})

#: SPEC_00 section 6.1. Which instance of kind T a file is.
THERMO_INSTANCES = frozenset({"source_profile", "retrieval", "model_output"})

#: SPEC_00 section 6.1. Required when `thermo_instance` is `source_profile`, and forbidden
#: for the other two instances, whose latitude is a coordinate rather than a stated scalar.
SOURCE_PROFILE_GLOBALS = (
    "spacecraft",
    "event",
    "observation_date",
    "frequency_bands",
    "latitude_definition",
    "longitude_deg",
    "longitude_system",
    "height_datum",
    "source_top_boundary",
    "source_gravity_citation",
    "source_rotation_system",
    "source_wind_citation",
    "raw_bundle",
)

ABSENT_MEANINGS = frozenset({"uniform", "point"})

ROLES = frozenset({"reduction", "forward"})

VERTICAL_COORDINATES = frozenset({"pressure_Pa", "geopotential_m2s2"})

POSITIVE_VALUES = frozenset({"up", "down"})

DIRECTION_VALUES = frozenset({"increasing", "decreasing"})

#: Coordinate names that carry a vertical sense, and so must declare `positive`. A harmonic
#: `degree`, a `species` index, a `surface` index and a latitude carry neither `positive` nor
#: `direction` (SPEC_00 section 5).
VERTICAL_COORDINATE_HINTS = ("pressure", "geopotential", "height", "radius")

#: Unit suffixes a variable name may carry. Used only to place `_uncertainty`; longest match
#: wins, so `GM_m3s2` resolves on `_m3s2` and not on `_m`.
UNIT_SUFFIXES = ("_m3s2", "_m2s2", "_rad_s", "_kg_mol", "_deg", "_m3", "_ms", "_Pa", "_K",
                 "_m", "_s")


def uncertainty_companion(name: str) -> str:
    """Return the name of the uncertainty companion of `name`.

    SPEC_00 section 5: the companion is named by inserting `_uncertainty` before the unit
    suffix of the variable it accompanies, so `pressure_Pa` has `pressure_uncertainty_Pa`,
    `GM_m3s2` has `GM_uncertainty_m3s2` and `u_total_ms` has `u_total_uncertainty_ms`. A
    variable with no unit suffix appends it, so `refractivity` has `refractivity_uncertainty`
    and `x_H2` has `x_H2_uncertainty`.
    """
    for suffix in sorted(UNIT_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return f"{name[: -len(suffix)]}_uncertainty{suffix}"
    return f"{name}_uncertainty"

#: A bare name carries no convention and no unit, so it cannot be read unambiguously.
FORBIDDEN_BARE_NAMES = frozenset({"latitude", "r", "z", "p"})

#: Globals `write` fills in; an author who sets them has them overwritten.
WRITER_FILLED_GLOBALS = (
    "Conventions",
    "casspian_kind",
    "casspian_schema_version",
    "created_by",
    "created_at",
    "casspian_git_commit",
    "codata_release",
)

#: Globals the caller must supply for every kind. `history` is created by `write` and
#: appended to thereafter, so it is not demanded of the caller.
AUTHOR_REQUIRED_GLOBALS = ("title", "profile_or_run", "role", "source")


class CasspianSchemaError(Exception):
    """A file, or a dataset about to become one, does not conform to its kind.

    SPEC_00 section 8 calls this refusal. The message always names the item at fault.
    """


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class VarSpec:
    """One variable a kind requires or allows."""

    name: str
    required: bool = True
    #: `None` means the dimensions are not checked (scalars, or a shape the kind leaves open).
    dims: tuple[str, ...] | None = None
    units: str | None = None
    #: Flag variables are int8 with `flag_values` and `flag_meanings`, CF style.
    is_flag: bool = False
    #: An uncertainty companion must exist even when the source states none, holding NaN.
    needs_uncertainty: bool = False
    #: An override, for a kind table that names a companion against the section 5 rule. The
    #: default derives the name from the rule and needs no override for any current kind.
    uncertainty_name: str | None = None

    def companion(self) -> str:
        """The name of this variable's uncertainty companion (SPEC_00 section 5)."""
        return self.uncertainty_name or uncertainty_companion(self.name)


@dataclass(frozen=True)
class KindSpec:
    """One file kind: what it requires, and what version of that requirement this is."""

    name: str
    schema_version: int
    #: Dimensions the kind defines. A declared dimension may be omitted only with an
    #: `<dimension>_absent_meaning` global (SPEC_00 section 5).
    dimensions: tuple[str, ...] = ()
    #: Dimensions that must be present and may not be omitted.
    required_dimensions: tuple[str, ...] = ()
    variables: tuple[VarSpec, ...] = ()
    #: Global attributes this kind requires beyond the section 5 set.
    globals_required: tuple[str, ...] = ()
    #: Groups the kind requires, by path relative to the root.
    groups_required: tuple[str, ...] = ()
    #: SPEC_00 section 3.1, return-type rule: `io.read` returns an `xarray.DataTree` for a
    #: grouped kind (C, N, raw) and an `xarray.Dataset` for the rest (T, D, G, R, W). The type
    #: is fixed by the kind, never by the content of a particular file.
    grouped: bool = False
    #: SPEC_00 section 5 requires `units`, `long_name` and `provenance` on every variable.
    #: The `raw` kind is exempt by section 2.2.1: it requires only the section 5 globals.
    enforce_variable_attributes: bool = True
    notes: str = ""
    conditional_globals: dict[str, tuple[str, ...]] = field(default_factory=dict)


def _uncertainty_names(specs: tuple[VarSpec, ...]) -> tuple[str, ...]:
    return tuple(f"{s.name}_uncertainty" for s in specs if s.needs_uncertainty)


_THERMO = KindSpec(
    name="thermo",
    schema_version=1,
    dimensions=("level", "latitude_planetocentric"),
    required_dimensions=("level",),
    variables=(
        VarSpec("temperature_K", dims=None, units="K", needs_uncertainty=True),
        VarSpec("pressure_Pa", required=False, units="Pa", needs_uncertainty=True),
        # Source profiles only (SPEC_00 v0.16 section 6.1): the pressure as printed, beside a
        # `pressure_Pa` placed on the grid the global `pressure_grid_rule` declares.
        VarSpec("pressure_printed_Pa", required=False, units="Pa"),
        VarSpec("geopotential_m2s2", required=False, units="m2 s-2"),
        # Source profiles only. The measured quantity, above the datum in `height_datum`.
        VarSpec("height_m", required=False, units="m", needs_uncertainty=True),
        VarSpec("latitude_planetocentric_deg", required=False, units="degrees_north"),
    ),
    globals_required=("vertical_coordinate", "thermo_instance"),
    conditional_globals={"source_profile": SOURCE_PROFILE_GLOBALS},
    notes="SPEC_00 section 6.1. No composition and no geodesy in this file.",
)

_COMPOSITION = KindSpec(
    name="composition",
    schema_version=1,
    dimensions=("level", "latitude_planetocentric"),
    required_dimensions=("level",),
    variables=(
        VarSpec("pressure_Pa", required=False, units="Pa"),
        VarSpec("geopotential_m2s2", required=False, units="m2 s-2"),
        VarSpec("latitude_planetocentric_deg", required=False, units="degrees_north"),
    ),
    globals_required=("composition_role", "vertical_coordinate"),
    groups_required=("species",),
    grouped=True,
    conditional_globals={"reduction": ("source_statement",)},
    notes=(
        "SPEC_00 section 6.2. At least one x_<species> variable is required, and the mole "
        "fractions must sum to one within 1e-9 at every point."
    ),
)

_GEODESY = KindSpec(
    name="geodesy",
    schema_version=1,
    dimensions=("surface",),
    required_dimensions=("surface",),
    variables=(
        VarSpec("surface_pressure_Pa", dims=("surface",), units="Pa"),
        VarSpec("radius_equatorial_m", dims=("surface",), units="m", needs_uncertainty=True),
        VarSpec("radius_polar_m", dims=("surface",), units="m", needs_uncertainty=True),
        VarSpec("oblateness", dims=("surface",), units="1", needs_uncertainty=True),
        VarSpec("fit_residual_m", required=False, dims=("surface",), units="m"),
    ),
    globals_required=("citation", "fit_inputs", "fit_latitude_convention"),
    notes="SPEC_00 section 6.3. A reduction phase input only; `forward` refuses it.",
)

_GRAVITY = KindSpec(
    name="gravity",
    schema_version=1,
    dimensions=("degree",),
    required_dimensions=("degree",),
    variables=(
        VarSpec("degree", dims=("degree",), units="1"),
        VarSpec("J", dims=("degree",), units="1", needs_uncertainty=True),
        VarSpec("J_status", dims=("degree",), units="1", is_flag=True),
        VarSpec("GM_m3s2", dims=(), units="m3 s-2", needs_uncertainty=True),
        VarSpec("normalization_radius_m", dims=(), units="m"),
    ),
    globals_required=("GM_scope", "epoch", "harmonic_convention"),
    notes="SPEC_00 section 6.4. J values are unscaled, the printed value times 1e-6.",
)

_ROTATION = KindSpec(
    name="rotation",
    schema_version=1,
    variables=(
        VarSpec("period_s", dims=(), units="s"),
        VarSpec("angular_rate_rad_s", dims=(), units="rad s-1"),
    ),
    globals_required=("system_name", "epoch", "citation"),
    notes="SPEC_00 section 6.5.",
)

_WIND = KindSpec(
    name="wind",
    schema_version=1,
    dimensions=("latitude_planetocentric", "pressure"),
    required_dimensions=("latitude_planetocentric", "pressure"),
    variables=(
        VarSpec("latitude_planetocentric_deg", dims=("latitude_planetocentric",), units="degrees_north"),
        VarSpec("pressure_Pa", dims=("pressure",), units="Pa"),
        VarSpec("u_total_ms", units="m s-1", needs_uncertainty=True),
        VarSpec("u_cylindrical_ms", units="m s-1"),
        VarSpec("u_shear_ms", units="m s-1"),
        VarSpec("value_provenance", units="1", is_flag=True),
        VarSpec("reference_level_pressure_Pa", dims=(), units="Pa"),
    ),
    globals_required=(
        "rotation_system_name",
        "rotation_rate_rad_s",
        "epoch",
        "method",
        "observation_level_Pa",
        "observation_level_justification",
        "source_latitude_convention",
        "vertical_structure",
        "coverage_pressure_Pa",
        "coverage_latitude_planetocentric_deg",
        "decomposition_geometry",
    ),
    notes=(
        "SPEC_00 section 6.6. The sum identity u_total = u_cylindrical + u_shear is checked "
        "on read. The constancy of Omega_abs on cylinders is a model load time check that "
        "needs the run geometry and so is not performed here."
    ),
)

_REFRACTIVITY = KindSpec(
    name="refractivity",
    schema_version=1,
    dimensions=("level",),
    required_dimensions=("level",),
    variables=(
        VarSpec("radius_m", dims=("level",), units="m", needs_uncertainty=True),
        VarSpec("height_above_anchor_isobar_m", dims=("level",), units="m",
                needs_uncertainty=True),
        VarSpec("number_density_m3", dims=("level",), units="m-3", needs_uncertainty=True),
        VarSpec("refractivity", dims=("level",), units="1", needs_uncertainty=True),
        # SPEC_00 section 6.7 v0.13: every consumer needs these two, with their companions.
        VarSpec("mean_refractivity_m3", dims=("level",), units="m3", needs_uncertainty=True),
        VarSpec("mean_molar_mass_kg_mol", dims=("level",), units="kg mol-1",
                needs_uncertainty=True),
        # "Scalars with uncertainties" (SPEC_00 section 6.7): each of the six has a companion.
        VarSpec("latitude_planetocentric_deg", dims=(), units="degrees_north", needs_uncertainty=True),
        VarSpec("psi_deg", dims=(), units="degrees", needs_uncertainty=True),
        VarSpec("latitude_planetographic_deg", dims=(), units="degrees_north",
                needs_uncertainty=True),
        VarSpec("anchor_isobar_pressure_Pa", dims=(), units="Pa", needs_uncertainty=True),
        VarSpec("anchor_isobar_radius_m", dims=(), units="m", needs_uncertainty=True),
        VarSpec("anchor_isobar_height_m", dims=(), units="m", needs_uncertainty=True),
    ),
    groups_required=(
        "inputs/thermo",
        "inputs/composition",
        "inputs/geodesy",
        "inputs/gravity",
        "inputs/rotation",
        "inputs/wind",
        # SPEC_00 section 6.7: both are groups with no variables, carrying attributes only.
        # `manifest` holds the reduction manifest verbatim in `text` and its hash in `sha256`.
        "manifest",
        "reduction_record",
    ),
    grouped=True,
    notes="SPEC_00 section 6.7. `role` is always `reduction` for this kind.",
)

_RAW = KindSpec(
    name="raw",
    schema_version=1,
    enforce_variable_attributes=False,
    grouped=True,
    notes=(
        "SPEC_00 section 2.2.1. No standard schema beyond the section 5 globals, by design: "
        "its purpose is a complete and checkable transcription, not a model input."
    ),
)

KINDS: dict[str, KindSpec] = {
    k.name: k
    for k in (
        _THERMO,
        _COMPOSITION,
        _GEODESY,
        _GRAVITY,
        _ROTATION,
        _WIND,
        _REFRACTIVITY,
        _RAW,
    )
}


def kind_spec(kind: str) -> KindSpec:
    """Return the specification for `kind`, refusing an unknown name."""
    try:
        return KINDS[kind]
    except KeyError:
        known = ", ".join(sorted(KINDS))
        raise CasspianSchemaError(
            f"unknown kind {kind!r}; the kinds of SPEC_00 section 6 are: {known}"
        ) from None


# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------


def check_names(dataset, where: str = "dataset") -> None:
    """Refuse a bare `latitude`, `r`, `z`, or `p` as a variable, dimension, or attribute.

    SPEC_00 section 5. Enforced for every kind, `raw` included: the rule exists so that no
    number anywhere in the system carries a name without its convention or its unit.
    """
    for name in list(dataset.variables) + list(dataset.dims):
        if str(name) in FORBIDDEN_BARE_NAMES:
            raise CasspianSchemaError(
                f"{where}: the name {name!r} carries no convention and no unit. "
                "SPEC_00 section 5 requires latitude_planetocentric_deg or "
                "latitude_planetographic_deg, radius_m or height_m, pressure_Pa."
            )
    for name in dataset.attrs:
        if str(name) in FORBIDDEN_BARE_NAMES:
            raise CasspianSchemaError(
                f"{where}: the global attribute {name!r} carries no convention and no unit "
                "(SPEC_00 section 5)."
            )
    for var_name, var in dataset.variables.items():
        for name in var.attrs:
            if str(name) in FORBIDDEN_BARE_NAMES:
                raise CasspianSchemaError(
                    f"{where}: variable {var_name!r} has an attribute named {name!r}, which "
                    "carries no convention and no unit (SPEC_00 section 5)."
                )


def _check_globals(dataset, spec: KindSpec, writer_filled: bool, where: str) -> None:
    attrs = dataset.attrs
    required = list(AUTHOR_REQUIRED_GLOBALS) + list(spec.globals_required)
    if writer_filled:
        required += list(WRITER_FILLED_GLOBALS)
        required.append("history")
    # `codata_release` is required only of a file that used a constant of nature, and
    # `input_hashes` only of a derived file (SPEC_00 section 5), so neither is demanded here.
    for name in required:
        if name == "codata_release":
            continue
        if name not in attrs:
            raise CasspianSchemaError(
                f"{where}: required global attribute {name!r} is missing "
                f"(kind {spec.name}, SPEC_00 section 5 and 6)."
            )
    role = attrs.get("role")
    if role is not None and role not in ROLES:
        raise CasspianSchemaError(
            f"{where}: role is {role!r}; SPEC_00 section 5 allows {sorted(ROLES)}."
        )
    if spec.name == "refractivity" and role != "reduction":
        raise CasspianSchemaError(
            f"{where}: kind N is always role 'reduction' (SPEC_00 section 6.7), not {role!r}."
        )
    if spec.name == "gravity":
        stated = attrs.get("harmonic_convention")
        if stated not in KNOWN_HARMONIC_CONVENTIONS:
            known = ", ".join(sorted(KNOWN_HARMONIC_CONVENTIONS))
            raise CasspianSchemaError(
                f"{where}: harmonic_convention is {stated!r}, which this reader does not know "
                f"(SPEC_00 section 6.4). The codes it knows are: {known}. A different "
                "convention gets a new code, never a reinterpretation of an old one."
            )
    if spec.name == "thermo":
        instance = attrs.get("thermo_instance")
        if instance not in THERMO_INSTANCES:
            raise CasspianSchemaError(
                f"{where}: thermo_instance is {instance!r}; SPEC_00 section 6.1 allows "
                f"{sorted(THERMO_INSTANCES)}."
            )
        if "pressure_printed_Pa" in dataset.variables and "pressure_grid_rule" not in attrs:
            raise CasspianSchemaError(
                f"{where}: pressure_printed_Pa is present, so the global attribute "
                "'pressure_grid_rule' must say which grid pressure_Pa was placed on "
                "(SPEC_00 section 6.1 v0.16)."
            )
        if instance == "source_profile":
            for name in SOURCE_PROFILE_GLOBALS:
                if name not in attrs:
                    raise CasspianSchemaError(
                        f"{where}: thermo_instance is 'source_profile', so the global "
                        f"attribute {name!r} is required (SPEC_00 section 6.1)."
                    )
            if not any(
                str(n).startswith("latitude_") and str(n).endswith("_deg") for n in attrs
            ):
                raise CasspianSchemaError(
                    f"{where}: a source profile states its latitude as a global attribute in "
                    "the source's convention, such as latitude_planetographic_deg "
                    "(SPEC_00 section 6.1). None is present."
                )
        else:
            present = [name for name in SOURCE_PROFILE_GLOBALS if name in attrs]
            if present:
                raise CasspianSchemaError(
                    f"{where}: thermo_instance is {instance!r}, so the source profile "
                    f"attributes {present} must be absent (SPEC_00 section 6.1). A "
                    "retrieval and a delivered field carry latitude as a coordinate."
                )
    vertical = attrs.get("vertical_coordinate")
    if vertical is not None and vertical not in VERTICAL_COORDINATES:
        raise CasspianSchemaError(
            f"{where}: vertical_coordinate is {vertical!r}; "
            f"SPEC_00 section 6 allows {sorted(VERTICAL_COORDINATES)}."
        )


def _check_dimensions(dataset, spec: KindSpec, where: str) -> None:
    present = set(dataset.dims)
    for name in spec.required_dimensions:
        if name not in present:
            raise CasspianSchemaError(
                f"{where}: required dimension {name!r} is missing (kind {spec.name})."
            )
    for name in spec.dimensions:
        if name in present or name in spec.required_dimensions:
            continue
        marker = f"{name}_absent_meaning"
        meaning = dataset.attrs.get(marker)
        if meaning is None:
            raise CasspianSchemaError(
                f"{where}: dimension {name!r} is absent and {marker!r} does not say why. "
                "SPEC_00 section 5: absent, assumed, and missing are three different things."
            )
        if meaning not in ABSENT_MEANINGS:
            raise CasspianSchemaError(
                f"{where}: {marker} is {meaning!r}; SPEC_00 section 5 allows "
                f"{sorted(ABSENT_MEANINGS)}."
            )
        if meaning == "point" and name == "latitude_planetocentric":
            companions = (
                "latitude_planetocentric_deg",
                "latitude_planetographic_deg",
            )
            if not any(c in dataset.attrs or c in dataset.variables for c in companions):
                raise CasspianSchemaError(
                    f"{where}: latitude is declared 'point' but the single coordinate value "
                    f"is stated nowhere. SPEC_00 section 5 requires one of {companions}."
                )


def _check_variables(dataset, spec: KindSpec, where: str) -> None:
    names = set(dataset.variables)
    for var in spec.variables:
        if var.required and var.name not in names:
            raise CasspianSchemaError(
                f"{where}: required variable {var.name!r} is missing (kind {spec.name})."
            )
        if var.needs_uncertainty and var.name in names:
            companion = var.companion()
            if companion not in names:
                raise CasspianSchemaError(
                    f"{where}: {companion!r} is missing. SPEC_00 section 5 requires the "
                    "uncertainty companion to exist even when the source states none, "
                    "holding NaN."
                )
    if spec.name == "thermo":
        vertical = dataset.attrs.get("vertical_coordinate")
        if vertical is not None and vertical not in names:
            raise CasspianSchemaError(
                f"{where}: vertical_coordinate declares {vertical!r} but no such variable "
                "is present."
            )
    if spec.name == "composition":
        if not any(str(n).startswith("x_") and not str(n).endswith("_uncertainty") for n in names):
            raise CasspianSchemaError(
                f"{where}: kind C requires at least one x_<species> variable "
                "(SPEC_00 section 6.2)."
            )


def _check_variable_attributes(dataset, spec: KindSpec, where: str) -> None:
    if not spec.enforce_variable_attributes:
        return
    for name, var in dataset.variables.items():
        for required in ("units", "long_name", "provenance"):
            if required not in var.attrs:
                raise CasspianSchemaError(
                    f"{where}: variable {name!r} has no {required!r} attribute "
                    "(SPEC_00 section 5)."
                )
        provenance = var.attrs["provenance"]
        if provenance not in PROVENANCE_VALUES:
            raise CasspianSchemaError(
                f"{where}: variable {name!r} has provenance {provenance!r}; SPEC_00 "
                f"section 5 allows {sorted(PROVENANCE_VALUES)}."
            )
        # The companion may end in `_uncertainty` (section 5) or carry the unit suffix after
        # it (`GM_uncertainty_m3s2`, section 6.4), so match on the token, not the ending.
        if "_uncertainty" in str(name):
            kind_attr = var.attrs.get("uncertainty_kind")
            if kind_attr is None:
                raise CasspianSchemaError(
                    f"{where}: {name!r} has no uncertainty_kind attribute "
                    "(SPEC_00 section 5)."
                )
            if kind_attr not in UNCERTAINTY_KINDS:
                raise CasspianSchemaError(
                    f"{where}: {name!r} has uncertainty_kind {kind_attr!r}; SPEC_00 "
                    f"section 5 allows {sorted(UNCERTAINTY_KINDS)}."
                )


def _check_coordinates(dataset, spec: KindSpec, where: str) -> None:
    if not spec.enforce_variable_attributes:
        return
    import numpy as np

    for name in dataset.dims:
        if name not in dataset.variables:
            continue
        values = np.asarray(dataset[name].values)
        if values.ndim != 1 or values.size < 2:
            continue
        if not np.issubdtype(values.dtype, np.number):
            continue
        increasing = bool(np.all(np.diff(values) > 0))
        decreasing = bool(np.all(np.diff(values) < 0))
        if not (increasing or decreasing):
            raise CasspianSchemaError(
                f"{where}: coordinate {name!r} is not strictly monotonic "
                "(SPEC_00 section 5)."
            )
        attrs = dataset[name].attrs
        # SPEC_00 section 5 asks every coordinate for a `positive` attribute. That is a CF
        # attribute for vertical coordinates, and it has no meaning on a harmonic degree, a
        # species index, or a surface index, so it is required only of the vertical ones.
        # See REPORT_01_step1.md section 6, item 5: this narrowing awaits a ruling.
        if not any(hint in str(name) for hint in VERTICAL_COORDINATE_HINTS):
            continue
        positive = attrs.get("positive")
        if positive is None or positive not in POSITIVE_VALUES:
            raise CasspianSchemaError(
                f"{where}: vertical coordinate {name!r} needs a 'positive' attribute of "
                f"{sorted(POSITIVE_VALUES)} so a reader never guesses (SPEC_00 section 5)."
            )
        if "pressure" in str(name):
            direction = attrs.get("direction")
            if direction is None or direction not in DIRECTION_VALUES:
                raise CasspianSchemaError(
                    f"{where}: pressure coordinate {name!r} needs a 'direction' attribute "
                    f"of {sorted(DIRECTION_VALUES)} (SPEC_00 section 5)."
                )


def check_mole_fractions(dataset, where: str, tolerance: float = 1.0e-9) -> None:
    """Refuse a kind C dataset whose mole fractions do not sum to one.

    SPEC_00 section 6.2 sets the tolerance at 1e-9 and section 8 makes the reader refuse.
    """
    import numpy as np

    fractions = [
        dataset[n]
        for n in dataset.variables
        if str(n).startswith("x_") and not str(n).endswith("_uncertainty")
    ]
    if not fractions:
        return
    total = sum(f for f in fractions)
    worst = float(np.max(np.abs(np.asarray(total.values) - 1.0)))
    if worst > tolerance:
        raise CasspianSchemaError(
            f"{where}: mole fractions depart from one by {worst:.3e}, which exceeds the "
            f"tolerance {tolerance:.1e} of SPEC_00 section 6.2."
        )


def check_wind_components(dataset, where: str) -> None:
    """Refuse a kind W dataset whose components do not sum to its total.

    SPEC_00 section 6.6: the identity holds at every point to round-off, or the file is
    refused. The tolerance scales with the magnitude of the field.
    """
    import numpy as np

    needed = ("u_total_ms", "u_cylindrical_ms", "u_shear_ms")
    if not all(n in dataset.variables for n in needed):
        return
    total = np.asarray(dataset["u_total_ms"].values, dtype=float)
    rebuilt = np.asarray(dataset["u_cylindrical_ms"].values, dtype=float) + np.asarray(
        dataset["u_shear_ms"].values, dtype=float
    )
    scale = max(float(np.nanmax(np.abs(total))), 1.0)
    worst = float(np.nanmax(np.abs(total - rebuilt)))
    if worst > 1.0e-12 * scale:
        raise CasspianSchemaError(
            f"{where}: u_total does not equal u_cylindrical + u_shear; worst departure "
            f"{worst:.3e} m/s against a field scale of {scale:.3e} m/s "
            "(SPEC_00 section 6.6)."
        )


def check_wind_poles(dataset, where: str) -> None:
    """Refuse a kind W dataset whose wind is not exactly zero at both poles.

    SPEC_00 section 6.6 (v0.8). A zonal wind is zero at a pole by definition, and a nonzero
    value there makes `Omega_abs` and the meridional gravity singular: the centrifugal part of
    `G_phi` carries `u^2 tan(phi) / r`, which diverges (REPORT_01_step7, finding 5). Both poles
    must be nodes of the latitude coordinate, and the tool that writes the file, not the model,
    is where a source that does not reach the pole is brought to zero.
    """
    import numpy as np

    name = "latitude_planetocentric_deg"
    if name not in dataset.variables or "u_total_ms" not in dataset.variables:
        return
    latitude = np.asarray(dataset[name].values, dtype="float64")
    for pole in (-90.0, 90.0):
        match = np.flatnonzero(np.abs(latitude - pole) <= 1e-9)
        if match.size == 0:
            raise CasspianSchemaError(
                f"{where}: {pole:+.0f} degrees is not a node of the latitude coordinate. "
                "SPEC_00 section 6.6 requires both poles to be nodes of a kind W file."
            )
        values = np.asarray(dataset["u_total_ms"].values, dtype="float64")[match]
        worst = float(np.max(np.abs(values)))
        if worst != 0.0:
            raise CasspianSchemaError(
                f"{where}: u_total_ms at {pole:+.0f} degrees is {worst:.3e} m/s, not exactly "
                "zero (SPEC_00 section 6.6). A zonal wind is zero at a pole by definition, and "
                "a nonzero value makes the meridional gravity singular there."
            )


def validate(dataset, kind: str, where: str = "dataset", writer_filled: bool = True) -> None:
    """Validate `dataset` against `kind`, raising `CasspianSchemaError` on the first fault.

    `writer_filled` is False when validating a dataset before `write` has filled in the
    globals that are the writer's job; it is True for anything read back from disk.
    """
    spec = kind_spec(kind)
    check_names(dataset, where)
    _check_globals(dataset, spec, writer_filled, where)
    _check_dimensions(dataset, spec, where)
    _check_variables(dataset, spec, where)
    _check_variable_attributes(dataset, spec, where)
    _check_coordinates(dataset, spec, where)
    if kind == "composition":
        check_mole_fractions(dataset, where)
    if kind == "wind":
        check_wind_components(dataset, where)
        check_wind_poles(dataset, where)
