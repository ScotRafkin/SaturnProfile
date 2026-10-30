"""The file conventions of SPEC_00 section 5 and the file kinds of section 6, as code.

This module is the single validation layer called for in SPEC_00 section 8. It declares what
every kind requires and it refuses anything that does not conform. It performs no file access;
`lib.io` opens and closes files and calls `validate` here.

Manuscript equations implemented: none. This module encodes file structure only.

SPEC_R1: a requirement here stays only if without it the code would compute a wrong number
without saying so, or fail with no clear message. So a kind requires the variables, dimensions,
groups and globals some code reads, and the checks that guard the arithmetic (monotonic
coordinates, mole fractions that sum to one, a known harmonic convention). Descriptive globals,
variable attributes, labels, naming rules and the season identifier are written by the tools and
required by nothing.
"""

from __future__ import annotations

from dataclasses import dataclass

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

#: SPEC_00 v0.17 section 5: the `epoch` of a forward product whose run declared a season but no date.
PROFILE_EPOCH_WITHOUT_DATE = "none: season declared as solar longitude only"

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
    #: Dimensions the kind defines. Only `required_dimensions` is checked (SPEC_R1).
    dimensions: tuple[str, ...] = ()
    #: Dimensions that must be present and may not be omitted.
    required_dimensions: tuple[str, ...] = ()
    variables: tuple[VarSpec, ...] = ()
    #: Global attributes this kind requires: the ones some code reads (SPEC_R1).
    globals_required: tuple[str, ...] = ()
    #: Groups the kind requires, by path relative to the root.
    groups_required: tuple[str, ...] = ()
    #: SPEC_00 section 3.1, return-type rule: `io.read` returns an `xarray.DataTree` for a
    #: grouped kind (C, N, raw) and an `xarray.Dataset` for the rest (T, D, G, R, W). The type
    #: is fixed by the kind, never by the content of a particular file.
    grouped: bool = False
    #: Whether coordinates are checked strictly monotonic. The `raw` kind is a transcription
    #: and is exempt (SPEC_00 section 2.2.1).
    enforce_variable_attributes: bool = True
    notes: str = ""


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
    notes="SPEC_00 section 6.1. No composition and no geodesy in this file.",
)

_COMPOSITION = KindSpec(
    name="composition",
    schema_version=1,
    dimensions=("level", "latitude_planetocentric"),
    required_dimensions=("level", "latitude_planetocentric"),
    variables=(
        VarSpec("pressure_Pa", required=False, units="Pa"),
        VarSpec("geopotential_m2s2", required=False, units="m2 s-2"),
        VarSpec("latitude_planetocentric_deg", dims=("latitude_planetocentric",),
                units="degrees_north"),
    ),
    groups_required=("species",),
    grouped=True,
    notes=(
        "SPEC_00 section 6.2 as SPEC_04 decision O amends it: kind C is one structure for "
        "every use, a field on (level, latitude_planetocentric) with pressure as the vertical "
        "coordinate, uniform in latitude when that is the hypothesis. The one-latitude form "
        "and its latitude_planetocentric_absent_meaning marker are retired; a reader that "
        "wants a column asks for it by latitude (lib.composition.column_at, decision L). At "
        "least one x_<species> variable is required, and the mole fractions must sum to one "
        "within 1e-9 at every point."
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
    globals_required=("harmonic_convention",),
    notes="SPEC_00 section 6.4. J values are unscaled, the printed value times 1e-6.",
)

_ROTATION = KindSpec(
    name="rotation",
    schema_version=1,
    variables=(
        VarSpec("period_s", dims=(), units="s"),
        VarSpec("angular_rate_rad_s", dims=(), units="rad s-1"),
    ),
    globals_required=("system_name",),
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
        VarSpec("u_total_ms", required=False, units="m s-1", needs_uncertainty=True),
        VarSpec("u_reference_ms", dims=("latitude_planetocentric",), units="m s-1"),
        VarSpec("u_shear_ms", required=False, units="m s-1"),
        VarSpec("value_provenance", units="1", is_flag=True),
        VarSpec("reference_level_pressure_Pa", dims=(), units="Pa"),
    ),
    globals_required=(
        "rotation_system_name",
        "rotation_rate_rad_s",
        "vertical_structure",
        "coverage_pressure_Pa",
        "coverage_latitude_planetocentric_deg",
    ),
    notes=(
        "SPEC_00 section 6.6 as SPEC_08 amends it. The file is a wind at a stated level, on "
        "the source's own latitude and pressure grid: the reference level wind "
        "u_reference_ms(latitude) at reference_level_pressure_Pa, both required, and, where the "
        "source gives them, the total u_total_ms(latitude, pressure) and the shear u_shear_ms, "
        "the change of the wind from the reference level along the local vertical. Nothing ties "
        "the optional parts to each other, and neither the poles nor a range of pressure is "
        "required on read; the model's loaders require u_total_ms and apply the polar rule "
        "(check_wind_poles). `decomposition`, `u_cylindrical_ms`, `decomposition_geometry` and "
        "the Omega_abs cylinder check are retired (SPEC_04 decision A): the decomposition on "
        "a given geometry is a later diagnostic tool, and the model reads the total only."
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
    globals_required=("profile_or_run",),
    grouped=True,
    notes="SPEC_00 section 6.7. `profile_or_run` names the figures (tools/plots/render.py).",
)

_PROFILE = KindSpec(
    name="profile",
    schema_version=1,
    dimensions=("level", "latitude_planetocentric"),
    required_dimensions=("level",),
    variables=(
        VarSpec("geopotential_m2s2", dims=("level",), units="m2 s-2"),
        VarSpec("radius_m", dims=("level",), units="m"),
        VarSpec("height_above_anchor_isobar_m", dims=("level",), units="m"),
        VarSpec("refractivity", dims=("level",), units="1"),
        VarSpec("mean_refractivity_m3", dims=("level",), units="m3"),
        VarSpec("mean_molar_mass_kg_mol", dims=("level",), units="kg mol-1"),
        VarSpec("number_density_m3", dims=("level",), units="m-3"),
        VarSpec("pressure_Pa", dims=("level",), units="Pa"),
        VarSpec("temperature_K", dims=("level",), units="K"),
        # Closure mode only, and then both (SPEC_03 Step 3 deliverable 3).
        VarSpec("pressure_tabulated_Pa", required=False, dims=("level",), units="Pa"),
        VarSpec("temperature_tabulated_K", required=False, dims=("level",), units="K"),
        # Transfer mode only (SPEC_04 Step 5 deliverable 3). Absent in closure mode, where the
        # levels are the anchor's own and there is nothing to label, carry or place a datum from.
        VarSpec("pressure_label_Pa", required=False, dims=("level",), units="Pa"),
        VarSpec("pressure_identity_residual", required=False, dims=("level",), units="1"),
        VarSpec("refractivity_gauge", required=False, dims=("level",), units="1"),
        VarSpec("geopotential_gauge_m2s2", required=False, dims=("level",), units="m2 s-2"),
        VarSpec("altitude_m", required=False, dims=("level",), units="m"),
        VarSpec("gauge_latitude_planetocentric_deg", required=False, dims=(),
                units="degrees_north"),
        VarSpec("datum_isobar_Pa", required=False, dims=(), units="Pa"),
        VarSpec("datum_geopotential_m2s2", required=False, dims=(), units="m2 s-2"),
        VarSpec("datum_radius_m", required=False, dims=(), units="m"),
        VarSpec("reference_surface_radius_m", required=False, dims=(), units="m"),
        VarSpec("latitude_planetocentric_deg", dims=(), units="degrees_north"),
        VarSpec("psi_deg", dims=(), units="degrees"),
        VarSpec("gauge_isobar_Pa", dims=(), units="Pa"),
        VarSpec("gauge_level_index", dims=(), units="1"),
        VarSpec("boundary_pressure_Pa", dims=(), units="Pa"),
        VarSpec("boundary_level_index", dims=(), units="1"),
    ),
    globals_required=("profile_or_run",),
    groups_required=(
        "inputs/composition",
        "inputs/gravity",
        "inputs/rotation",
        "inputs/wind",
        "namelist",
    ),
    grouped=True,
    notes=(
        "SPEC_00 v0.16 section 6.8, SPEC_03 Step 3 deliverable 3, SPEC_04 Step 5 deliverable 3. "
        "Written by forward, role forward, one latitude; every anchor verbatim under "
        "anchors/<slug>; the run's record in production_record in closure mode and in "
        "transfer_record in transfer mode, exactly one of the two."
    ),
)

#: SPEC_04 Step 5 deliverable 3: the groups a transfer-mode product carries beyond the closure's.
#: `reference_surface` and `isobars` hold fields on the mesh's latitude nodes rather than on the
#: profile's levels, so they are groups and not variables of the root.
PROFILE_TRANSFER_GROUPS = ("estimate", "isobars", "reference_surface", "transfer_record")

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
        _PROFILE,
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


def _check_globals(dataset, spec: KindSpec, where: str) -> None:
    attrs = dataset.attrs
    for name in spec.globals_required:
        if name not in attrs:
            raise CasspianSchemaError(
                f"{where}: required global attribute {name!r} is missing (kind {spec.name})."
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


def _check_profile(dataset, where: str) -> None:
    """The rules of kind `profile` beyond its variable table (SPEC_03 Step 3 deliverable 3)."""
    import numpy as np

    names = set(dataset.variables)
    tabulated = {"pressure_tabulated_Pa", "temperature_tabulated_K"}
    present = tabulated & names
    if len(present) == 1:
        missing = (tabulated - present).pop()
        raise CasspianSchemaError(
            f"{where}: {present.pop()!r} is present without {missing!r}; the closure comparison "
            "carries both tabulated quantities or neither (SPEC_03 Step 3 deliverable 3)."
        )
    values = np.asarray(dataset["geopotential_m2s2"].values, dtype="float64")
    if values.ndim != 1 or not (np.all(np.diff(values) > 0) or np.all(np.diff(values) < 0)):
        raise CasspianSchemaError(
            f"{where}: the coordinate 'geopotential_m2s2' is not strictly monotonic (SPEC_03 Step 3 "
            "deliverable 3)."
        )


def _check_dimensions(dataset, spec: KindSpec, where: str) -> None:
    present = set(dataset.dims)
    for name in spec.required_dimensions:
        if name not in present:
            raise CasspianSchemaError(
                f"{where}: required dimension {name!r} is missing (kind {spec.name})."
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
    if spec.name == "profile":
        _check_profile(dataset, where)
    if spec.name == "composition":
        if not any(str(n).startswith("x_") and not str(n).endswith("_uncertainty") for n in names):
            raise CasspianSchemaError(
                f"{where}: kind C requires at least one x_<species> variable "
                "(SPEC_00 section 6.2)."
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


def check_wind_poles(dataset, where: str) -> None:
    """Refuse a kind W dataset whose wind is not exactly zero at both poles.

    SPEC_00 section 6.6 (v0.8), applied where the model takes a wind (the loaders of
    `lib.control`), not on read (SPEC_08). A zonal wind is zero at a pole by definition, and a nonzero
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


def validate(dataset, kind: str, where: str = "dataset") -> None:
    """Validate `dataset` against `kind`, raising `CasspianSchemaError` on the first fault."""
    spec = kind_spec(kind)
    _check_globals(dataset, spec, where)
    _check_dimensions(dataset, spec, where)
    _check_variables(dataset, spec, where)
    _check_coordinates(dataset, spec, where)
    if kind == "composition":
        check_mole_fractions(dataset, where)
