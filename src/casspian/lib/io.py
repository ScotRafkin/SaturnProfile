"""netCDF read and write for every kind, with the provenance rules of SPEC_00 sections 5 and 8.

`write` fills in the globals that are the writer's job, validates against the kind through
`lib.schema`, and writes netCDF-4. `read` opens, checks the declared kind against the kind
asked for, checks the schema version, validates, and returns the dataset. Neither guesses and
neither repairs: a file that does not conform is refused with a message naming the item at
fault.

Manuscript equations implemented: none. This module moves bytes and enforces conventions.

Groups. Kinds C (the `species` group), N (the `inputs/*`, `manifest` and `reduction_record`
groups) and `raw` (the transcription bundle) are group bearing. `write` takes them as a mapping
of group path to Dataset. The return-type rule of SPEC_00 section 3.1: `read` returns an
`xarray.DataTree` for those three kinds and an `xarray.Dataset` for T, D, G, R and W. The type
is fixed by the kind, never by the content of a particular file, so a caller can write against
it without opening the file first.
"""

from __future__ import annotations

import hashlib
import os
import subprocess
import warnings
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import xarray as xr

from casspian.lib import constants
from casspian.lib.schema import (
    CONVENTIONS,
    WRITER_FILLED_GLOBALS,
    CasspianSchemaError,
    kind_spec,
    validate,
)

_HASH_BLOCK = 1 << 20


# ---------------------------------------------------------------------------
# Provenance helpers
# ---------------------------------------------------------------------------


def sha256(path) -> str:
    """Return the SHA-256 of the bytes on disk at `path`, as lower case hex.

    SPEC_00 section 6, principle 6: provenance is by content hash. The hash is of the file as
    it sits on disk, which is why SPEC_01 Step 1 item 0 pins text files to LF.
    """
    digest = hashlib.sha256()
    with open(Path(path), "rb") as handle:
        for block in iter(lambda: handle.read(_HASH_BLOCK), b""):
            digest.update(block)
    return digest.hexdigest()


def recorded_path(path, product) -> str:
    """The path of `path` as the file at `product` records it (SPEC_00 v0.18 section 5).

    Relative to the directory of `product`, with `/` as the separator, so that a file and the
    inputs it names can be moved together and a record never carries a machine's paths. Every
    writer records every path it names through this function. A path that cannot be expressed
    relative to the product (another drive) is refused rather than recorded absolute.
    """
    target = Path(path).resolve()
    base = Path(product).resolve().parent
    try:
        relative = os.path.relpath(target, base)
    except ValueError:
        raise CasspianSchemaError(
            f"{target} cannot be recorded relative to {Path(product).name}: the two are on "
            "different drives. SPEC_00 v0.18 section 5 records every path relative to the file "
            "that names it; an absolute path is a defect."
        ) from None
    return Path(relative).as_posix()


def input_hash_entry(path, product) -> str:
    """One `input_hashes` entry, `<path relative to product> sha256:<hex>` (SPEC_00 section 5)."""
    return f"{recorded_path(path, product)} sha256:{sha256(path)}"


def input_hashes(paths, product) -> str:
    """An `input_hashes` global for the file at `product`: one entry per input, newline separated."""
    return "\n".join(input_hash_entry(p, product) for p in paths)


def hash_entries(text) -> dict:
    """Parse an `input_hashes` value into `{recorded path: hex digest}`."""
    entries = {}
    for line in str(text).splitlines():
        if " sha256:" in line:
            name, digest = line.rsplit(" sha256:", 1)
            entries[name.strip()] = digest.strip()
    return entries


def warn_changed_inputs(attrs, path) -> list:
    """Warn for every `input_hashes` entry that names an existing file whose hash differs.

    SPEC_00 v0.18 section 8, for every derived kind. Each recorded path is resolved against the
    directory of the file being read. A path that does not resolve to a file is silent, since a
    copied file is not a defect; the embedded copy or the recorded hash remains authoritative.
    Returns the messages, so a caller can report them.
    """
    path = Path(path)
    base = path.resolve().parent
    messages = []
    for name, digest in hash_entries(attrs.get("input_hashes", "")).items():
        candidate = base / name
        if candidate.is_file() and sha256(candidate) != digest:
            message = (
                f"{path.name}: input_hashes entry {name} does not match the file now at that path "
                "beside it. The embedded copy and the recorded hash remain authoritative (SPEC_00 "
                "section 8)."
            )
            warnings.warn(message, stacklevel=3)
            messages.append(message)
    return messages


def history_append(dataset, line: str):
    """Append one line to the append-only `history` global (SPEC_00 section 5).

    The dataset is modified in place and returned, so calls chain.
    """
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    entry = f"{stamp} {line}"
    existing = dataset.attrs.get("history")
    dataset.attrs["history"] = f"{existing}\n{entry}" if existing else entry
    return dataset


def package_version() -> str:
    """Return the installed package version, refusing if it disagrees with `__version__`.

    SPEC_01 v0.3 Step 1, the version rule: `pyproject.toml` and `casspian.__version__` must
    agree, `write` obtains the version through `importlib.metadata` so that there is one
    source at run time, and the literal in `__init__.py` is for humans.
    """
    from importlib.metadata import PackageNotFoundError
    from importlib.metadata import version as _installed_version

    import casspian

    literal = casspian.__version__
    try:
        installed = _installed_version("casspian")
    except PackageNotFoundError:
        raise CasspianSchemaError(
            "casspian is not installed, so no version can be recorded in created_by. "
            "Install it with 'pip install -e .'."
        ) from None
    if installed != literal:
        raise CasspianSchemaError(
            f"version disagreement: importlib.metadata reports {installed!r} but "
            f"casspian.__version__ is {literal!r}. SPEC_01 Step 1 requires them to agree; "
            "reinstall after changing pyproject.toml."
        )
    return installed


def git_commit(start=None) -> str:
    """Return the commit of the package checkout, for `casspian_git_commit`.

    A dirty working tree is reported as `<sha>-dirty`. That is not cosmetic: products are
    built before a step is accepted and therefore before its code is committed, so the
    distinction is what tells a reader whether the recorded commit fully describes the code
    that wrote the file.
    """
    here = Path(start) if start is not None else Path(__file__).resolve().parent
    try:
        sha = subprocess.run(
            ["git", "-C", str(here), "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(here), "status", "--porcelain"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return "unknown (not a git checkout)"
    return f"{sha}-dirty" if dirty else sha


# ---------------------------------------------------------------------------
# Type rules, SPEC_00 section 5
# ---------------------------------------------------------------------------

_ALLOWED_ATTR_TYPES = (str, int, float, np.integer, np.floating, np.ndarray, list, tuple)


def _check_attribute_types(dataset, where: str) -> None:
    for name, value in dataset.attrs.items():
        if isinstance(value, bool) or not isinstance(value, _ALLOWED_ATTR_TYPES):
            raise CasspianSchemaError(
                f"{where}: global attribute {name!r} has type {type(value).__name__}, which "
                "netCDF cannot store. SPEC_00 section 5 allows strings and numbers; write a "
                "flag as an int8 variable with flag_values and flag_meanings."
            )


def _check_types(dataset, kind: str, where: str) -> None:
    spec = kind_spec(kind)
    flags = {v.name for v in spec.variables if v.is_flag}
    for name, var in dataset.variables.items():
        dtype = var.dtype
        if str(name) in flags:
            if dtype != np.int8:
                raise CasspianSchemaError(
                    f"{where}: {name!r} is a status flag and must be int8 with flag_values "
                    f"and flag_meanings (SPEC_00 section 5); it is {dtype}."
                )
            for required in ("flag_values", "flag_meanings"):
                if required not in var.attrs:
                    raise CasspianSchemaError(
                        f"{where}: flag variable {name!r} has no {required!r} attribute "
                        "(SPEC_00 section 5, CF style)."
                    )
            continue
        if np.issubdtype(dtype, np.floating) and dtype != np.float64:
            raise CasspianSchemaError(
                f"{where}: {name!r} is {dtype}; SPEC_00 section 5 requires float64 for all "
                "physical quantities."
            )


# ---------------------------------------------------------------------------
# Write
# ---------------------------------------------------------------------------


def write(path, dataset, kind: str, groups=None, created_by: str = "casspian"):
    """Validate `dataset` against `kind`, fill the writer's globals, and write netCDF-4.

    `groups` is an optional mapping of group path (for example `"table1"` or
    `"inputs/thermo"`) to Dataset. Group datasets are written as they are given, after the
    naming rule of SPEC_00 section 5 has been applied to them.

    Refuses on any missing required item, naming it. Returns the resolved `Path` written.
    """
    path = Path(path)
    spec = kind_spec(kind)
    dataset = dataset.copy()

    version = package_version()
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

    filled = {
        "Conventions": CONVENTIONS,
        "casspian_kind": kind,
        "casspian_schema_version": int(spec.schema_version),
        "created_by": f"{created_by} {version}",
        "created_at": stamp,
        "casspian_git_commit": git_commit(),
        "codata_release": constants.CODATA_RELEASE,
    }
    for name, value in filled.items():
        dataset.attrs[name] = value

    history_append(dataset, f"written by {created_by} {version} as kind {kind}")

    _check_attribute_types(dataset, f"{path.name}")
    _check_types(dataset, kind, f"{path.name}")
    validate(dataset, kind, where=f"{path.name}", writer_filled=True)

    group_datasets = dict(groups or {})
    for group_path, group_ds in group_datasets.items():
        from casspian.lib.schema import check_names

        check_names(group_ds, f"{path.name}:{group_path}")
        _check_attribute_types(group_ds, f"{path.name}:{group_path}")

    for name in spec.groups_required:
        if name not in group_datasets:
            raise CasspianSchemaError(
                f"{path.name}: kind {kind} requires the group {name!r} (SPEC_00 section 6)."
            )
    if kind == "profile":
        _check_profile_anchor_groups(list(group_datasets), path.name)

    encoding = {}
    for name, var in dataset.data_vars.items():
        if np.issubdtype(var.dtype, np.floating):
            encoding[str(name)] = {"_FillValue": np.nan}

    path.parent.mkdir(parents=True, exist_ok=True)
    dataset.to_netcdf(path, mode="w", format="NETCDF4", engine="netcdf4", encoding=encoding)
    for group_path, group_ds in group_datasets.items():
        group_encoding = {
            str(n): {"_FillValue": np.nan}
            for n, v in group_ds.data_vars.items()
            if np.issubdtype(v.dtype, np.floating)
        }
        group_ds.to_netcdf(
            path,
            mode="a",
            format="NETCDF4",
            engine="netcdf4",
            group=group_path,
            encoding=group_encoding,
        )
    return path


# ---------------------------------------------------------------------------
# Read
# ---------------------------------------------------------------------------


def group_paths(path) -> list[str]:
    """List every group path in a netCDF file, depth first, as posix style paths."""
    import netCDF4

    found: list[str] = []

    def walk(node, prefix: str) -> None:
        for name, child in node.groups.items():
            here = f"{prefix}/{name}" if prefix else name
            found.append(here)
            walk(child, here)

    with netCDF4.Dataset(Path(path), "r") as root:
        walk(root, "")
    return found


def read(path, kind: str):
    """Open `path`, check it is `kind`, validate it, and return it.

    Returns an `xarray.DataTree` for the grouped kinds (C, N, raw) and an `xarray.Dataset` for
    the rest (T, D, G, R, W), per the return-type rule of SPEC_00 section 3.1. The type follows
    the kind, not the file. Refuses per SPEC_00 section 8: a declared kind that is not asked for, a
    schema version newer than this reader, a missing required item, a non monotonic
    coordinate, mole fractions that do not sum to one, an unknown harmonic convention, or
    wind components that do not sum to their total.
    """
    path = Path(path)
    spec = kind_spec(kind)
    where = path.name

    root = xr.open_dataset(path, engine="netcdf4")
    try:
        declared = root.attrs.get("casspian_kind")
        if declared is None:
            raise CasspianSchemaError(
                f"{where}: no casspian_kind attribute, so this is not a CASSPIAN file "
                "(SPEC_00 section 5)."
            )
        if declared != kind:
            raise CasspianSchemaError(
                f"{where}: declared kind is {declared!r} but {kind!r} was asked for "
                "(SPEC_00 section 8)."
            )
        stated_version = root.attrs.get("casspian_schema_version")
        if stated_version is None:
            raise CasspianSchemaError(
                f"{where}: no casspian_schema_version attribute (SPEC_00 section 5)."
            )
        stated_version = int(stated_version)
        if stated_version > spec.schema_version:
            raise CasspianSchemaError(
                f"{where}: schema version {stated_version} is newer than this reader, which "
                f"knows version {spec.schema_version} of kind {kind} (SPEC_00 section 8)."
            )
        for name in WRITER_FILLED_GLOBALS:
            if name == "codata_release":
                continue
            if name not in root.attrs:
                raise CasspianSchemaError(
                    f"{where}: required global attribute {name!r} is missing "
                    "(SPEC_00 section 5)."
                )
        validate(root, kind, where=where, writer_filled=True)
        present = group_paths(path)
        for name in spec.groups_required:
            if name not in present:
                raise CasspianSchemaError(
                    f"{where}: kind {kind} requires the group {name!r}, which is not in the "
                    "file (SPEC_00 section 6)."
                )
        if kind == "profile":
            _check_profile_anchor_groups(present, where)
    except Exception:
        root.close()
        raise

    if not spec.grouped:
        warn_changed_inputs(root.attrs, path)
        return root
    root.close()
    tree = xr.open_datatree(path, engine="netcdf4")
    if kind == "refractivity":
        try:
            _check_refractivity_hashes(tree, path)
        except Exception:
            tree.close()
            raise
    warn_changed_inputs(tree.attrs, path)
    return tree


def _check_profile_anchor_groups(groups, where: str) -> None:
    """Kind `profile` carries each anchor verbatim under `anchors/<slug>` (SPEC_03 Step 3)."""
    anchors = [g for g in groups if g.startswith("anchors/") and g.count("/") == 1]
    if not anchors:
        raise CasspianSchemaError(
            f"{where}: kind profile carries every anchor it used as a group anchors/<slug>, the "
            "kind N file verbatim; none is present (SPEC_03 Step 3 deliverable 3)."
        )


#: The inputs whose hashes kind N records, in `reduction_record` and in `input_hashes`.
_REFRACTIVITY_HASHED = ("thermo", "composition", "geodesy", "gravity", "rotation", "wind",
                        "manifest")


def _check_refractivity_hashes(tree, path: Path) -> None:
    """Refuse a kind N file whose `input_hashes` disagrees with its reduction record.

    SPEC_02 Step 4. `refrac` records the SHA-256 of each input twice: in the global
    `input_hashes` (SPEC_00 section 5) and as `input_sha256_<input>` in `reduction_record`. A
    global entry that no longer lists the recorded hash means the file was altered after it was
    written, and it is refused, naming the input. An entry that no longer matches the file now
    beside it only warns, by the rule `read` applies to every derived kind
    (`warn_changed_inputs`, SPEC_00 v0.18 section 8): the embedded copy is authoritative.
    """
    where = path.name
    record = tree["reduction_record"].attrs
    entries = hash_entries(tree.attrs.get("input_hashes", ""))
    listed = set(entries.values())
    for key in _REFRACTIVITY_HASHED:
        recorded = record.get(f"input_sha256_{key}")
        if recorded is None:
            raise CasspianSchemaError(
                f"{where}: reduction_record carries no input_sha256_{key} (SPEC_02 Step 4)."
            )
        if str(recorded) not in listed:
            raise CasspianSchemaError(
                f"{where}: input_hashes does not list the {key} hash the reduction recorded "
                f"({str(recorded)[:16]}...). The file has been altered since refrac wrote it."
            )
    if len(entries) != len(_REFRACTIVITY_HASHED):
        raise CasspianSchemaError(
            f"{where}: input_hashes lists {len(entries)} entries; kind N lists the six inputs "
            "and the manifest (SPEC_00 section 6.7)."
        )
