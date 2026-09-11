"""Control file reading for the gravity and rotation tools.

SPEC_00 section 7: a control file holds pointers and declared choices only. A physical value in
a control file is a defect, and an unknown key is an error rather than a warning, because the
parser knows which keys exist. A relative path is resolved against the directory holding the
control file, and the resolved absolute form is what gets recorded in the product.

Manuscript equations implemented: none.

This is deliberately local to the tools. SPEC_00 section 3.1 lists `lib.control` as the home
for control file parsing once the reduction manifest and the run namelist exist; those are
Steps 9 and later, and building the general parser now would be implementing ahead. When
`lib.control` arrives it absorbs this module.
"""

from __future__ import annotations

import tomllib
from pathlib import Path


class ControlFileError(Exception):
    """A control file is malformed, incomplete, or carries a key the parser does not know."""


def load_section(path, section: str, allowed: dict[str, bool], path_keys=()) -> dict:
    """Load one section of a TOML control file and check its keys against `allowed`.

    `allowed` maps key name to whether it is required. Any key outside `allowed` is an error
    (SPEC_00 section 7).

    `path_keys` names the keys whose values are paths; each is resolved against the control
    file's own directory and returned as an absolute `Path`, per SPEC_00 section 7. The keys
    are named explicitly rather than guessed from the file extension, because guessing was
    wrong the first time a tool pointed at a `.csv`.
    """
    path = Path(path).resolve()
    if not path.exists():
        raise ControlFileError(f"{path}: control file does not exist")
    with open(path, "rb") as handle:
        document = tomllib.load(handle)

    if section not in document:
        raise ControlFileError(
            f"{path}: no [{section}] section. The sections this tool reads are named in "
            "SPEC_01 Step 3."
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
        raise ControlFileError(
            f"{path}: [{section}] is missing required key(s) {missing}."
        )

    resolved = {}
    for key, value in table.items():
        if key in set(path_keys) and isinstance(value, str):
            resolved[key] = (path.parent / value).resolve()
        else:
            resolved[key] = value
    return resolved


def reject_physical_values(table: dict, path, section: str) -> None:
    """Refuse a control file section that carries a bare number.

    SPEC_00 principle 2: data lives in netCDF, choices live in TOML, and a physical value in a
    control file is a defect. Every value this tool's sections accept is a path or a name.
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
