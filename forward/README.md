# Forward run directories

One directory per forward run, named by the user, as specified in SPEC_00 section 2.3 [D2a].

```
forward/<run>/
    <run>.toml             the run namelist (SPEC_00 section 7.2)
    <run>_build.toml       the run's build control file: pointers and declared choices for the tools
    inputs/                the run's own input files, written by the tools
        <run>_composition.nc    kind C
        <run>_gravity.nc        kind G
        <run>_rotation.nc       kind R
        <run>_wind.nc           kind W
    output/                written by casspian.forward, and nothing else writes there
        <run>_profile.nc        kind profile
        figures/                the standard diagnostics, when the namelist asks for them
```

**Every input a run reads is its own file under `inputs/`**, built by the tools from
`<run>_build.toml` with `role = "forward"` and the run prefix, never a reduction file under
`occul_data/` and never a copy embedded in kind N, even when the content is identical (SPEC_00
section 2.3, SPEC_03 decision 7). The reduction's files record what a source assumed; a run's
inputs record what the run assumes. The namelist points at `inputs/` and at exactly one kind of
thing outside the run directory: the registered refractivity products it anchors on.

The four inputs are made by one command:

```
casspian-run-inputs forward/<run>/<run>_build.toml
```

Every path a file records is relative to that file's own directory (SPEC_00 v0.18 section 5), so a
run directory can be archived, copied or diffed as a unit. An ensemble is a directory of such
directories.

The namelist and the build file are committed. The netCDF under `inputs/` and everything under
`output/` are excluded by `.gitignore` and regenerable; the hashes of the inputs are recorded in
the report of the step that built them.

The runs here:

- `lindal_closure/`: the hydrostatic closure of the Lindal reduction at its own latitude, under
  forward inputs that match the source's assumptions (SPEC_03).
