# REVIEW 01, Step 0

Review of `REPORT_01_step0.md` against SPEC_01 v0.2 Step 0. 10 September 2026.

**Disposition: accepted.** Both commits and the tag are as specified; `git ls-tree -r
phase1-prototype` listing the force-added netCDF settles the point that caused the revision.

**On the deviations.**

1. `src/casspian/__init__.py`: an omission in the spec's enumeration, correctly read.
2. `__version__` in that file: keep it. Rule from here on: `pyproject.toml` and `__version__`
   must agree, and Step 1's `write()` reads the package version through `importlib.metadata`
   so there is one source at run time; the literal in `__init__.py` is for humans. If the two
   ever disagree, `write()` refuses.
3. No attribution trailer: correct, and it stays that way.
4. `data_static/` and `occul_data/` inside the `phase1-prototype` tag: acceptable. The tag
   records the working tree on the day; that it also holds current transcriptions is harmless
   and is now written down here.
5. `.pytest_cache/` removal: fine.

**On the questions.**

1. Line endings: agreed, and it is done before Step 1 writes anything. Step 1 now opens with a
   `.gitattributes` commit (SPEC_01 v0.3, Step 1, item 0) pinning every text transcription to
   LF and marking binaries, followed by `git add --renormalize .` so the committed bytes and
   the on-disk bytes agree on every platform. Provenance hashes are of the bytes on disk.
2. `raw/notes.md`: the author's deliverable, not the coding agent's. It is now in
   `occul_data/lindal/raw/`. Step 2 treats it as documentation: not read by the tool, its hash
   recorded in the raw bundle's `raw_sources` alongside the two data files.
3. Species table and GM values: both literature items were closed on 10 September (every
   number the Lindal chain needs has a value and a status; primary-source confirmations of
   values already read from a secondary source remain, none of which changes a value).
   SPEC_00 §9 is updated to say so.
4. Noted.

**Step 1 may begin.** STATE.md row 0 is set to `accepted` by this review.
