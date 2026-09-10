# REVIEW 01, Step 1

Review of `REPORT_01_step1.md` against SPEC_01 v0.3 Step 1. 10 September 2026.

**Disposition: accepted, with four changes to make before committing.** The code is read and is
clean: naming rule enforced, the three-state dimension rule implemented as specified, refusals
name the item and cite the section, no dashes, no attribution. The self-consistency numbers for
the constants are the right kind of thing to report and the right conclusion (truncation, not
error). The three specification defects were real. Rulings follow, and SPEC_00 v0.4 and SPEC_01
v0.4 carry them.

**Changes before commit.**

1. `HARMONIC_CONVENTION` becomes a code, `"CASSPIAN-J1"`, with the sentence defining it in the
   comment beside it (SPEC_00 §6.4 v0.4). A long descriptive string is a mismatch waiting to
   happen; a code is compared exactly and its meaning lives in one place. The kind G reader
   refuses a code it does not know.
2. `UNCERTAINTY_KINDS` returns to `{"1sigma", "range", "stated"}`. A bin standard deviation is a
   1σ estimate; the method goes in a free-text `uncertainty_method` attribute, which the
   registry allows on any uncertainty variable and never validates (SPEC_00 §5 v0.4).
3. `PROVENANCE_VALUES` gains `modeled`: produced by `casspian.forward`, appearing only in files
   the model writes (SPEC_00 §5 v0.4). Its purpose is to keep a delivered field from ever being
   mistaken for a retrieval when the two are laid side by side.
4. Kind T gains a required global `thermo_instance` in `{"source_profile", "retrieval",
   "model_output"}`; the source-profile globals are enforced when it is `source_profile` and
   forbidden otherwise (SPEC_00 §6.1 v0.4). This is the discriminator the report asked for;
   the guess it declined to make (point in latitude) was rightly declined.

**Decisions on the other items.**

5. `read` returning a Dataset for group-less kinds and a DataTree for kinds with groups is
   accepted, made a rule by kind rather than by file (SPEC_00 §3.1 v0.4): T, D, G, R, W give a
   Dataset; C, N, raw give a DataTree.
6. Kind N: `manifest` and `reduction_record` are groups with no variables and attributes only,
   both required (SPEC_00 §6.7 v0.4).
7. Tolerances as implemented stand. The tool's 1e-12 and the reader's 1e-9 are different
   numbers on purpose: the tool must be exact, the reader must not refuse a file over round-off
   from a different platform.
8. Type enforcement on write: accepted and welcome. Stricter than §8's list, and §8 is the
   floor, not the ceiling.
9. Acceptance scripts live in `reports/<step>/` at the repository root, ignored by git, one
   subdirectory per step holding the script and its printed output; kept while a step is under
   review, deletable afterward (SPEC_00 §2 v0.4). Add `reports/` to `.gitignore` in this
   step's commit. The physically meaningful checks of Steps 4 to 6 will be promoted somewhere
   permanent when the test discussion happens; that is not now.
10. `-dirty`: it is git's word for a working tree with uncommitted changes, so a hash with the
    suffix says "the code that ran was this commit plus something not recorded." Keeping the
    suffix is correct; recording a bare hash would be false. The rule (SPEC_00 §5 and §8 v0.4,
    SPEC_01 §0 v0.4): products written during a step may carry `-dirty`; after the acceptance
    commit the step's tool is rerun so every product on disk is clean; `refrac` and `forward`
    refuse any input carrying `-dirty`, tools do not.

**Then:** commit Step 1 (code plus the `.gitignore` line), rerun nothing (Step 1 writes no
product), set STATE row 1 to `accepted`, and begin Step 2.

---

## Addendum, after the changes were applied

Both further defects are real and both proposed readings are adopted (SPEC_00 v0.5).

11. `positive` is required only of vertical coordinates, defined by name (`pressure`,
    `geopotential`, `height`, `radius` in the coordinate name); other coordinates carry neither
    `positive` nor `direction`. Strict monotonicity stays on every coordinate.
12. Uncertainty companions are named by inserting `_uncertainty` before the unit suffix, and a
    variable with no unit suffix appends it. The section 6 tables were already written that way;
    section 5 now says it. The `VarSpec.uncertainty_name` mechanism is fine as the
    implementation, but with the rule stated the default should produce the right name without
    per-variable overrides; keeping the override for the two cases already set is acceptable.

Step 1 stands accepted. Commit it, set STATE row 1 to `accepted`, begin Step 2.
