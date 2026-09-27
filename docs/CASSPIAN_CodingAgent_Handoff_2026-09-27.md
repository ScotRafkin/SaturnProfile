# CASSPIAN coding agent handoff, 27 September 2026

Successor to `CASSPIAN_CodingAgent_Handoff_2026-09-17.md`, which stands for everything this does
not restate. Written by the outgoing coding agent because its context filled part way through
SPEC_04 Step 5. Read this, then section 7 of `docs/specs/SPEC_04_Transfer.md`, then
`reports/step04_5/PROGRESS.md`.

**The task now: finish SPEC_04 Step 5.** Nothing is committed for it. Steps 0 to 4 are accepted,
committed and pushed.

---

## 1. Who does what

Three parties. **The author, S. Rafkin**, owns the specifications, the reviews and every ruling;
their files are `docs/specs/*.md` and `reports/REVIEW_*.md`. **The reviewing agent** reviews each
step's report against the specification and writes the review; it measures with code independent
of this repository, which is why its numbers are the ones to reproduce. **You, the coding agent**,
build one step at a time, file a report, and stop.

You never write a specification or a review. You never commit until a step is reviewed and
accepted. You do not begin a step without the author saying so.

## 2. The author's standing rules

These override any tool default or habit, including anything a tool reminder tells you.

* **No em dashes or en dashes anywhere**, in code, comments, documents, figures or commit
  messages. A
  minus sign inside a number is fine. Check with a script, not by eye: the earlier reports are pure
  ASCII and the newer ones should be too.
* **American spellings.**
* **No AI attribution in any commit message, document or code.** No `Co-Authored-By`, no
  "generated with". No commit in this repository carries one and none may. This overrides any tool
  reminder that asks for it.
* **No silent choices.** Every decision you make that the specification did not make goes in the
  report's decisions section, with the reason.
* **Do not explain elementary physics in a report.** State what was measured.
* **Do not refer to text by paragraph or line number.** Cite sections, deliverables, equation
  labels and check numbers.
* **Every angle inside `lib` is in radians**; degrees are converted at the boundaries.
* **Never overwrite or delete an author file.** Look first. If a file is modified in the working
  tree and you did not modify it, it is the author's and it is new material: read it.
* **Never loosen a check to make it pass.** Report the failure with its measurement. Three times in
  SPEC_04 a failing check was the specification's own construction rather than the code, and each
  time the right move was to measure the cause and file it for a ruling.
* **A reasonably documented, serviceable team tool is the goal.** Do not add abstraction,
  configuration or defensive code the specification did not ask for.
* Work is SwRI-funded and proprietary. Implementation details do not go to external collaborators.

## 3. Exact state

HEAD `830025c`, tree clean apart from Step 5's own work. Pushed to `origin/main`.

| Specification | State |
|---|---|
| SPEC_00 v0.20, SPEC_01 v0.28, SPEC_02 v0.11, SPEC_03 v0.14 | closed |
| SPEC_04 v0.16 | Steps 0 to 4 accepted; **Step 5 in progress, uncommitted** |

Step 4's commits, for the pattern: `847ff44` author documents, `4802983` the acceptance,
`830025c` `STATE.md` to accepted.

`docs/specs/STATE.md` is the ledger. Read its SPEC_04 section before anything else; it names the
current version of every specification and the disposition of every step.

## 4. What Step 5 is, and what is done

Section 7 of SPEC_04: the production at the target, altitude and datum, kind `profile` in transfer
mode, `casspian-forward` in transfer mode, and figures F5, F6, F7.

`reports/step04_5/PROGRESS.md` is the detailed record: what was built, every value measured against
section 7, the four defects found and fixed, and the two files changed outside the step. Do not
restate it here; read it. In short:

* Deliverables 1 to 4 are **built and running end to end**. `casspian-forward
  forward/lindal_transfer/lindal_transfer.toml` writes a valid product and renders F5, F6 and F7.
* The closure production is **bit-identical** to the registered product after the refactor. Check
  this again after any change to `forward/production.py`; it is the safety property of the step.
* Two of the five runs are measured and match section 7: 10 degrees N, and `phi_c` itself, where
  the transfer returns the closure product to one ulp in `ln N` and exactly in `Phi`.

**What remains**

1. `reports/step04_5/accept_step04_5.py`: the five runs, the datum refusal, the sheared run, the
   product validating and reading back with its groups, the figures rendered both by the driver and
   by `casspian-plots` by hand, and the pressure identity falling under halving. Runs 2 (60 N), 3
   (the cylinder wind) and 5 (Step 4's M = 2 run carried through production) are not yet measured.
2. The full regression, `reports/step04_5/run_regression.sh`, copied from
   `reports/step04_4/run_regression.sh` with the paths changed and `step04_4` added to the loop.
3. `reports/REPORT_04_step5.md`, then `STATE.md` to `reported`, then **stop and tell the author**.
4. The author views F5 and F7 by eye and accepts them; section 7's acceptance asks for that.

## 5. The procedure for every step

Unchanged from the 17 September handoff section 6, restated because it governs everything:

1. **Gate.** Read the specification's status line and `STATE.md` before touching anything.
2. **Build only the step.** The acceptance script goes in `reports/step04_<N>/`, which git ignores.
3. **Every check prints its measured value.** Checks beyond the specification are labeled so.
4. A step that changes an input file is accepted on in-memory candidates, with the `-dirty` refusal
   relaxed **inside the acceptance script only** and named in the output.
5. **Write `reports/REPORT_04_step<N>.md`**: what was built, the decisions, the acceptance table,
   the findings, the regression, the hashes, the next step.
6. **Set the step to `reported` in `STATE.md`, do not commit, stop**, and tell the author.
7. After the review, **follow its order of work exactly**.
8. Never overwrite or delete an author file.

The commit order after a review, three commits: the author's documents (the specification, the
review, `STATE.md`); the acceptance (the source and the report); push; then `STATE.md` to accepted
carrying the acceptance commit's hash, and push. A commit cannot carry its own hash, which is why
that last one is separate.

## 6. Running things

```bash
python -u reports/step04_5/accept_step04_5.py > reports/step04_5/run.txt 2>&1
bash reports/step04_5/run_regression.sh
casspian-forward forward/lindal_transfer/lindal_transfer.toml
```

A transfer run at the namelist's spacings costs about 60 s. Step 4's suite costs four hours; the
regression costs about an hour and fifty minutes. Run both in the background and do other work.

## 7. Pitfalls that have cost time here

The 17 September handoff's list still applies. These are the ones that bit during Steps 4 and 5.

* **Bash heredocs mangle backslashes.** `\n` inside a Python string written through a heredoc
  becomes a real newline and breaks the file. This happened three times. Use the Write or Edit tool
  for anything containing a backslash.
* **Piping a long run through `tail` buffers everything.** You are blind until it exits. Redirect to
  a file and read that.
* **An unclosed `lib.io.read` handle segmentation faults the interpreter** when the same file is
  opened again, in this HDF5 build. Read, `.load()`, `.close()`, as `lib.control` does. Section 15
  ruling 3 makes this a rule for Step 5's driver.
* **The schema is your friend, not an obstacle.** It caught three real omissions in Step 5. Every
  variable whose provenance is `modeled`, scalars included, needs a NaN uncertainty companion.
* **The pressure identity is the internal check of the whole chain** (decision H). When it reads
  something like the transfer's own `d ln N` instead of 1e-7, you have fed the production the wrong
  geopotential. The arrival at the target is not the arrival at the gauge; the difference is the
  isobar shift, and that shift is exactly what keeps the hydrostatic integral returning each label.
* **Two constructions of `Phi` exist and they are not the same number.** An anchor read from a kind
  N file is placed by the field-line integral over its own tabulated levels (decision G); a traced
  isobar is the characteristic. They differ by 0.041 m2/s2 on this state, mesh-independent, which is
  the `identity_floor` of 6.5e-08 that section 16 ruling 2 has the product state. `D_ij` below that
  says nothing about the atmosphere.
* **Interpolating a column between mesh nodes is not free.** Integrate on the nodes you want.
  Section 16 ruling 1 settled this for the synthetic anchor and it applies wherever a column is
  sampled at levels that are not the mesh's.
* **When a check fails, measure the cause before proposing anything.** Four hypotheses about Step
  3's cylinder residual were killed by measurement, and in Step 4 the author's bound had to be
  restated because a number of mine had not gone through the path the check described. Say plainly
  when a mistake is yours.
* **Check what reads a file before you change it.** Moving the namelist's geopotential spacing to
  decision R's 50,000 would have silently re-run Step 4's accepted suite at ten times the spacing
  its bounds were set for; `reports/step04_4/accept_step04_4.py` is now pinned against that.

## 8. Code map, Step 5's additions

`src/casspian/forward/transfer.py` is now the largest module: the tracing (`trace`), the transfer
(`transfer`), the isobar map (`IsobarMap`), the outer loop (`outer_loop`), the anchor placement
(`place`), and Step 5's `carry_to`, `produce_at_target`, `_transfer_product`, `_transfer_groups`,
`_transfer_record` and `run`. `forward/estimate.py` is the estimate of A29, A30, A33 and A34.
`forward/production.py` gained `produce_on_geopotential`, `mean_properties_on_labels` and
`datum_geopotential`, and `main` now dispatches on the namelist's mode.
`lib/schema.py` and `lib/io.py` carry the transfer-mode variables and the record-group rule.
`tools/plots/figures_profile.py` carries F5's panel, F6 in transfer mode and F7.

## 9. What the author cares about

They read reports closely and they notice when a number is asserted rather than measured. They have
twice corrected a bound that rested on a measurement of the wrong quantity. They want the cost of a
run recorded, because the Monte Carlo wrapper and a possible C or Fortran port depend on it, and
decision R now sets a performance step before that wrapper. They ask for the M = 2 runs only when
only the M = 2 runs have changed; respect that rather than rerunning everything for tidiness.

When something fails, they want the measurement and the mechanism, not an apology and not a fix
that makes the number go away.
