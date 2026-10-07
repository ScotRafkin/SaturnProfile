# REVIEW 10, Step 1. The fit of the IRIS temperatures and their gradient

Reviewing agent, 7 October 2026. Report: `reports/REPORT_10_step1.md`. Specification:
`docs/specs/SPEC_10_IRIS_Fit.md` v0.3.

**Accepted.** `fit` reproduces linear and quadratic fields to machine precision, and its standard
error matches the scatter of the slopes for the linear fits.

**The author's choice (7 October): (a), local linear, at least 3 effective points.** These become
`fit`'s defaults. (c) invents structure at the gap edges. (a) and (b) agree where the data are dense,
and (a) keeps the sparse structure with its standard error. The reviewing agent's trials of a local
quadratic with a higher minimum and of a quadratic and linear blend gave fuller peaks but nearly the
same gradient; they can be revisited.

The extra `fwhm_deg` field is kept.

Order of work: set the defaults; rerun `step09_1` check 2 with its overlay closed; the author's go;
the commit of the code, the suite, the report and this review. No regression.
