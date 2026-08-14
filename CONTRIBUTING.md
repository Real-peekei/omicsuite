# Contributing to omicsuite

omicsuite is meant to grow past its current four modules. This document is
the contract every pipeline follows, so a new one -- whatever the analysis --
fits in without a redesign.

## The pipeline contract

Every pipeline in this package is a single entry-point function that:

1. **Fits a model.** Named `fit_<method>_pipeline()` (e.g.
   `fit_coxph_pipeline()`), `simulate_<method>()` for stochastic simulation
   (e.g. `simulate_gillespie_epidemic()`), or `integrate_<method>()` for
   data-integration methods (e.g. `integrate_multiomics()`). Pick whichever
   verb matches what the function actually does.

2. **Runs the diagnostics a careful methods section would include** --
   not just whether the model fit, but whether its assumptions hold, whether
   any single observation is driving the result, and whether the fitted
   values are internally consistent. Build these on top of the base
   R/Bioconductor/CRAN package doing the fitting; don't assume that
   package's own convenience diagnostics are enough to build a verdict on
   without checking what they actually compute.

3. **Returns one S3 object** (`structure(list(...), class = "<name>_pipeline")`)
   containing at minimum:
   - the fitted model object itself
   - a `plots` named list of `ggplot2` objects (use `theme_omicsuite()` from
     `R/utils.R`)
   - a `verdicts` data.frame -- see below

4. **Implements `print()`, `summary()`, and `plot()` S3 methods** for that
   class. `print()` should call `print_verdicts()` on `x$verdicts` as its
   last step. `plot()` should take a `which` argument defaulting to
   `names(x$plots)`.

## The verdicts data.frame

Every verdict row has the same six columns: `check`, `verdict`, `statistic`,
`p_value`, `note`, `plot`. `verdict` is one of:

| Label | Built with | Means |
|---|---|---|
| `"pass"` / `"flagged"` | `make_verdict(check, passed, ...)` | An actual pass/fail assumption check (e.g. proportional hazards, MCMC convergence) |
| `"info"` | `make_note(check, "info", ...)` | A fitted-value note that isn't a check (e.g. dispersion parameter, R0) |
| `"review"` | `make_note(check, "review", ...)` | Needs visual/manual inspection, not auto-scored (e.g. functional form, posterior predictive check) |
| `"interpretation"` | `make_note(check, "interpretation", ...)` | A plain-language restatement of a fitted effect and its uncertainty |

`plot` is an optional character column (`NA` by default) naming the entry
in the pipeline's `plots` list that a row explains, e.g. `plot = "ph_plot"`
for the proportional hazards rows. Set it when a verdict is genuinely about
one specific plot; leave it `NA` otherwise. **Plots themselves never carry
captions or embedded interpretation text** -- the `plot` column is the only
link between a figure and its explanation, so a plot always renders clean
and the explanation always lives in one place (`verdicts`), not duplicated
onto the figure.

**On interpretation rows specifically:** these restate what the model
estimated (an effect size, a CI, a significance threshold) in plain
language -- they are not a biological, clinical, or causal claim about *why*
the effect exists. "A one-unit increase in `age` is associated with a 3.2%
higher hazard" is an interpretation row. "Older patients have worse outcomes
because of X" is not something a verdict should ever say, because the model
doesn't know that. Stick to restating the statistic; leave the causal story
to the person reading it.

Prefer `check` names in the form `check_type[term]` (e.g.
`"proportional_hazards[age]"`, `"interpretation[armtreatment]"`) so a
consumer of the verdicts table can `grepl()` for a specific term across
checks.

## Dependency policy

- **The core four's dependencies** (`survival`, `ggplot2`, `stats`, `rlang`,
  `RGCCA`) stay hard `Imports` -- they're either foundational (used by
  multiple modules) or were added when the package was small enough that
  hard-importing an ordinary CRAN package was harmless.
- **Every module beyond those four goes in `Suggests`, regardless of
  toolchain weight.** This changed from the original policy (which allowed
  hard `Imports` for any toolchain-free CRAN package): with dozens of
  planned modules (see `ROADMAP.md`), hard-importing each one's specific
  packages would eventually force everyone to install the entire roadmap
  just to use `fit_coxph_pipeline()`. `survminer` and `flexsurv` are pure R
  with no compiler toolchain, and they're still `Suggests` -- the rule is
  about keeping the base install lean as the module count scales, not about
  toolchain difficulty specifically.
- The function guards on availability at call time via a small internal
  `<pkg>_is_available()` wrapper around `requireNamespace()` (see
  `brms_is_available()` in `R/utils.R` for the pattern -- note the
  wrapper-function trick needed to make it mockable in tests). If a module
  needs more than one optional package, guard on each independently so a
  partial install still gets as far as possible before erroring.
- Either way: the survival module must keep working for someone who installs
  omicsuite without any of the heavier or module-specific dependencies.

## Testing

- Structural tests (return type, expected fields, error messages) should run
  fast and always. Use small toy data.
- Tests that fit a real, slow model (Stan, bootstrap-heavy) should be gated
  with `skip_on_cran()` and, where relevant, `skip_if_not_installed()`.
- Write a regression test for every real bug you find, named after what it
  guards against -- see `test-coxph_pipeline.R`'s NA-row-alignment test for
  the pattern.

## Vignettes

Each module gets its own vignette under `vignettes/`, added to the navbar
menu in `_pkgdown.yml`. If the fitting method is too slow to run on every
package build (e.g. Stan), set `eval = FALSE` on the model-fitting chunks
and say so explicitly in the vignette text -- see
`vignettes/rnaseq-nb-pipeline.Rmd`.
