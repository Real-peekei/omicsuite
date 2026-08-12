# omicsuite 0.7.0

* New: every plot across every module now carries its matching verdict
  note as a wrapped plot caption (via a new shared `add_caption()` helper
  in `R/utils.R`), so a figure alone -- not just the `verdicts` table --
  carries its own interpretation. Multiple matching verdicts (e.g. a KM
  plot's median-survival and log-rank notes) are joined with a blank line.
  `theme_omicsuite()` now styles captions left-aligned and italic, since
  ggplot2's right-aligned default reads poorly for multi-line paragraph
  text. Applied to all six modules: `fit_coxph_pipeline()` (`ph_plot`,
  `influence_plot`, `functional_form_<var>`), `fit_km_pipeline()`
  (`km_plot`, `parametric_overlay`), `fit_rnaseq_nb_pipeline()`
  (`pp_check`, `shrinkage_plot`, `rhat_plot`), `simulate_gillespie_epidemic()`
  (`trajectory_plot`, `final_size_hist`, `peak_time_hist`),
  `integrate_multiomics()` (`block_scores`, `variance_explained`,
  `stability`), and `fit_competing_risks_pipeline()` (`cif_plot`).

# omicsuite 0.6.0

* New: `fit_competing_risks_pipeline()`, the sixth module -- completes the
  "Survival analysis" roadmap row (Cox PH, Kaplan-Meier, competing risks):
  - nonparametric cumulative incidence functions (CIFs) per event type via
    `cmprsk::cuminc()`, optionally grouped
  - Gray's test for CIF differences across groups, per event type
  - a Fine-Gray subdistribution hazard model (`cmprsk::crr()`) with
    subdistribution hazard ratios reported the same way
    `fit_coxph_pipeline()` reports ordinary hazard ratios, but explicitly
    labeled so the two quantities don't get conflated
  - a standing `subdistribution_vs_cause_specific` verdict explaining why
    this differs from treating competing events as censoring in a
    standard Cox model
  - `event_var` here encodes *which* event occurred (0 = censored, 1 =
    event of interest, 2+ = competing events), documented prominently since
    it differs from the binary 0/1 convention used elsewhere in the package
  - `print()`, `summary()`, and `plot()` S3 methods
  - `cmprsk` is a `Suggests`, per the dependency policy
* Bug fix: `fit_competing_risks_pipeline()` errored with "arguments imply
  differing number of rows" whenever `group_var` was omitted. Root cause:
  `cmprsk::cuminc()` checks `missing(group)` internally to decide whether
  to fall back to its own default, but the code always passed `group = `
  (even when its value was `NULL`), so `missing()` was always `FALSE` and
  `as.factor(NULL)` produced a length-0 factor against the full data.
  Fixed by building the argument list dynamically via `do.call()`, so
  `group` is only present in the call at all when `group_var` is supplied.

# omicsuite 0.5.0

* New: `fit_km_pipeline()`, the fifth module -- Kaplan-Meier survival curve
  estimation, the descriptive complement to `fit_coxph_pipeline()`:
  - nonparametric survival curves, optionally stratified, with a log-rank
    test for group differences
  - median survival time with confidence interval per stratum, explicitly
    reporting "not reached" rather than a bare `NA` under heavy censoring
  - an optional parametric overlay (Weibull, exponential, etc. via
    `flexsurv`) with an AIC comparison and a visual-review verdict
  - a `survminer`-based KM plot with a risk table when `survminer` is
    installed, falling back to a plain `ggplot2` step-function plot
    (no risk table, linear rather than stair-stepped confidence bands)
    when it isn't
  - `print()`, `summary()`, and `plot()` S3 methods
* Changed dependency policy (see `CONTRIBUTING.md`): every module beyond the
  original four now goes in `Suggests` regardless of toolchain weight, not
  just Stan-style heavy dependencies. `survminer` and `flexsurv` are pure R
  with no compiler toolchain, and they're still `Suggests` -- with dozens of
  modules planned (`ROADMAP.md`), hard-`Imports`-ing each one's specific
  packages doesn't scale.
* Added `ROADMAP.md`, cataloging every pipeline/package from the user's
  biostatistics/bioinformatics/computational-biology reference list against
  the pipeline contract, with a status per entry.

# omicsuite 0.4.2

* Replicated the `"interpretation"` verdict pattern (introduced in the
  survival module in v0.4.1) to the other three modules:
  - `fit_rnaseq_nb_pipeline()`: a new `interpretation[<condition_coef>]` row
    restates the population-level `condition_var` effect as a fold-change
    with its 95% credible interval, notes whether that interval excludes no
    change, and names the genes with the strongest partially pooled
    evidence.
  - `simulate_gillespie_epidemic()`: the existing R0 and extinction-
    probability notes were relabeled from `"info"` to `"interpretation"`
    (they already restated a fitted quantity in plain language, just under
    the wrong category), and a new `interpretation[final_size]` row
    summarizes the final-size distribution (median, IQR) among realizations
    that became a sustained outbreak.
  - `integrate_multiomics()`: a new `alpha_group_separation` argument and,
    when `group` is supplied, a per-block `interpretation[group_separation_<block>]`
    row -- a Kruskal-Wallis test of whether each block's shared component
    actually separates the supplied group labels, not just a plot you have
    to eyeball.
  All three follow the same restraint documented in `CONTRIBUTING.md`:
  mechanical restatement of a fitted estimate and its uncertainty, never a
  causal or biological claim about why the effect exists.

# omicsuite 0.4.1

* Reframed the package as an extensible collection rather than fixed at four
  modules -- `Title`/`Description` now describe the shared pipeline contract
  (fit + diagnose + plot + verdict) so future additions across biostatistics,
  bioinformatics, and computational biology fit the same shape.
* New: `fit_coxph_pipeline()` verdicts now include an `"interpretation"` row
  per model term -- a plain-language restatement of the hazard ratio, its
  95% CI, and significance (e.g. "a one-unit increase in `age` is associated
  with a 3.2% higher hazard... statistically significant, p = 0.01"),
  distinguishing numeric covariates ("one-unit increase") from factor levels
  ("this level vs. the reference"). This is the template for adding the same
  interpretation layer to the other three modules.
* Fixed: `print_verdicts()` labeled every non-`"pass"` verdict as `"FLAG"`,
  including the pre-existing `"info"`/`"review"` rows in the RNA-seq and
  epidemic modules -- misleading, since those aren't failed checks. Now
  labeled distinctly (`INFO`, `RVW`, `INTP`, `FLAG`, `OK`).
* Added `make_note()`, an internal helper for info/review/interpretation
  verdict rows, replacing the `make_verdict(passed = NA)` + manual override
  pattern used ad hoc in earlier modules.
* Fixed: the pkgdown site failed to build (`reference[1].contents[2]
  (print.coxph_pipeline) must be a known topic name or alias`) because
  `print.*`/`summary.*` S3 methods across all four modules only had a bare
  `@export` tag with no title/description, so roxygen2 never generated a
  `.Rd` page for them -- `_pkgdown.yml` was referencing topics that didn't
  exist. Added full docblocks (title, `@param`, `@return`) to all eight
  `print.*`/`summary.*` methods, matching the documentation level already
  given to `plot.*`.
* Added a pkgdown documentation site (https://real-peekei.github.io/omicsuite/),
  auto-built and deployed to GitHub Pages via GitHub Actions on every push to
  `main`. Reference pages are grouped by module and all four vignettes are
  linked from the site navbar. CI deliberately does not install `brms` --
  see the comment in `.github/workflows/pkgdown.yaml`.

# omicsuite 0.4.0

* New: `integrate_multiomics()`, a multi-omics data integration pipeline via
  regularized generalized canonical correlation analysis (`RGCCA`), with:
  - automatic sample alignment across blocks by row name, with a verdict
    flagging any dropped samples
  - per-block average variance explained (AVE) per component
  - a case-resampling bootstrap loading-stability check (sign-agreement
    rate for each block's top-loading features), built independently of
    `RGCCA::rgcca_bootstrap()`'s internals
  - sample-score plots faceted by block, optionally colored by a supplied
    group/phenotype label
  - `print()`, `summary()`, and `plot()` S3 methods
  - `RGCCA` is a hard `Imports` (CRAN-only, no compiler toolchain needed,
    same risk profile as `survival`)

# omicsuite 0.3.0

* New: `simulate_gillespie_epidemic()`, a stochastic SIR/SEIR epidemic
  simulation pipeline via Gillespie's direct method, with:
  - exact stochastic simulation (no time-discretization approximation)
  - an ensemble summary (median + 5-95% envelope) across realizations
  - final outbreak size and peak-timing distributions across realizations
  - the implied basic reproduction number (`beta / gamma`)
  - the proportion of realizations that fade out early by chance --
    a distinction a deterministic ODE model can't show
  - `print()`, `summary()`, and `plot()` S3 methods
  - pure base R + ggplot2; no additional hard or suggested dependency

# omicsuite 0.2.0

* Bug fix: `fit_rnaseq_nb_pipeline()`'s `seed` argument defaulted to `NULL`,
  which `brms::brm()` cannot coerce to a numeric value (it expects `NA` or a
  real number). Default changed to `NA`, matching `brms`'s own default.
* New: `fit_rnaseq_nb_pipeline()`, a Bayesian hierarchical negative binomial
  mixed model for longitudinal RNA-seq counts (via `brms`), with:
  - gene-level random slope for the effect of interest (partial pooling
    across genes rather than one GLM per gene)
  - optional subject-level random intercept for repeated-measures designs
  - a library-size offset computed automatically, or supplied directly
  - MCMC convergence diagnostics (Rhat, effective-sample-size ratio)
  - a dispersion check on the negative binomial shape parameter
  - a posterior predictive check plot
  - a shrinkage plot showing each gene's partially pooled effect estimate
  - `print()`, `summary()`, and `plot()` S3 methods
  - `brms` is a `Suggests`, not a hard dependency, since it requires a Stan
    toolchain -- the survival module works without installing it.

# omicsuite 0.1.0.9000

* Bug fix: `fit_coxph_pipeline()` errored (`arguments imply differing number
  of rows`) when any model variable contained `NA` values, because `coxph()`
  silently drops incomplete rows but the influence and functional-form
  diagnostics were indexing against the original `data`. Diagnostics now
  track exactly which rows the model was fitted on (`fit$influence$used_rows`,
  `fit$influence$n_used`), and `flagged_rows` is reported in terms of the
  original `data`'s row numbers so `data[fit$influence$flagged_rows, ]` is
  always valid.

# omicsuite 0.1.0

* First release: `fit_coxph_pipeline()`, a high-level Cox proportional
  hazards pipeline with:
  - unadjusted-vs-adjusted model comparison
  - global and per-covariate proportional hazards testing (`survival::cox.zph()`)
  - influence diagnostics via dfbeta/dfbetas with a reported cutoff and flagged rows
  - functional form checks via martingale residuals for continuous covariates
  - `ggplot2` diagnostic plots and a structured verdict table
  - `print()`, `summary()`, and `plot()` S3 methods
* Planned for future releases: Bayesian hierarchical negative binomial mixed
  models for longitudinal RNA-seq, stochastic SIR/SEIR epidemic simulation via
  the Gillespie algorithm, and multi-omics data integration.
