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
