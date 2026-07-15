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
