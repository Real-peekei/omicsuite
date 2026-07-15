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
