# omicsuite

Portfolio-grade analytical pipelines for biostatistics and bioinformatics.
Each pipeline in this package returns a fitted model *plus* the diagnostic
checks, plots, and a written verdict you'd want before reporting the
result -- the same level of rigor as a methods section, not just a point
estimate.

**Status:** v0.2.0 ships the survival and RNA-seq modules. Stochastic
SIR/SEIR epidemic simulation and multi-omics integration are planned for
subsequent releases (see `NEWS.md`).

## Installation

```r
# install.packages("devtools")
devtools::install_github("Real-peekei/omicsuite")
```

Or, from a local clone:

```r
devtools::document()   # regenerates NAMESPACE and man/ pages from roxygen comments
devtools::install()
```

## Quick example

```r
library(omicsuite)

set.seed(1)
n <- 300
dat <- data.frame(
  time  = rexp(n, rate = 0.05),
  event = rbinom(n, 1, 0.7),
  arm   = factor(sample(c("control", "treatment"), n, replace = TRUE)),
  age   = rnorm(n, 55, 10)
)

fit <- fit_coxph_pipeline(
  data       = dat,
  time_var   = "time",
  event_var  = "event",
  covariates = "arm",
  adjust_for = "age"
)

print(fit)                 # models + full verdict table
plot(fit, which = "ph_plot")
plot(fit, which = "functional_form_age")
fit$verdicts                # data.frame you can knit straight into a report
```

## What `fit_coxph_pipeline()` checks

| Check | Method | Where to look |
|---|---|---|
| Overall effect | Unadjusted vs. adjusted `coxph` fit | `fit$model_unadjusted`, `fit$model_adjusted` |
| Proportional hazards | `survival::cox.zph()`, global + per-term | `fit$ph_test`, `fit$plots$ph_plot` |
| Influence | dfbeta/dfbetas, flagged rows beyond `2/sqrt(n)` | `fit$influence`, `fit$plots$influence_plot` |
| Functional form | Martingale residuals vs. each continuous covariate | `fit$functional_form`, `fit$plots$functional_form_<var>` |

All of the above are also summarized in `fit$verdicts`, a data.frame with one
row per check, a pass/flagged/review status, the relevant statistic, and a
plain-language note -- written to be pasted directly into an article or
report.

## Development notes

This package was scaffolded without a local R installation available, so the
`NAMESPACE` was hand-authored to match expected `roxygen2` output. Before
your first `R CMD check` or `devtools::check()`, run:

```r
devtools::document()
devtools::check()
```

and treat any diff in `NAMESPACE` or `man/` as the source of truth.

## RNA-seq: quick example

Requires `install.packages("brms")` (a Stan toolchain, not bundled with
`omicsuite` since it's a `Suggests` dependency -- the survival module
doesn't need it).

```r
library(omicsuite)

set.seed(1)
n_genes <- 20
n_samples <- 16
dat <- expand.grid(
  gene = paste0("gene_", seq_len(n_genes)),
  sample = paste0("sample_", seq_len(n_samples)),
  stringsAsFactors = FALSE
)
dat$condition <- factor(rep(rep(c("control", "treated"), each = n_samples / 2), n_genes))
dat$count <- rnbinom(nrow(dat), mu = 40, size = 4)

fit <- fit_rnaseq_nb_pipeline(
  data          = dat,
  gene_var      = "gene",
  count_var     = "count",
  condition_var = "condition",
  sample_var    = "sample"
)

fit$verdicts
plot(fit, which = "shrinkage_plot")
```

## What `fit_rnaseq_nb_pipeline()` checks

| Check | Method | Where to look |
|---|---|---|
| Convergence | Rhat + effective-sample-size ratio across all parameters | `fit$convergence`, `fit$plots$rhat_plot` |
| Dispersion | Posterior summary of the NB shape parameter | `fit$dispersion` |
| Posterior predictive fit | `brms::pp_check()` overlay (visual review) | `fit$plots$pp_check` |
| Per-gene effect | Partially pooled estimate vs. population-level estimate | `fit$shrinkage`, `fit$plots$shrinkage_plot` |

