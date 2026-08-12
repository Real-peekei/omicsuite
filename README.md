# omicsuite

[![pkgdown site](https://img.shields.io/badge/docs-pkgdown-3B6E8F)](https://real-peekei.github.io/omicsuite/)

An extensible collection of analytical pipelines for biostatistics,
bioinformatics, and computational biology. Every pipeline follows the same
contract: fit a model, run the diagnostic checks a careful methods section
would include, render `ggplot2` plots, and return a structured verdict
table -- assumption checks, informational notes, and a plain-language
interpretation of the fitted effect -- rather than a bare point estimate.
Every plot also carries its matching verdict note as a caption, so a
figure alone carries its own interpretation, not just the verdicts table.

**Status:** v0.7.0 ships six modules (survival Cox PH, survival
Kaplan-Meier, survival competing risks, RNA-seq, epidemic, multi-omics),
each with plain-language interpretation verdicts and captioned plots, plus
a documentation site. See `NEWS.md` for what's shipped and
[`ROADMAP.md`](ROADMAP.md) for what's next -- omicsuite is meant to grow
past these; see [`CONTRIBUTING.md`](CONTRIBUTING.md) for the shared
pipeline contract new modules follow. The
[package website](https://real-peekei.github.io/omicsuite/) has the full
reference and vignettes.

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

## Kaplan-Meier: quick example

```r
library(omicsuite)

set.seed(1)
n <- 200
dat <- data.frame(
  time  = rexp(n, rate = 0.05),
  event = rbinom(n, 1, 0.7),
  arm   = factor(sample(c("control", "treatment"), n, replace = TRUE))
)

fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm")
fit$verdicts
plot(fit, which = "km_plot")
```

## What `fit_km_pipeline()` checks

| Check | Method | Where to look |
|---|---|---|
| Median survival per stratum | `summary(survfit)$table`, explicit "not reached" handling | `fit$median_survival`, `fit$verdicts` |
| Group difference | Log-rank test | `fit$logrank_test`, `fit$verdicts` |
| Censoring | Proportion censored vs. observed | `fit$verdicts` (`censoring_summary`) |
| Parametric fit (optional) | AIC + visual overlay against the KM curve | `fit$parametric_fit`, `fit$plots$parametric_overlay` |

## Competing risks: quick example

```r
library(omicsuite)

set.seed(1)
n <- 300
dat <- data.frame(
  time  = rexp(n, rate = 0.1),
  # 0 = censored, 1 = event of interest, 2 = competing event
  event = sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3, 0.5, 0.2)),
  arm   = factor(sample(c("control", "treatment"), n, replace = TRUE))
)

fit <- fit_competing_risks_pipeline(
  dat, "time", "event", group_var = "arm", covariates = "arm"
)
fit$verdicts
plot(fit, which = "cif_plot")
```

## What `fit_competing_risks_pipeline()` checks

| Check | Method | Where to look |
|---|---|---|
| Event-type breakdown | Counts/proportions per event code | `fit$event_summary`, `fit$verdicts` |
| Group difference in CIF | Gray's test, per event type | `fit$cuminc_fit`, `fit$verdicts` |
| Subdistribution effect | Fine-Gray model (`cmprsk::crr()`) | `fit$crr_fit`, `fit$verdicts` |
| Cause-specific vs. subdistribution | Standing interpretive note | `fit$verdicts` (`subdistribution_vs_cause_specific`) |

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

## Epidemic: quick example

```r
library(omicsuite)

fit <- simulate_gillespie_epidemic(
  model = "SIR",
  initial_state = c(S = 999, I = 1, R = 0),
  params = list(beta = 0.4, gamma = 0.1),
  t_max = 100, n_sim = 100, seed = 1
)

print(fit)                        # R0, extinction rate, verdicts
plot(fit, which = "trajectory_plot")
plot(fit, which = "final_size_hist")
```

For SEIR, add an `E` compartment to `initial_state` and a `sigma` rate
(1/mean incubation period) to `params`:

```r
fit_seir <- simulate_gillespie_epidemic(
  model = "SEIR",
  initial_state = c(S = 999, E = 0, I = 1, R = 0),
  params = list(beta = 0.4, sigma = 0.2, gamma = 0.1),
  t_max = 150, n_sim = 100, seed = 1
)
```

## What `simulate_gillespie_epidemic()` checks

| Check | Method | Where to look |
|---|---|---|
| Basic reproduction number | `beta / gamma` (same threshold quantity for SIR and SEIR) | `fit$r0` |
| Stochastic extinction | Proportion of realizations with final size below a threshold fraction of the population | `fit$prop_extinct`, `fit$is_extinct` |
| Outbreak trajectory | Median + 5-95% envelope across realizations, with individual realizations overlaid | `fit$summary`, `fit$plots$trajectory_plot` |
| Final size / peak timing | Distributions across realizations | `fit$final_size`, `fit$peak_time`, `fit$plots$final_size_hist`, `fit$plots$peak_time_hist` |

## Multi-omics: quick example

```r
library(omicsuite)

set.seed(1)
n <- 60
sample_ids <- paste0("patient_", seq_len(n))
transcriptomics <- matrix(rnorm(n * 30), nrow = n, dimnames = list(sample_ids, paste0("gene_", 1:30)))
proteomics <- matrix(rnorm(n * 15), nrow = n, dimnames = list(sample_ids, paste0("protein_", 1:15)))

fit <- integrate_multiomics(
  blocks = list(transcriptomics = transcriptomics, proteomics = proteomics),
  ncomp = 2, n_boot = 30, seed = 1
)

fit$verdicts
plot(fit, which = "block_scores")
plot(fit, which = "stability")
```

## What `integrate_multiomics()` checks

| Check | Method | Where to look |
|---|---|---|
| Sample alignment | Row-name intersection across blocks | `fit$verdicts` (`sample_alignment`) |
| Variance explained | Per-block AVE from `RGCCA::rgcca()` | `fit$variance_explained`, `fit$plots$variance_explained` |
| Loading stability | Case-resampling bootstrap sign-agreement rate, top loadings per block | `fit$stability`, `fit$plots$stability` |
| Sample structure | Scores on shared components, faceted by block, optionally colored by group | `fit$scores`, `fit$plots$block_scores` |

