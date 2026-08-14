# Fit a Kaplan-Meier Survival Pipeline

Fits a nonparametric Kaplan-Meier survival curve (optionally stratified
by a grouping variable), runs a log-rank test for group differences, and
reports median survival time with its confidence interval per stratum –
the descriptive complement to
[`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md)'s
modeling of covariate effects. Optionally overlays a parametric fit (via
[`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md))
for a visual and AIC-based check of whether a parametric distribution
(Weibull, exponential, etc.) is a reasonable summary of the same data.

## Usage

``` r
fit_km_pipeline(
  data,
  time_var,
  event_var,
  strata_var = NULL,
  parametric_dist = NULL,
  conf_type = "log-log",
  alpha = 0.05
)
```

## Arguments

- data:

  A data.frame containing the survival data.

- time_var:

  Character. Name of the time-to-event column.

- event_var:

  Character. Name of the event indicator column (1 = event occurred, 0 =
  censored).

- strata_var:

  Optional character. Name of a grouping variable for separate KM curves
  and a log-rank test. If `NULL` (the default), a single overall curve
  is fit with no group comparison.

- parametric_dist:

  Optional character naming a distribution accepted by
  [`flexsurv::flexsurvreg()`](http://chjackson.github.io/flexsurv-dev/reference/flexsurvreg.md)'s
  `dist` argument (e.g. `"weibull"`, `"exponential"`, `"gompertz"`,
  `"gengamma"`). If supplied, fits a parametric model of that family and
  overlays it on the KM curve.

- conf_type:

  Confidence interval type passed to
  [`survival::survfit()`](https://rdrr.io/pkg/survival/man/survfit.html)'s
  `conf.type`. Default `"log-log"` (bounded to `[0, 1]`, generally
  preferred over the plain `"plain"` linear CI).

- alpha:

  Significance threshold for the log-rank test verdict. Default `0.05`.

## Value

An object of class `"km_pipeline"`, a list with elements:

- fit:

  The `survfit` object.

- parametric_fit:

  The `flexsurvreg` object, or `NULL` if `parametric_dist` was not
  supplied.

- median_survival:

  A data.frame: one row per stratum (or one row named `"Overall"` if
  `strata_var` was not supplied), with median survival time and its
  confidence interval.

- logrank_test:

  The `survdiff` object, or `NULL` if `strata_var` was not supplied.

- plots:

  A named list of `ggplot`-family objects: `km_plot` (a `ggsurvplot`
  object if `survminer` is installed, otherwise a plain `ggplot`), and
  `parametric_overlay` if `parametric_dist` was supplied and the overlay
  could be built (see Details).

- verdicts:

  A data.frame: median survival per stratum and the log-rank test as
  `"interpretation"` rows, censoring summary and parametric-fit AIC as
  `"info"` rows, and the parametric overlay as a `"review"` row. A
  `plot` column names which entry in `plots` each row explains
  (`"km_plot"` or `"parametric_overlay"`).

## Details

`survminer` and `flexsurv` are `Suggests`, not hard dependencies (see
`CONTRIBUTING.md`'s dependency policy). Without `survminer` installed,
the KM plot falls back to a plain `ggplot2` step-function plot (no risk
table). `flexsurv` is only required when `parametric_dist` is supplied.
If the parametric overlay plot can't be built (e.g. an unexpected
`summary.flexsurvreg()` output shape on some `flexsurv` version), it's
silently omitted from `plots` with a warning rather than raising an
error – `parametric_fit` and its AIC verdict are unaffected either way.

## Examples

``` r
set.seed(1)
n <- 200
dat <- data.frame(
  time = rexp(n, rate = 0.05),
  event = rbinom(n, 1, 0.7),
  arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
)
fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm")
print(fit)
#> <omicsuite Kaplan-Meier pipeline>
#> 
#> Call: survfit(formula = survival::Surv(time, event) ~ arm, data = data, 
#>     conf.type = conf_type)
#> 
#>                 n events median 0.95LCL 0.95UCL
#> arm=control    94     69   22.7    15.0    25.9
#> arm=treatment 106     69   20.7    15.7    26.3
#> 
#> Diagnostic verdicts:
#> [INTP] interpretation[median_survival_control] Median survival for `control`: 22.74 (95% CI [15.03, 25.90]) among 94 subjects, 69 events. [see plots$km_plot]
#> [INTP] interpretation[median_survival_treatment] Median survival for `treatment`: 20.70 (95% CI [15.70, 26.32]) among 106 subjects, 69 events. [see plots$km_plot]
#> [INTP] interpretation[log_rank_test] Log-rank test across `arm`: chi-squared = 0.44 (df = 1), p = 0.5083. No significant evidence that survival differs across groups at this alpha. [see plots$km_plot]
#> [INFO] censoring_summary            62 of 200 subjects (31.0%) were censored rather than observed to have the event. [see plots$km_plot]
```
