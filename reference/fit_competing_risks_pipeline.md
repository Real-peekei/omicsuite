# Fit a Competing Risks Pipeline

Fits nonparametric cumulative incidence functions (CIFs) for each event
type, optionally compares them across groups with Gray's test, and fits
a Fine-Gray subdistribution hazard model for the event type of primary
interest via [`cmprsk::crr()`](https://rdrr.io/pkg/cmprsk/man/crr.html).

## Usage

``` r
fit_competing_risks_pipeline(
  data,
  time_var,
  event_var,
  failcode = 1,
  cencode = 0,
  group_var = NULL,
  covariates = NULL,
  alpha = 0.05
)
```

## Arguments

- data:

  A data.frame containing the survival data.

- time_var:

  Character. Name of the time-to-event column.

- event_var:

  Character. Name of the event-type column (see Details above for its
  encoding – not the same as `event_var` elsewhere in `omicsuite`).

- failcode:

  The event code (as it appears in `data[[event_var]]`) representing the
  event of primary interest. Default `1`.

- cencode:

  The event code representing censoring. Default `0`.

- group_var:

  Optional character. Grouping variable for separate CIF curves and
  Gray's test. If `NULL` (the default), a single CIF per event type is
  estimated with no group comparison.

- covariates:

  Optional character vector of covariate names for the Fine-Gray model
  (`cov1` in
  [`cmprsk::crr()`](https://rdrr.io/pkg/cmprsk/man/crr.html)). Factors
  are dummy-coded automatically via
  [`stats::model.matrix()`](https://rdrr.io/r/stats/model.matrix.html).
  If `NULL` (the default), only the CIF/Gray's-test part of the pipeline
  runs – no subdistribution hazard model is fit.

- alpha:

  Significance threshold for Gray's test and the Fine-Gray model's
  verdicts. Default `0.05`.

## Value

An object of class `"competing_risks_pipeline"`, a list with elements:

- cuminc_fit:

  The `cuminc` object from
  [`cmprsk::cuminc()`](https://rdrr.io/pkg/cmprsk/man/cuminc.html).

- crr_fit:

  The `crr` object from
  [`cmprsk::crr()`](https://rdrr.io/pkg/cmprsk/man/crr.html), or `NULL`
  if `covariates` was not supplied.

- cif_summary:

  A long-format data.frame: time, cumulative incidence estimate, event
  type, and group (if `group_var` was supplied) – the tidy form of
  `cuminc_fit` used for plotting.

- event_summary:

  A data.frame: counts and proportions for each event code observed in
  `data[[event_var]]`, including censoring.

- plots:

  A named list of `ggplot` objects: `cif_plot`.

- verdicts:

  A data.frame: event-type breakdown as `"info"`, Gray's test as
  `"interpretation"` (when `group_var` is supplied), and subdistribution
  hazard ratios as `"interpretation"` per covariate term (when
  `covariates` is supplied). A `plot` column names which entry in
  `plots` each row explains.

## Details

**`event_var` means something different here than in
[`fit_km_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_km_pipeline.md)
or
[`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).**
Those treat `event_var` as a binary indicator (1 = event, 0 = censored).
Here it must encode *which* event type occurred: `0` = censored (or
whatever `cencode` is set to), `1` = the event of primary interest, `2`,
`3`, ... = competing event types that preclude the event of interest
from ever happening (e.g. death from another cause when the event of
interest is relapse). Reusing a plain 0/1 binary indicator here silently
treats every non-event as a competing risk of type... there is none,
which just reduces to the single-event case – not wrong, but usually not
what you meant to ask.

The reason a dedicated model matters: fitting a standard Cox model that
treats competing events as ordinary censoring overestimates the event of
interest's cumulative incidence, because it implicitly assumes subjects
who experienced a competing event would have gone on to experience the
event of interest had they lived. The Fine-Gray model avoids that
assumption by keeping competing-event subjects in the risk set on the
subdistribution timescale instead of removing them.

## Examples

``` r
# \donttest{
if (requireNamespace("cmprsk", quietly = TRUE)) {
  set.seed(1)
  n <- 300
  dat <- data.frame(
    time = rexp(n, rate = 0.1),
    # 0 = censored, 1 = event of interest, 2 = competing event
    event = sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3, 0.5, 0.2)),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
  fit <- fit_competing_risks_pipeline(
    dat, "time", "event", group_var = "arm", covariates = "arm"
  )
  print(fit)
}
#> <omicsuite competing risks pipeline>
#> 
#> Event summary:
#>   code             label   n proportion
#> 1    0          censored  92  0.3066667
#> 2    1 event of interest 146  0.4866667
#> 3    2   competing event  62  0.2066667
#> 
#> Fine-Gray subdistribution hazard model:
#> convergence:  TRUE 
#> coefficients:
#> armtreatment 
#>       0.2271 
#> standard errors:
#> [1] 0.1631
#> two-sided p-values:
#> armtreatment 
#>         0.16 
#> 
#> Diagnostic verdicts:
#> [INFO] event_type_breakdown         code 0 (censored): n = 92 (30.7%); code 1 (event of interest): n = 146 (48.7%); code 2 (competing event): n = 62 (20.7%)
#> [INFO] subdistribution_vs_cause_specific This pipeline models the subdistribution hazard (Fine-Gray), which is the right quantity for predicting cumulative incidence in the presence of competing risks. It answers a different question than a cause-specific hazard model (an ordinary Cox model with competing events treated as censored) -- the two can even point in different directions for the same covariate. Report which one you used and why. [see plots$cif_plot]
#> [INTP] interpretation[gray_test_event_1] Gray's test for event code 1 across `arm`: statistic = 1.70, p = 0.1925. No significant evidence that cumulative incidence differs across groups for this event type at this alpha. [see plots$cif_plot]
#> [INTP] interpretation[gray_test_event_2] Gray's test for event code 2 across `arm`: statistic = 0.79, p = 0.3746. No significant evidence that cumulative incidence differs across groups for this event type at this alpha. [see plots$cif_plot]
#> [INTP] interpretation[subdistribution_hr_armtreatment] Subdistribution HR for `armtreatment` = 1.255, p = 0.1600. Associated with a 25.5% higher subdistribution hazard of the event of interest (not statistically significant at this alpha). This is not the same quantity as a cause-specific hazard ratio from a standard Cox model -- see this function's documentation for why. [see plots$cif_plot]
# }
```
