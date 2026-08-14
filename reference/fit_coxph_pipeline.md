# Fit a Cox Proportional Hazards Pipeline with Full Diagnostics

Fits a Cox proportional hazards model and runs the diagnostic suite
you'd want in a methods section before trusting the hazard ratios: an
unadjusted-versus-adjusted comparison, global and covariate-level
proportional hazards testing, influence diagnostics, functional form
checks for continuous covariates via martingale residuals, and a
plain-language interpretation of every fitted hazard ratio.

## Usage

``` r
fit_coxph_pipeline(
  data,
  time_var,
  event_var,
  covariates,
  adjust_for = NULL,
  id_var = NULL,
  frailty = FALSE,
  alpha = 0.05,
  influence_cutoff_sd = 1
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

- covariates:

  Character vector of covariate names for the primary model of interest.

- adjust_for:

  Optional character vector of additional covariates to adjust for. When
  supplied, both an unadjusted model (`covariates` only) and an adjusted
  model (`covariates` + `adjust_for`) are fit, so you can see how much
  the estimates of interest move.

- id_var:

  Optional character. Cluster/subject identifier, for robust sandwich
  variance (`cluster()`) or a shared frailty term.

- frailty:

  Logical. If `TRUE` and `id_var` is supplied, fits a gamma frailty term
  instead of a robust cluster-adjusted variance. Default `FALSE`.

- alpha:

  Significance threshold used to flag diagnostic checks. Default `0.05`.

- influence_cutoff_sd:

  Numeric. Observations with a scaled dfbeta beyond this many "typical"
  units (`2 / sqrt(n)` scaled by this factor) are flagged as
  influential. Default `1` (i.e. the standard `2/sqrt(n)` cutoff).

## Value

An object of class `"coxph_pipeline"`, a list with elements:

- model_unadjusted:

  The unadjusted `coxph` fit.

- model_adjusted:

  The adjusted `coxph` fit, or `NULL` if `adjust_for` was not supplied.

- primary_model:

  Whichever of the two above is the model of record (adjusted if
  available, else unadjusted) – used for diagnostics.

- ph_test:

  The
  [`survival::cox.zph()`](https://rdrr.io/pkg/survival/man/cox.zph.html)
  result for `primary_model`.

- influence:

  A list with the dfbeta/dfbetas matrices (indexed to `used_rows`, in
  fitted-model order), `flagged_rows` (indices into the original `data`
  you passed in – so `data[fit$influence$flagged_rows, ]` always works,
  even if `coxph()` silently dropped rows with missing values),
  `used_rows` (which original rows the model was actually fitted on),
  and `n_used`.

- functional_form:

  A named list of data.frames (martingale residual vs. covariate value)
  for each continuous covariate in `covariates`.

- plots:

  A named list of `ggplot` objects: `ph_plot`, `influence_plot`, and one
  `functional_form_<var>` entry per continuous covariate.

- verdicts:

  A data.frame summarizing every diagnostic check
  (`"pass"`/`"flagged"`), each functional-form check (`"review"`), and a
  plain-language hazard-ratio interpretation for every model term
  (`"interpretation"`). A `plot` column names which entry in `plots` (if
  any) each row explains, e.g. `"ph_plot"` for the proportional hazards
  rows.

- alpha:

  The significance threshold used.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 200
dat <- data.frame(
  time = rexp(n, rate = 0.05),
  event = rbinom(n, 1, 0.7),
  age = rnorm(n, 55, 10),
  arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
)
fit <- fit_coxph_pipeline(
  data = dat, time_var = "time", event_var = "event",
  covariates = "arm", adjust_for = "age"
)
print(fit)
#> <omicsuite Cox PH pipeline>
#> 
#> Unadjusted model:
#> Call:
#> survival::coxph(formula = formula_unadjusted, data = data, x = TRUE)
#> 
#>                 coef exp(coef) se(coef)     z   p
#> armtreatment 0.04359   1.04455  0.17172 0.254 0.8
#> 
#> Likelihood ratio test=0.06  on 1 df, p=0.7996
#> n= 200, number of events= 138 
#> 
#> Adjusted model:
#> Call:
#> survival::coxph(formula = formula_adjusted, data = data, x = TRUE)
#> 
#>                  coef exp(coef) se(coef)     z     p
#> armtreatment 0.018610  1.018784 0.172798 0.108 0.914
#> age          0.010377  1.010431 0.008487 1.223 0.221
#> 
#> Likelihood ratio test=1.57  on 2 df, p=0.4551
#> n= 200, number of events= 138 
#> 
#> Diagnostic verdicts:
#> [OK  ] proportional_hazards[arm]    No evidence against proportional hazards for this term. [see plots$ph_plot]
#> [OK  ] proportional_hazards[age]    No evidence against proportional hazards for this term. [see plots$ph_plot]
#> [OK  ] proportional_hazards[GLOBAL] Global test finds no overall violation of the proportional hazards assumption. [see plots$ph_plot]
#> [FLAG] influence[dfbetas]           18 of 200 observations used in the fit (9.0%) exceed the |dfbetas| > 0.141 cutoff. A non-trivial share of observations are individually influential -- inspect flagged rows before reporting hazard ratios. [see plots$influence_plot]
#> [RVW ] functional_form[age]         Martingale residuals plotted against this covariate; a loess smooth that is flat suggests the linear (log-hazard) form is adequate. Requires visual review -- not auto-scored. [see plots$functional_form_age]
#> [INTP] interpretation[armtreatment] HR = 1.019 (95% CI [0.726, 1.429]). this level of `arm` (vs. the reference level) is associated with a 1.9% higher hazard at any given time, holding other model terms fixed (not statistically significant at this alpha, p = 0.9142).
#> [INTP] interpretation[age]          HR = 1.010 (95% CI [0.994, 1.027]). a one-unit increase in `age` is associated with a 1.0% higher hazard at any given time, holding other model terms fixed (not statistically significant at this alpha, p = 0.2214).
# }
```
