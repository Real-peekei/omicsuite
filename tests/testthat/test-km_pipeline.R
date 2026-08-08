test_that("fit_km_pipeline returns the expected structure (no strata)", {
  set.seed(1)
  n <- 150
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7)
  )
  fit <- fit_km_pipeline(dat, "time", "event")

  expect_s3_class(fit, "km_pipeline")
  expect_s3_class(fit$fit, "survfit")
  expect_null(fit$logrank_test)
  expect_null(fit$parametric_fit)
  expect_identical(nrow(fit$median_survival), 1L)
  expect_identical(fit$median_survival$stratum, "Overall")
  expect_true(is.data.frame(fit$verdicts))
  expect_true(any(grepl("^interpretation\\[median_survival_Overall\\]$", fit$verdicts$check)))
  expect_true(any(fit$verdicts$check == "censoring_summary"))
  expect_true("km_plot" %in% names(fit$plots))
})

test_that("fit_km_pipeline with strata_var produces a log-rank test and per-stratum medians", {
  set.seed(2)
  n <- 200
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
  fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm")

  expect_identical(nrow(fit$median_survival), 2L)
  expect_setequal(fit$median_survival$stratum, c("control", "treatment"))
  expect_s3_class(fit$logrank_test, "survdiff")
  expect_true(any(fit$verdicts$check == "interpretation[log_rank_test]"))
  logrank_row <- fit$verdicts[fit$verdicts$check == "interpretation[log_rank_test]", ]
  expect_false(is.na(logrank_row$p_value))
})

test_that("censoring_summary reports the correct proportion", {
  dat <- data.frame(
    time = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
    event = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0)
  )
  fit <- fit_km_pipeline(dat, "time", "event")
  cens_row <- fit$verdicts[fit$verdicts$check == "censoring_summary", ]
  expect_equal(cens_row$statistic, 0.7, tolerance = 1e-8)
  expect_true(grepl("7 of 10 subjects \\(70\\.0%\\) were censored", cens_row$note))
})

test_that("median survival is reported as not-reached rather than a numeric NA when appropriate", {
  set.seed(3)
  n <- 100
  # heavy censoring: only ~10% of subjects ever have the event, so the KM
  # curve should never cross 0.5 survival probability
  dat <- data.frame(
    time = stats::runif(n, 1, 50),
    event = stats::rbinom(n, 1, 0.1)
  )
  fit <- fit_km_pipeline(dat, "time", "event")
  median_row <- fit$verdicts[fit$verdicts$check == "interpretation[median_survival_Overall]", ]
  if (is.na(fit$median_survival$median)) {
    expect_true(grepl("was not reached", median_row$note))
  }
})

test_that("fit_km_pipeline errors informatively on missing columns", {
  dat <- data.frame(time = 1:10, event = rep(c(0, 1), 5))
  expect_error(
    fit_km_pipeline(dat, "time", "event", strata_var = "arm"),
    "missing required column"
  )
})

test_that("fit_km_pipeline errors clearly when parametric_dist is requested but flexsurv is unavailable", {
  testthat::local_mocked_bindings(
    flexsurv_is_available = function() FALSE,
    .package = "omicsuite"
  )
  dat <- data.frame(time = 1:10, event = rep(c(0, 1), 5))
  expect_error(
    fit_km_pipeline(dat, "time", "event", parametric_dist = "weibull"),
    "requires the 'flexsurv' package|not installed"
  )
})

test_that("fit_km_pipeline falls back to a plain ggplot KM plot when survminer is unavailable", {
  testthat::local_mocked_bindings(
    survminer_is_available = function() FALSE,
    .package = "omicsuite"
  )
  set.seed(4)
  n <- 100
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7)
  )
  fit <- fit_km_pipeline(dat, "time", "event")
  expect_s3_class(fit$plots$km_plot, "ggplot")
  expect_false(inherits(fit$plots$km_plot, "ggsurvplot"))
})

test_that("parametric_dist fits a flexsurv model, adds the corresponding verdicts, and the overlay plot actually renders (no strata)", {
  # Regression test: ggplot2 builds aesthetics lazily, so merely
  # *constructing* the ggplot object (which the fit_km_pipeline() call
  # alone does) does not exercise summary.flexsurvreg()'s column-name
  # handling -- only rendering it does. This is exactly the case that broke:
  # summary.flexsurvreg() returns a single data.frame directly when there's
  # no strata, but a named list of data.frames when there is one.
  skip_if_not_installed("flexsurv")
  set.seed(5)
  n <- 150
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7)
  )
  fit <- fit_km_pipeline(dat, "time", "event", parametric_dist = "weibull")

  expect_false(is.null(fit$parametric_fit))
  expect_true(any(fit$verdicts$check == "parametric_fit[weibull]"))
  expect_true(any(fit$verdicts$check == "parametric_overlay[weibull]"))
  overlay_row <- fit$verdicts[fit$verdicts$check == "parametric_overlay[weibull]", ]
  expect_identical(overlay_row$verdict, "review")
  expect_true("parametric_overlay" %in% names(fit$plots))

  grDevices::pdf(NULL)
  expect_no_error(print(fit$plots$parametric_overlay))
  grDevices::dev.off()
})

test_that("the parametric overlay plot also renders with strata_var supplied", {
  skip_if_not_installed("flexsurv")
  set.seed(6)
  n <- 150
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
  fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm", parametric_dist = "weibull")

  grDevices::pdf(NULL)
  expect_no_error(print(fit$plots$parametric_overlay))
  grDevices::dev.off()
})

test_that("print.km_pipeline and plot.km_pipeline run without error", {
  set.seed(7)
  n <- 120
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
  fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm")
  expect_output(print(fit), "omicsuite Kaplan-Meier pipeline")
  grDevices::pdf(NULL)
  plots_shown <- plot(fit, which = "km_plot")
  grDevices::dev.off()
  expect_true("km_plot" %in% names(plots_shown))
})
