test_that("fit_coxph_pipeline returns the expected structure (unadjusted only)", {
  set.seed(42)
  n <- 150
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )

  fit <- fit_coxph_pipeline(
    data = dat, time_var = "time", event_var = "event",
    covariates = "arm"
  )

  expect_s3_class(fit, "coxph_pipeline")
  expect_s3_class(fit$model_unadjusted, "coxph")
  expect_null(fit$model_adjusted)
  expect_identical(fit$primary_model, fit$model_unadjusted)
  expect_true("ph_test" %in% names(fit))
  expect_true(is.data.frame(fit$verdicts))
  expect_true(all(c("check", "verdict", "note") %in% names(fit$verdicts)))
  expect_true("ph_plot" %in% names(fit$plots))
  expect_true("influence_plot" %in% names(fit$plots))
})

test_that("fit_coxph_pipeline builds both unadjusted and adjusted models when adjust_for is supplied", {
  set.seed(1)
  n <- 200
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.04),
    event = stats::rbinom(n, 1, 0.8),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE)),
    age = stats::rnorm(n, 55, 10)
  )

  fit <- fit_coxph_pipeline(
    data = dat, time_var = "time", event_var = "event",
    covariates = "arm", adjust_for = "age"
  )

  expect_s3_class(fit$model_adjusted, "coxph")
  expect_identical(fit$primary_model, fit$model_adjusted)
  expect_true(length(stats::coef(fit$model_adjusted)) == 2)
  expect_true(any(grepl("^functional_form\\[age\\]$", fit$verdicts$check)))
  expect_true("functional_form_age" %in% names(fit$plots))
})

test_that("fit_coxph_pipeline errors informatively on missing columns", {
  dat <- data.frame(time = 1:10, event = rep(c(0, 1), 5))
  expect_error(
    fit_coxph_pipeline(dat, "time", "event", covariates = "arm"),
    "missing required column"
  )
})

test_that("influence diagnostics flag rows consistently with the reported cutoff", {
  set.seed(7)
  n <- 100
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.06),
    event = stats::rbinom(n, 1, 0.75),
    x1 = stats::rnorm(n)
  )
  fit <- fit_coxph_pipeline(dat, "time", "event", covariates = "x1")
  n_flagged_from_matrix <- length(unique(which(
    abs(fit$influence$dfbetas) > fit$influence$cutoff, arr.ind = TRUE
  )[, 1]))
  expect_identical(n_flagged_from_matrix, length(fit$influence$flagged_rows))
})

test_that("print.coxph_pipeline and plot.coxph_pipeline run without error", {
  set.seed(3)
  n <- 120
  dat <- data.frame(
    time = stats::rexp(n, rate = 0.05),
    event = stats::rbinom(n, 1, 0.7),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
  fit <- fit_coxph_pipeline(dat, "time", "event", covariates = "arm")
  expect_output(print(fit), "omicsuite Cox PH pipeline")
  expect_silent(invisible(grDevices::pdf(NULL)))
  plots_shown <- plot(fit, which = "ph_plot")
  grDevices::dev.off()
  expect_true("ph_plot" %in% names(plots_shown))
})
