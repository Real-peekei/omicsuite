skip_if_not_installed("cmprsk")

make_toy_cr_data <- function(n = 300, seed = 1) {
  set.seed(seed)
  data.frame(
    time = stats::rexp(n, rate = 0.1),
    # 0 = censored, 1 = event of interest, 2 = competing event
    event = sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3, 0.5, 0.2)),
    arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
  )
}

test_that("fit_competing_risks_pipeline returns the expected structure (no group, no covariates)", {
  dat <- make_toy_cr_data()
  fit <- fit_competing_risks_pipeline(dat, "time", "event")

  expect_s3_class(fit, "competing_risks_pipeline")
  expect_null(fit$crr_fit)
  expect_true(is.data.frame(fit$event_summary))
  expect_true(is.data.frame(fit$cif_summary))
  expect_setequal(fit$cif_summary$group, "Overall")
  expect_setequal(fit$event_summary$label, c("censored", "event of interest", "competing event"))
  expect_true(any(fit$verdicts$check == "event_type_breakdown"))
  expect_true(any(fit$verdicts$check == "subdistribution_vs_cause_specific"))
  expect_true("cif_plot" %in% names(fit$plots))
})

test_that("fit_competing_risks_pipeline with group_var produces Gray's test verdicts", {
  dat <- make_toy_cr_data(seed = 2)
  fit <- fit_competing_risks_pipeline(dat, "time", "event", group_var = "arm")

  expect_true(any(grepl("^interpretation\\[gray_test_event_", fit$verdicts$check)))
  gray_rows <- fit$verdicts[grepl("^interpretation\\[gray_test_event_", fit$verdicts$check), ]
  expect_true(all(!is.na(gray_rows$p_value)))
  expect_setequal(unique(fit$cif_summary$group), c("control", "treatment"))
})

test_that("fit_competing_risks_pipeline with covariates fits a Fine-Gray model and reports subdistribution HRs", {
  dat <- make_toy_cr_data(seed = 3)
  fit <- fit_competing_risks_pipeline(dat, "time", "event", covariates = "arm")

  expect_false(is.null(fit$crr_fit))
  expect_true(any(grepl("^interpretation\\[subdistribution_hr_", fit$verdicts$check)))
  hr_row <- fit$verdicts[grepl("^interpretation\\[subdistribution_hr_", fit$verdicts$check), ]
  expect_identical(nrow(hr_row), 1L)
  expect_true(hr_row$statistic > 0)  # subdistribution HR is always positive
  expect_true(grepl("Subdistribution HR", hr_row$note))
})

test_that("event_type_breakdown accounts for all subjects", {
  dat <- make_toy_cr_data(n = 250, seed = 4)
  fit <- fit_competing_risks_pipeline(dat, "time", "event")
  expect_equal(sum(fit$event_summary$n), 250L)
  expect_equal(sum(fit$event_summary$proportion), 1, tolerance = 1e-8)
})

test_that("fit_competing_risks_pipeline errors informatively on missing columns and bad failcode", {
  dat <- make_toy_cr_data(seed = 5)
  expect_error(
    fit_competing_risks_pipeline(dat, "time", "event", group_var = "nonexistent"),
    "missing required column"
  )
  expect_error(
    fit_competing_risks_pipeline(dat, "time", "event", failcode = 99),
    "does not appear in data"
  )
})

test_that("fit_competing_risks_pipeline errors clearly when cmprsk is unavailable", {
  testthat::local_mocked_bindings(
    cmprsk_is_available = function() FALSE,
    .package = "omicsuite"
  )
  dat <- make_toy_cr_data(n = 20, seed = 6)
  expect_error(
    fit_competing_risks_pipeline(dat, "time", "event"),
    "requires the 'cmprsk' package"
  )
})

test_that("cif_plot carries a non-empty caption combining Gray's test and subdistribution HR notes", {
  dat <- make_toy_cr_data(seed = 8)
  fit <- fit_competing_risks_pipeline(dat, "time", "event", group_var = "arm", covariates = "arm")
  caption <- fit$plots$cif_plot$labels$caption
  expect_false(is.null(caption))
  expect_true(grepl("Gray's test", caption) || grepl("gray", tolower(caption)))
  expect_true(grepl("Subdistribution HR", caption))
})

test_that("print.competing_risks_pipeline and plot.competing_risks_pipeline run without error", {
  dat <- make_toy_cr_data(seed = 7)
  fit <- fit_competing_risks_pipeline(dat, "time", "event", group_var = "arm", covariates = "arm")
  expect_output(print(fit), "omicsuite competing risks pipeline")
  grDevices::pdf(NULL)
  expect_no_error(plot(fit, which = "cif_plot"))
  grDevices::dev.off()
})
