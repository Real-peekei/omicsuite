# These tests fit a real (if tiny) Stan model via brms, so they're slow
# relative to the rest of the suite. They're gated on brms being installed
# and skipped on CRAN -- brms/Stan compilation time is exactly why it's a
# Suggests, not an Imports, dependency of omicsuite. A single fit is reused
# across expectations rather than re-fitting per test_that() block, since
# each fit pays Stan's compilation cost.

skip_if_not_installed("brms")
skip_on_cran()

make_toy_rnaseq_fit <- function() {
  set.seed(1)
  n_genes <- 6
  n_samples <- 10
  dat <- expand.grid(
    gene = paste0("gene_", seq_len(n_genes)),
    sample = paste0("sample_", seq_len(n_samples)),
    stringsAsFactors = FALSE
  )
  dat$condition <- factor(rep(
    rep(c("control", "treated"), each = n_samples / 2), n_genes
  ))
  dat$count <- stats::rnbinom(nrow(dat), mu = 40, size = 4)

  suppressWarnings(fit_rnaseq_nb_pipeline(
    data = dat, gene_var = "gene", count_var = "count",
    condition_var = "condition", sample_var = "sample",
    chains = 1, iter = 60, warmup = 30, cores = 1, seed = 1,
    refresh = 0
  ))
}

test_that("fit_rnaseq_nb_pipeline returns the expected structure", {
  fit <- make_toy_rnaseq_fit()

  expect_s3_class(fit, "rnaseq_nb_pipeline")
  expect_s3_class(fit$model, "brmsfit")
  expect_true(all(c("max_rhat", "min_neff_ratio") %in% names(fit$convergence)))
  expect_true(is.data.frame(fit$shrinkage))
  expect_true(all(c("gene", "pooled_estimate", "population_estimate") %in% names(fit$shrinkage)))
  expect_identical(nrow(fit$shrinkage), 6L)
  expect_true(is.data.frame(fit$verdicts))
  expect_true(all(c("mcmc_convergence", "dispersion[shape]", "posterior_predictive_check") %in%
    fit$verdicts$check))
  expect_true("shrinkage_plot" %in% names(fit$plots))
  expect_true("rhat_plot" %in% names(fit$plots))
})

test_that("fit_rnaseq_nb_pipeline errors clearly when brms is unavailable", {
  testthat::local_mocked_bindings(
    brms_is_available = function() FALSE,
    .package = "omicsuite"
  )
  expect_error(
    fit_rnaseq_nb_pipeline(
      data = data.frame(gene = "g1", sample = "s1", condition = "a", count = 1),
      gene_var = "gene", count_var = "count",
      condition_var = "condition", sample_var = "sample"
    ),
    "requires the 'brms' package"
  )
})

test_that("fit_rnaseq_nb_pipeline errors informatively on missing columns", {
  dat <- data.frame(gene = "g1", sample = "s1", count = 1)
  expect_error(
    fit_rnaseq_nb_pipeline(
      data = dat, gene_var = "gene", count_var = "count",
      condition_var = "condition", sample_var = "sample"
    ),
    "missing required column"
  )
})

test_that("print.rnaseq_nb_pipeline and plot.rnaseq_nb_pipeline run without error", {
  fit <- make_toy_rnaseq_fit()
  expect_output(suppressWarnings(print(fit)), "omicsuite RNA-seq")
  grDevices::pdf(NULL)
  plots_shown <- plot(fit, which = "shrinkage_plot")
  grDevices::dev.off()
  expect_true("shrinkage_plot" %in% names(plots_shown))
})
