skip_if_not_installed("RGCCA")

make_toy_blocks <- function(n = 40, seed = 1) {
  set.seed(seed)
  sample_ids <- paste0("s", seq_len(n))
  block1 <- matrix(stats::rnorm(n * 15), nrow = n, dimnames = list(sample_ids, paste0("f1_", 1:15)))
  block2 <- matrix(stats::rnorm(n * 10), nrow = n, dimnames = list(sample_ids, paste0("f2_", 1:10)))
  shared <- stats::rnorm(n)
  block1[, 1] <- block1[, 1] + shared
  block2[, 1] <- block2[, 1] + shared
  list(omics_a = block1, omics_b = block2)
}

test_that("integrate_multiomics returns the expected structure", {
  blocks <- make_toy_blocks()
  fit <- integrate_multiomics(blocks, ncomp = 2, n_boot = 10, seed = 1)

  expect_s3_class(fit, "multiomics_pipeline")
  expect_identical(names(fit$scores), c("omics_a", "omics_b"))
  expect_true(all(c("comp1", "comp2", "sample_id") %in% names(fit$scores$omics_a)))
  expect_true(is.data.frame(fit$variance_explained))
  expect_true(is.data.frame(fit$stability))
  expect_true(is.data.frame(fit$verdicts))
  expect_true(any(grepl("^sample_alignment$", fit$verdicts$check)))
  expect_true(any(grepl("^variance_explained\\[omics_a\\]$", fit$verdicts$check)))
  expect_true(any(grepl("^loading_stability\\[omics_a\\]$", fit$verdicts$check)))
  expect_true(all(c("block_scores", "variance_explained", "stability") %in% names(fit$plots)))
})

test_that("integrate_multiomics detects the injected shared signal in variance explained", {
  blocks <- make_toy_blocks()
  fit <- integrate_multiomics(blocks, ncomp = 2, n_boot = 5, seed = 1)
  comp1_ave <- fit$variance_explained$variance_explained[
    fit$variance_explained$component == 1
  ]
  # both blocks share injected signal on component 1, so AVE should be
  # comfortably above what pure noise blocks would produce
  expect_true(all(comp1_ave > 0.05))
})

test_that("integrate_multiomics aligns samples and reports dropped rows", {
  blocks <- make_toy_blocks(n = 40)
  # drop a few samples from one block only
  blocks$omics_b <- blocks$omics_b[1:35, , drop = FALSE]

  fit <- integrate_multiomics(blocks, ncomp = 2, n_boot = 5, seed = 1)
  alignment_row <- fit$verdicts[fit$verdicts$check == "sample_alignment", ]

  expect_identical(alignment_row$verdict, "flagged")
  expect_identical(alignment_row$statistic, 35)
  expect_identical(nrow(fit$scores$omics_a), 35L)
})

test_that("integrate_multiomics errors informatively on unnamed or single-block input", {
  blocks <- make_toy_blocks()
  names(blocks) <- NULL
  expect_error(integrate_multiomics(blocks), "fully named list")

  blocks2 <- make_toy_blocks()
  expect_error(integrate_multiomics(blocks2[1]), "at least two omics layers")
})

test_that("integrate_multiomics errors informatively when row names are missing", {
  blocks <- make_toy_blocks()
  rownames(blocks$omics_a) <- NULL
  expect_error(integrate_multiomics(blocks), "row names")
})

test_that("print.multiomics_pipeline and plot.multiomics_pipeline run without error", {
  blocks <- make_toy_blocks()
  fit <- integrate_multiomics(blocks, ncomp = 2, n_boot = 5, seed = 1)
  expect_output(print(fit), "omicsuite multi-omics")
  grDevices::pdf(NULL)
  plots_shown <- plot(fit, which = "variance_explained")
  grDevices::dev.off()
  expect_true("variance_explained" %in% names(plots_shown))
})
