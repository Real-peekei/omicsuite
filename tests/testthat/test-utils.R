test_that("make_verdict() defaults plot to NA and accepts an explicit value", {
  v_default <- make_verdict(check = "some_check", passed = TRUE)
  expect_true("plot" %in% names(v_default))
  expect_true(is.na(v_default$plot))

  v_linked <- make_verdict(check = "some_check", passed = TRUE, plot = "some_plot")
  expect_identical(v_linked$plot, "some_plot")
})

test_that("make_note() defaults plot to NA and accepts an explicit value", {
  n_default <- make_note(check = "some_check", label = "info")
  expect_true("plot" %in% names(n_default))
  expect_true(is.na(n_default$plot))

  n_linked <- make_note(check = "some_check", label = "interpretation", plot = "some_plot")
  expect_identical(n_linked$plot, "some_plot")
})

test_that("rows from make_verdict() and make_note() rbind cleanly with a mix of linked and unlinked plots", {
  verdicts <- rbind(
    make_verdict(check = "a", passed = TRUE, plot = "plot_a"),
    make_note(check = "b", label = "info"),
    make_note(check = "c", label = "interpretation", plot = "plot_c")
  )
  expect_identical(nrow(verdicts), 3L)
  expect_identical(verdicts$plot, c("plot_a", NA_character_, "plot_c"))
})

test_that("print_verdicts appends a plot pointer only for rows that have one", {
  verdicts <- rbind(
    make_verdict(check = "a", passed = TRUE, note = "note a", plot = "plot_a"),
    make_note(check = "b", label = "info", note = "note b")
  )
  out <- testthat::capture_output(print_verdicts(verdicts))
  expect_true(grepl("\\[see plots\\$plot_a\\]", out))
  expect_false(grepl("\\[see plots\\$NA\\]", out))
})
