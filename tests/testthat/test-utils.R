test_that("add_caption returns the plot unchanged when notes is empty, NA, or blank", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()

  expect_null(add_caption(p, character(0))$labels$caption)
  expect_null(add_caption(p, NA_character_)$labels$caption)
  expect_null(add_caption(p, "")$labels$caption)
  expect_null(add_caption(p, c(NA_character_, ""))$labels$caption)
})

test_that("add_caption attaches wrapped note text as a plot caption", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()

  short_note <- "HR = 0.72, p = 0.01."
  captioned <- add_caption(p, short_note)
  expect_identical(captioned$labels$caption, short_note)

  long_note <- paste(rep("word", 40), collapse = " ")
  captioned_long <- add_caption(p, long_note, width = 20)
  expect_true(grepl("\n", captioned_long$labels$caption))
  # wrapping shouldn't drop or duplicate any words
  expect_identical(
    gsub("\n", " ", captioned_long$labels$caption),
    long_note
  )
})

test_that("add_caption joins multiple notes with a blank line and drops empty entries", {
  p <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point()

  captioned <- add_caption(p, c("First note.", "", NA_character_, "Second note."))
  expect_identical(captioned$labels$caption, "First note.\n\nSecond note.")
})
