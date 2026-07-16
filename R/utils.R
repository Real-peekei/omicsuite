#' Minimal ggplot2 theme used across all omicsuite plots
#'
#' @return A ggplot2 theme object.
#' @keywords internal
#' @noRd
theme_omicsuite <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

#' Build a single verdict record
#'
#' @param check Character. Name of the diagnostic check.
#' @param passed Logical. Whether the check passed.
#' @param statistic Optional numeric test statistic.
#' @param p_value Optional numeric p-value.
#' @param note Character. Free-text interpretation.
#' @return A one-row data.frame.
#' @keywords internal
#' @noRd
make_verdict <- function(check, passed, statistic = NA_real_, p_value = NA_real_, note = "") {
  data.frame(
    check = check,
    verdict = ifelse(passed, "pass", "flagged"),
    statistic = statistic,
    p_value = p_value,
    note = note,
    stringsAsFactors = FALSE
  )
}

#' Assert that required columns exist in a data.frame
#'
#' @param data A data.frame.
#' @param cols Character vector of required column names.
#' @param data_name Character. Name to use in the error message.
#' @return Invisibly TRUE if all columns are present; otherwise throws an error.
#' @keywords internal
#' @noRd
assert_columns <- function(data, cols, data_name = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "%s is missing required column(s): %s",
      data_name, paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }
  invisible(TRUE)
}

#' Print a verdict table in a readable form
#'
#' @param verdicts A data.frame produced by `make_verdict()` rows.
#' @return Invisibly returns `verdicts`.
#' @keywords internal
#' @noRd
print_verdicts <- function(verdicts) {
  for (i in seq_len(nrow(verdicts))) {
    v <- verdicts[i, ]
    flag <- if (v$verdict == "pass") "OK  " else "FLAG"
    cat(sprintf("[%s] %-28s %s\n", flag, v$check, v$note))
  }
  invisible(verdicts)
}

#' Check whether the brms package is installed
#'
#' Small internal wrapper around [requireNamespace()] so it exists as a real
#' binding inside the omicsuite namespace and can be mocked in tests --
#' `requireNamespace()` itself, called unqualified, is only inherited via
#' lexical scoping from base, so testthat's mocking can't intercept it
#' directly.
#'
#' @return Logical.
#' @keywords internal
#' @noRd
brms_is_available <- function() {
  requireNamespace("brms", quietly = TRUE)
}
