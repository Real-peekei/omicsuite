#' Fit a Competing Risks Pipeline
#'
#' Fits nonparametric cumulative incidence functions (CIFs) for each event
#' type, optionally compares them across groups with Gray's test, and fits
#' a Fine-Gray subdistribution hazard model for the event type of primary
#' interest via [cmprsk::crr()].
#'
#' **`event_var` means something different here than in [fit_km_pipeline()]
#' or [fit_coxph_pipeline()].** Those treat `event_var` as a binary
#' indicator (1 = event, 0 = censored). Here it must encode *which* event
#' type occurred: `0` = censored (or whatever `cencode` is set to), `1` =
#' the event of primary interest, `2`, `3`, ... = competing event types
#' that preclude the event of interest from ever happening (e.g. death from
#' another cause when the event of interest is relapse). Reusing a plain
#' 0/1 binary indicator here silently treats every non-event as a
#' competing risk of type... there is none, which just reduces to the
#' single-event case -- not wrong, but usually not what you meant to ask.
#'
#' The reason a dedicated model matters: fitting a standard Cox model that
#' treats competing events as ordinary censoring overestimates the event of
#' interest's cumulative incidence, because it implicitly assumes subjects
#' who experienced a competing event would have gone on to experience the
#' event of interest had they lived. The Fine-Gray model avoids that
#' assumption by keeping competing-event subjects in the risk set on the
#' subdistribution timescale instead of removing them.
#'
#' @param data A data.frame containing the survival data.
#' @param time_var Character. Name of the time-to-event column.
#' @param event_var Character. Name of the event-type column (see Details
#'   above for its encoding -- not the same as `event_var` elsewhere in
#'   `omicsuite`).
#' @param failcode The event code (as it appears in `data[[event_var]]`)
#'   representing the event of primary interest. Default `1`.
#' @param cencode The event code representing censoring. Default `0`.
#' @param group_var Optional character. Grouping variable for separate CIF
#'   curves and Gray's test. If `NULL` (the default), a single CIF per
#'   event type is estimated with no group comparison.
#' @param covariates Optional character vector of covariate names for the
#'   Fine-Gray model (`cov1` in [cmprsk::crr()]). Factors are dummy-coded
#'   automatically via [stats::model.matrix()]. If `NULL` (the default),
#'   only the CIF/Gray's-test part of the pipeline runs -- no subdistribution
#'   hazard model is fit.
#' @param alpha Significance threshold for Gray's test and the Fine-Gray
#'   model's verdicts. Default `0.05`.
#'
#' @return An object of class `"competing_risks_pipeline"`, a list with
#'   elements:
#' \describe{
#'   \item{cuminc_fit}{The `cuminc` object from [cmprsk::cuminc()].}
#'   \item{crr_fit}{The `crr` object from [cmprsk::crr()], or `NULL` if
#'     `covariates` was not supplied.}
#'   \item{cif_summary}{A long-format data.frame: time, cumulative
#'     incidence estimate, event type, and group (if `group_var` was
#'     supplied) -- the tidy form of `cuminc_fit` used for plotting.}
#'   \item{event_summary}{A data.frame: counts and proportions for each
#'     event code observed in `data[[event_var]]`, including censoring.}
#'   \item{plots}{A named list of `ggplot` objects: `cif_plot`.}
#'   \item{verdicts}{A data.frame: event-type breakdown as `"info"`, Gray's
#'     test as `"interpretation"` (when `group_var` is supplied), and
#'     subdistribution hazard ratios as `"interpretation"` per covariate
#'     term (when `covariates` is supplied). A `plot` column names which
#'     entry in `plots` each row explains.}
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("cmprsk", quietly = TRUE)) {
#'   set.seed(1)
#'   n <- 300
#'   dat <- data.frame(
#'     time = rexp(n, rate = 0.1),
#'     # 0 = censored, 1 = event of interest, 2 = competing event
#'     event = sample(c(0, 1, 2), n, replace = TRUE, prob = c(0.3, 0.5, 0.2)),
#'     arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
#'   )
#'   fit <- fit_competing_risks_pipeline(
#'     dat, "time", "event", group_var = "arm", covariates = "arm"
#'   )
#'   print(fit)
#' }
#' }
#'
#' @export
fit_competing_risks_pipeline <- function(data,
                                          time_var,
                                          event_var,
                                          failcode = 1,
                                          cencode = 0,
                                          group_var = NULL,
                                          covariates = NULL,
                                          alpha = 0.05) {

  assert_columns(data, c(time_var, event_var, group_var, covariates), "data")

  if (!cmprsk_is_available()) {
    stop(
      "fit_competing_risks_pipeline() requires the 'cmprsk' package, which is not installed. ",
      "Install it with install.packages(\"cmprsk\") and try again.",
      call. = FALSE
    )
  }

  event_codes <- data[[event_var]]
  if (!(failcode %in% event_codes)) {
    stop(sprintf(
      "`failcode = %s` does not appear in data[[event_var]]. Observed codes: %s.",
      failcode, paste(sort(unique(event_codes)), collapse = ", ")
    ), call. = FALSE)
  }

  # --- event-type breakdown ---------------------------------------------------------
  event_counts <- table(event_codes)
  event_summary <- data.frame(
    code = names(event_counts),
    n = as.integer(event_counts),
    proportion = as.numeric(event_counts) / sum(event_counts),
    stringsAsFactors = FALSE
  )
  event_summary$label <- ifelse(
    event_summary$code == as.character(cencode), "censored",
    ifelse(event_summary$code == as.character(failcode), "event of interest", "competing event")
  )
  event_verdict <- make_note(
    check = "event_type_breakdown",
    label = "info",
    statistic = event_summary$proportion[event_summary$code == as.character(failcode)],
    p_value = NA_real_,
    note = paste(
      sprintf("code %s (%s): n = %d (%.1f%%)", event_summary$code, event_summary$label,
              event_summary$n, 100 * event_summary$proportion),
      collapse = "; "
    )
  )

  # --- cumulative incidence functions + Gray's test ----------------------------------
  # cuminc() checks missing(group) internally to decide whether to use its
  # own default (rep(1, length(ftime))). Passing group = NULL explicitly
  # (rather than omitting the argument) makes missing(group) FALSE, so
  # cuminc() tries as.factor(NULL) -- a length-0 factor -- against the full
  # ftime/fstatus vectors, which is exactly the "differing number of rows"
  # error this produced. Build the argument list dynamically instead, so
  # `group` is only present in the call at all when group_var is supplied.
  cuminc_args <- list(ftime = data[[time_var]], fstatus = event_codes, cencode = cencode)
  if (!is.null(group_var)) {
    cuminc_args$group <- data[[group_var]]
  }
  cuminc_fit <- do.call(cmprsk::cuminc, cuminc_args)

  # cuminc() returns a list: one element per group/event-type combination
  # (list names like "<group> <event code>", or just "<event code>" with no
  # group), each a list with $time/$est/$var, PLUS a $Tests element (a
  # matrix: rows = event types, columns include "stat"/"pv"/"df") when
  # group has more than one level. Split those apart defensively rather
  # than assuming which elements are curves vs. the test.
  is_curve <- vapply(cuminc_fit, function(x) is.list(x) && all(c("time", "est") %in% names(x)), logical(1))
  curve_elements <- cuminc_fit[is_curve]

  cif_summary <- do.call(rbind, lapply(names(curve_elements), function(nm) {
    curve <- curve_elements[[nm]]
    parts <- strsplit(trimws(nm), "\\s+")[[1]]
    event_code <- utils::tail(parts, 1)
    group_label <- if (is.null(group_var)) "Overall" else paste(utils::head(parts, -1), collapse = " ")
    data.frame(
      time = curve$time, estimate = curve$est,
      event_code = event_code, group = group_label,
      stringsAsFactors = FALSE
    )
  }))

  gray_verdict <- NULL
  if (!is.null(group_var) && "Tests" %in% names(cuminc_fit)) {
    gray_table <- cuminc_fit$Tests
    gray_verdict <- do.call(rbind, lapply(seq_len(nrow(gray_table)), function(i) {
      row <- gray_table[i, ]
      event_code <- rownames(gray_table)[i]
      make_note(
        check = sprintf("interpretation[gray_test_event_%s]", event_code),
        label = "interpretation",
        statistic = row["stat"], p_value = row["pv"],
        note = sprintf(
          "Gray's test for event code %s across `%s`: statistic = %.2f, p = %.4f. %s",
          event_code, group_var, row["stat"], row["pv"],
          if (row["pv"] < alpha) {
            "Cumulative incidence differs significantly across groups for this event type."
          } else {
            "No significant evidence that cumulative incidence differs across groups for this event type at this alpha."
          }
        ),
        plot = "cif_plot"
      )
    }))
  }

  # --- Fine-Gray subdistribution hazard model (optional) -----------------------------
  crr_fit <- NULL
  crr_verdicts <- NULL
  if (!is.null(covariates)) {
    cov_matrix <- stats::model.matrix(stats::reformulate(covariates), data = data)[, -1, drop = FALSE]
    crr_fit <- cmprsk::crr(
      ftime = data[[time_var]], fstatus = event_codes, cov1 = cov_matrix,
      failcode = failcode, cencode = cencode
    )
    crr_summary <- summary(crr_fit)
    # summary.crr() returns a list including $coef, a matrix with columns
    # "coef", "exp(coef)", "se(coef)", "z", "p-value" (one row per
    # covariate term) -- following the same coef/se/p-value table shape
    # used throughout the survival ecosystem (coxph, etc.), not verified by
    # execution here.
    coef_table <- crr_summary$coef
    crr_verdicts <- do.call(rbind, lapply(rownames(coef_table), function(term) {
      shr <- coef_table[term, "exp(coef)"]
      p_val <- coef_table[term, "p-value"]
      pct_change <- (shr - 1) * 100
      direction <- if (shr < 1) "lower" else "higher"
      sig <- if (!is.na(p_val) && p_val < alpha) "statistically significant" else "not statistically significant at this alpha"
      make_note(
        check = sprintf("interpretation[subdistribution_hr_%s]", term),
        label = "interpretation",
        statistic = shr, p_value = p_val,
        note = sprintf(
          "Subdistribution HR for `%s` = %.3f, p = %.4f. Associated with a %.1f%% %s subdistribution hazard of the event of interest (%s). This is not the same quantity as a cause-specific hazard ratio from a standard Cox model -- see this function's documentation for why.",
          term, shr, p_val, abs(pct_change), direction, sig
        ),
        plot = "cif_plot"
      )
    }))
  }

  conceptual_note <- make_note(
    check = "subdistribution_vs_cause_specific",
    label = "info",
    statistic = NA_real_, p_value = NA_real_,
    note = "This pipeline models the subdistribution hazard (Fine-Gray), which is the right quantity for predicting cumulative incidence in the presence of competing risks. It answers a different question than a cause-specific hazard model (an ordinary Cox model with competing events treated as censored) -- the two can even point in different directions for the same covariate. Report which one you used and why.",
    plot = "cif_plot"
  )

  verdicts <- rbind(event_verdict, conceptual_note, gray_verdict, crr_verdicts)
  rownames(verdicts) <- NULL

  # --- plots -------------------------------------------------------------------
  cif_summary$event_code <- factor(cif_summary$event_code)
  plots <- list(
    cif_plot = ggplot2::ggplot(
      cif_summary,
      ggplot2::aes(x = .data$time, y = .data$estimate, color = .data$group, linetype = .data$event_code)
    ) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::labs(
        title = "Cumulative incidence by event type",
        subtitle = if (is.null(group_var)) NULL else sprintf("Grouped by %s", group_var),
        x = "Time", y = "Cumulative incidence",
        color = if (is.null(group_var)) NULL else group_var, linetype = "Event code"
      ) +
      theme_omicsuite()
  )

  structure(
    list(
      call = match.call(),
      cuminc_fit = cuminc_fit,
      crr_fit = crr_fit,
      cif_summary = cif_summary,
      event_summary = event_summary,
      plots = plots,
      verdicts = verdicts
    ),
    class = "competing_risks_pipeline"
  )
}

#' Print a `competing_risks_pipeline` object
#'
#' @param x A `competing_risks_pipeline` object, as returned by
#'   [fit_competing_risks_pipeline()].
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.competing_risks_pipeline <- function(x, ...) {
  cat("<omicsuite competing risks pipeline>\n\n")
  cat("Event summary:\n")
  print(x$event_summary[, c("code", "label", "n", "proportion")])
  if (!is.null(x$crr_fit)) {
    cat("\nFine-Gray subdistribution hazard model:\n")
    print(x$crr_fit)
  }
  cat("\nDiagnostic verdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' Summarize a `competing_risks_pipeline` object
#'
#' @param object A `competing_risks_pipeline` object, as returned by
#'   [fit_competing_risks_pipeline()].
#' @param ... Ignored.
#' @return A list with elements `event_summary`, `cif_summary`, and
#'   `verdicts`.
#' @export
summary.competing_risks_pipeline <- function(object, ...) {
  list(
    event_summary = object$event_summary,
    cif_summary = object$cif_summary,
    verdicts = object$verdicts
  )
}

#' Plot a `competing_risks_pipeline` object
#'
#' @param x A `competing_risks_pipeline` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.competing_risks_pipeline <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    if (!is.null(x$plots[[nm]])) print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
