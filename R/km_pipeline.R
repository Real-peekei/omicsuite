#' Fit a Kaplan-Meier Survival Pipeline
#'
#' Fits a nonparametric Kaplan-Meier survival curve (optionally stratified
#' by a grouping variable), runs a log-rank test for group differences, and
#' reports median survival time with its confidence interval per stratum --
#' the descriptive complement to [fit_coxph_pipeline()]'s modeling of
#' covariate effects. Optionally overlays a parametric fit (via
#' [flexsurv::flexsurvreg()]) for a visual and AIC-based check of whether a
#' parametric distribution (Weibull, exponential, etc.) is a reasonable
#' summary of the same data.
#'
#' `survminer` and `flexsurv` are `Suggests`, not hard dependencies (see
#' `CONTRIBUTING.md`'s dependency policy). Without `survminer` installed,
#' the KM plot falls back to a plain `ggplot2` step-function plot (no risk
#' table). `flexsurv` is only required when `parametric_dist` is supplied.
#' If the parametric overlay plot can't be built (e.g. an unexpected
#' `summary.flexsurvreg()` output shape on some `flexsurv` version), it's
#' silently omitted from `plots` with a warning rather than raising an
#' error -- `parametric_fit` and its AIC verdict are unaffected either way.
#'
#' @param data A data.frame containing the survival data.
#' @param time_var Character. Name of the time-to-event column.
#' @param event_var Character. Name of the event indicator column
#'   (1 = event occurred, 0 = censored).
#' @param strata_var Optional character. Name of a grouping variable for
#'   separate KM curves and a log-rank test. If `NULL` (the default), a
#'   single overall curve is fit with no group comparison.
#' @param parametric_dist Optional character naming a distribution accepted
#'   by [flexsurv::flexsurvreg()]'s `dist` argument (e.g. `"weibull"`,
#'   `"exponential"`, `"gompertz"`, `"gengamma"`). If supplied, fits a
#'   parametric model of that family and overlays it on the KM curve.
#' @param conf_type Confidence interval type passed to
#'   [survival::survfit()]'s `conf.type`. Default `"log-log"` (bounded to
#'   `[0, 1]`, generally preferred over the plain `"plain"` linear CI).
#' @param alpha Significance threshold for the log-rank test verdict.
#'   Default `0.05`.
#'
#' @return An object of class `"km_pipeline"`, a list with elements:
#' \describe{
#'   \item{fit}{The `survfit` object.}
#'   \item{parametric_fit}{The `flexsurvreg` object, or `NULL` if
#'     `parametric_dist` was not supplied.}
#'   \item{median_survival}{A data.frame: one row per stratum (or one row
#'     named `"Overall"` if `strata_var` was not supplied), with median
#'     survival time and its confidence interval.}
#'   \item{logrank_test}{The `survdiff` object, or `NULL` if `strata_var`
#'     was not supplied.}
#'   \item{plots}{A named list of `ggplot`-family objects: `km_plot`
#'     (a `ggsurvplot` object if `survminer` is installed, otherwise a plain
#'     `ggplot`), and `parametric_overlay` if `parametric_dist` was supplied
#'     and the overlay could be built (see Details).}
#'   \item{verdicts}{A data.frame: median survival per stratum and the
#'     log-rank test as `"interpretation"` rows, censoring summary and
#'     parametric-fit AIC as `"info"` rows, and the parametric overlay as a
#'     `"review"` row.}
#' }
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(
#'   time = rexp(n, rate = 0.05),
#'   event = rbinom(n, 1, 0.7),
#'   arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
#' )
#' fit <- fit_km_pipeline(dat, "time", "event", strata_var = "arm")
#' print(fit)
#'
#' @export
fit_km_pipeline <- function(data,
                             time_var,
                             event_var,
                             strata_var = NULL,
                             parametric_dist = NULL,
                             conf_type = "log-log",
                             alpha = 0.05) {

  assert_columns(data, c(time_var, event_var, strata_var), "data")

  if (!is.null(parametric_dist) && !flexsurv_is_available()) {
    stop(
      "`parametric_dist` was supplied but the 'flexsurv' package is not installed. ",
      "Install it with install.packages(\"flexsurv\") and try again, or omit `parametric_dist`.",
      call. = FALSE
    )
  }

  build_formula <- function(rhs) {
    stats::as.formula(sprintf("survival::Surv(%s, %s) ~ %s", time_var, event_var, rhs))
  }
  rhs <- if (is.null(strata_var)) "1" else strata_var
  km_formula <- build_formula(rhs)

  fit <- survival::survfit(km_formula, data = data, conf.type = conf_type)
  # survfit() (like most modeling functions) captures its call unevaluated.
  # Fit literally at the top level (survfit(Surv(t, e) ~ arm, data = dat)),
  # fit$call$formula already holds a real formula. Built and passed through
  # a variable inside a function, as here, fit$call$formula instead holds
  # the bare symbol `km_formula` -- survminer::ggsurvplot() (via
  # surv_pvalue() -> .extract.survfit()) does fit$call$formula %>%
  # as.formula(), which fails on a symbol with "object of type 'symbol' is
  # not subsettable". Patch the call so it holds the actual formula object.
  fit$call$formula <- km_formula

  # --- median survival table -----------------------------------------------------
  surv_table <- summary(fit)$table
  if (is.null(dim(surv_table))) {
    surv_table <- matrix(surv_table, nrow = 1, dimnames = list("Overall", names(surv_table)))
  }
  strata_names <- if (is.null(strata_var)) {
    "Overall"
  } else {
    sub(sprintf("^%s=", strata_var), "", rownames(surv_table))
  }
  median_survival <- data.frame(
    stratum = strata_names,
    records = surv_table[, "records"],
    events = surv_table[, "events"],
    median = surv_table[, "median"],
    lower = surv_table[, "0.95LCL"],
    upper = surv_table[, "0.95UCL"],
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  median_verdicts <- do.call(rbind, lapply(seq_len(nrow(median_survival)), function(i) {
    row <- median_survival[i, ]
    note <- if (is.na(row$median)) {
      sprintf(
        "Median survival was not reached for `%s` (fewer than half of the %d subjects at risk had an event) -- report the follow-up period explicitly rather than a median.",
        row$stratum, row$records
      )
    } else {
      sprintf(
        "Median survival for `%s`: %.2f (95%% CI [%s, %s]) among %d subjects, %d events.",
        row$stratum, row$median,
        if (is.na(row$lower)) "NA" else sprintf("%.2f", row$lower),
        if (is.na(row$upper)) "NA" else sprintf("%.2f", row$upper),
        row$records, row$events
      )
    }
    make_note(
      check = sprintf("interpretation[median_survival_%s]", row$stratum),
      label = "interpretation",
      statistic = row$median, p_value = NA_real_,
      note = note
    )
  }))

  # --- log-rank test ---------------------------------------------------------------
  logrank_test <- NULL
  logrank_verdict <- NULL
  if (!is.null(strata_var) && length(unique(data[[strata_var]])) > 1) {
    logrank_test <- survival::survdiff(km_formula, data = data)
    logrank_df <- length(logrank_test$n) - 1
    logrank_p <- stats::pchisq(logrank_test$chisq, df = logrank_df, lower.tail = FALSE)
    logrank_verdict <- make_note(
      check = "interpretation[log_rank_test]",
      label = "interpretation",
      statistic = logrank_test$chisq, p_value = logrank_p,
      note = sprintf(
        "Log-rank test across `%s`: chi-squared = %.2f (df = %d), p = %.4f. %s",
        strata_var, logrank_test$chisq, logrank_df, logrank_p,
        if (logrank_p < alpha) {
          "Survival differs significantly across groups -- this is a global test, not a pairwise one, so a significant result with more than two groups doesn't say which pair differs."
        } else {
          "No significant evidence that survival differs across groups at this alpha."
        }
      )
    )
  }

  # --- censoring summary ------------------------------------------------------------
  n_total <- nrow(data)
  n_events <- sum(data[[event_var]] == 1, na.rm = TRUE)
  censoring_verdict <- make_note(
    check = "censoring_summary",
    label = "info",
    statistic = 1 - (n_events / n_total), p_value = NA_real_,
    note = sprintf(
      "%d of %d subjects (%.1f%%) were censored rather than observed to have the event.",
      n_total - n_events, n_total, 100 * (1 - n_events / n_total)
    )
  )

  # --- optional parametric overlay --------------------------------------------------
  parametric_fit <- NULL
  parametric_verdicts <- NULL
  parametric_overlay_error <- NULL
  if (!is.null(parametric_dist)) {
    parametric_fit <- flexsurv::flexsurvreg(km_formula, data = data, dist = parametric_dist)
    aic_verdict <- make_note(
      check = sprintf("parametric_fit[%s]", parametric_dist),
      label = "info",
      statistic = parametric_fit$AIC, p_value = NA_real_,
      note = sprintf(
        "%s parametric fit: AIC = %.1f, log-likelihood = %.1f. Compare AIC across candidate distributions if trying to choose one; lower is better.",
        tools::toTitleCase(parametric_dist), parametric_fit$AIC, parametric_fit$loglik
      )
    )
    overlay_verdict <- make_note(
      check = sprintf("parametric_overlay[%s]", parametric_dist),
      label = "review",
      statistic = NA_real_, p_value = NA_real_,
      note = "Compare the parametric curve against the KM step function in plots$parametric_overlay -- systematic divergence suggests this distributional family is a poor fit, regardless of what the AIC alone suggests."
    )
    parametric_verdicts <- rbind(aic_verdict, overlay_verdict)
  }

  verdicts <- rbind(median_verdicts, logrank_verdict, censoring_verdict, parametric_verdicts)
  rownames(verdicts) <- NULL

  # --- plots -------------------------------------------------------------------
  plots <- list()

  if (survminer_is_available()) {
    km_plot <- survminer::ggsurvplot(
      fit, data = data, conf.int = TRUE, risk.table = TRUE,
      legend.title = if (is.null(strata_var)) "" else strata_var
    )
    km_plot$plot <- km_plot$plot + theme_omicsuite()
    plots$km_plot <- km_plot
  } else {
    # Fallback with no survminer installed: a plain ggplot2 step-function
    # plot. Uses linear (not stair-stepped) ribbon interpolation between KM
    # jump points as a documented simplification -- survminer's stepped
    # ribbon construction is more visually correct but not worth
    # reimplementing here just for the no-survminer fallback path.
    n_pts <- length(fit$time)
    strata_labels <- if (!is.null(fit$strata)) {
      rep(sub(sprintf("^%s=", strata_var), "", names(fit$strata)), fit$strata)
    } else {
      rep("Overall", n_pts)
    }
    km_df <- data.frame(
      time = fit$time, surv = fit$surv, lower = fit$lower, upper = fit$upper,
      stratum = strata_labels
    )
    plots$km_plot <- ggplot2::ggplot(
      km_df, ggplot2::aes(x = .data$time, y = .data$surv, color = .data$stratum, fill = .data$stratum)
    ) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower, ymax = .data$upper), alpha = 0.15, color = NA) +
      ggplot2::geom_step(linewidth = 0.8) +
      ggplot2::labs(
        title = "Kaplan-Meier survival curve",
        subtitle = "Install 'survminer' for a risk table and stair-stepped confidence bands",
        x = "Time", y = "Survival probability", color = strata_var, fill = strata_var
      ) +
      theme_omicsuite()
  }

  if (!is.null(parametric_fit)) {
    # Root cause, confirmed by inspecting str(summary(flexsurv_fit, ...))
    # directly: summary.flexsurvreg(..., ci = FALSE) returns an *unnamed*
    # list of one data.frame per stratum (class "summary.flexsurvreg"),
    # each data.frame with columns "time" and "est" -- not a bare
    # data.frame (the first guess) and not a *named* list (the second). An
    # unnamed list means names(parametric_pred) is NULL, so the original
    # lapply(names(parametric_pred), ...) iterated zero times, silently
    # producing an empty parametric_df -- which is what surfaced downstream
    # as "column 'time' not found" during ggplot rendering. The unnamed-list
    # branch below assigns placeholder names before iterating, which fixes
    # this. The tryCatch stays as insurance against other flexsurv versions
    # behaving differently still -- on any failure, the overlay plot is
    # omitted (with a warning carrying the actual structure seen) instead
    # of raising an error that would halt package checks or vignette builds.
    # parametric_fit and its AIC verdict are unaffected either way.
    overlay_plot <- tryCatch({
      grid_times <- seq(0, max(data[[time_var]], na.rm = TRUE), length.out = 200)
      parametric_pred <- summary(parametric_fit, t = grid_times, ci = FALSE)

      # Normalize to a named list of data.frames regardless of whether
      # summary.flexsurvreg() returned one directly (no covariates/strata)
      # or a list of them (one per stratum).
      if (is.data.frame(parametric_pred)) {
        parametric_pred <- stats::setNames(list(parametric_pred), "Overall")
      }
      if (!is.list(parametric_pred)) {
        stop(sprintf(
          "summary.flexsurvreg() returned an object of class %s, not a data.frame or list -- can't build the overlay plot from this.",
          paste(class(parametric_pred), collapse = "/")
        ))
      }
      if (is.null(names(parametric_pred)) || any(names(parametric_pred) == "")) {
        names(parametric_pred) <- paste0("stratum_", seq_along(parametric_pred))
      }

      parametric_df <- do.call(rbind, lapply(names(parametric_pred), function(nm) {
        d <- parametric_pred[[nm]]
        if (!is.data.frame(d) || !all(c("time", "est") %in% names(d))) {
          stop(sprintf(
            "element '%s' of summary.flexsurvreg()'s output is not a data.frame with 'time'/'est' columns (class: %s; names: %s) -- flexsurv's output shape may differ from what this code assumes on your installed version.",
            nm, paste(class(d), collapse = "/"),
            if (is.null(names(d))) "(none)" else paste(names(d), collapse = ", ")
          ))
        }
        d$stratum <- if (is.null(strata_var)) "Overall" else sub(sprintf("^%s=", strata_var), "", nm)
        d
      }))

      n_pts <- length(fit$time)
      strata_labels <- if (!is.null(fit$strata)) {
        rep(sub(sprintf("^%s=", strata_var), "", names(fit$strata)), fit$strata)
      } else {
        rep("Overall", n_pts)
      }
      km_df <- data.frame(time = fit$time, surv = fit$surv, stratum = strata_labels)

      ggplot2::ggplot() +
        ggplot2::geom_step(
          data = km_df,
          ggplot2::aes(x = .data$time, y = .data$surv, color = .data$stratum),
          linewidth = 0.8
        ) +
        ggplot2::geom_line(
          data = parametric_df,
          ggplot2::aes(x = .data$time, y = .data$est, color = .data$stratum),
          linetype = "dashed", linewidth = 0.8
        ) +
        ggplot2::labs(
          title = sprintf("KM curve (solid) vs. %s fit (dashed)", parametric_dist),
          x = "Time", y = "Survival probability", color = strata_var
        ) +
        theme_omicsuite()
    }, error = function(e) {
      parametric_overlay_error <<- conditionMessage(e)
      NULL
    })

    if (!is.null(overlay_plot)) {
      plots$parametric_overlay <- overlay_plot
    } else {
      warning(
        "Could not build the parametric overlay plot (fit_km_pipeline()): ", parametric_overlay_error,
        " The parametric fit itself (result$parametric_fit) and its AIC verdict are unaffected.",
        call. = FALSE
      )
    }
  }

  structure(
    list(
      call = match.call(),
      fit = fit,
      parametric_fit = parametric_fit,
      median_survival = median_survival,
      logrank_test = logrank_test,
      plots = plots,
      verdicts = verdicts
    ),
    class = "km_pipeline"
  )
}

#' Print a `km_pipeline` object
#'
#' @param x A `km_pipeline` object, as returned by [fit_km_pipeline()].
#' @param ... Ignored.
#' @return Invisibly returns `x`.
#' @export
print.km_pipeline <- function(x, ...) {
  cat("<omicsuite Kaplan-Meier pipeline>\n\n")
  print(x$fit)
  cat("\nDiagnostic verdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' Summarize a `km_pipeline` object
#'
#' @param object A `km_pipeline` object, as returned by [fit_km_pipeline()].
#' @param ... Ignored.
#' @return A list with elements `median_survival`, `logrank_test`, and
#'   `verdicts`.
#' @export
summary.km_pipeline <- function(object, ...) {
  list(
    median_survival = object$median_survival,
    logrank_test = object$logrank_test,
    verdicts = object$verdicts
  )
}

#' Plot a `km_pipeline` object
#'
#' @param x A `km_pipeline` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.km_pipeline <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    if (!is.null(x$plots[[nm]])) print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
