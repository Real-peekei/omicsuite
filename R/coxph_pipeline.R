#' Fit a Cox Proportional Hazards Pipeline with Full Diagnostics
#'
#' Fits a Cox proportional hazards model and runs the diagnostic suite you'd
#' want in a methods section before trusting the hazard ratios: an
#' unadjusted-versus-adjusted comparison, global and covariate-level
#' proportional hazards testing, influence diagnostics, and functional form
#' checks for continuous covariates via martingale residuals.
#'
#' @param data A data.frame containing the survival data.
#' @param time_var Character. Name of the time-to-event column.
#' @param event_var Character. Name of the event indicator column
#'   (1 = event occurred, 0 = censored).
#' @param covariates Character vector of covariate names for the primary
#'   model of interest.
#' @param adjust_for Optional character vector of additional covariates to
#'   adjust for. When supplied, both an unadjusted model (`covariates` only)
#'   and an adjusted model (`covariates` + `adjust_for`) are fit, so you can
#'   see how much the estimates of interest move.
#' @param id_var Optional character. Cluster/subject identifier, for
#'   robust sandwich variance (`cluster()`) or a shared frailty term.
#' @param frailty Logical. If `TRUE` and `id_var` is supplied, fits a gamma
#'   frailty term instead of a robust cluster-adjusted variance. Default
#'   `FALSE`.
#' @param alpha Significance threshold used to flag diagnostic checks.
#'   Default `0.05`.
#' @param influence_cutoff_sd Numeric. Observations with a scaled dfbeta
#'   beyond this many "typical" units (`2 / sqrt(n)` scaled by this factor)
#'   are flagged as influential. Default `1` (i.e. the standard `2/sqrt(n)`
#'   cutoff).
#'
#' @return An object of class `"coxph_pipeline"`, a list with elements:
#' \describe{
#'   \item{model_unadjusted}{The unadjusted `coxph` fit.}
#'   \item{model_adjusted}{The adjusted `coxph` fit, or `NULL` if `adjust_for`
#'     was not supplied.}
#'   \item{primary_model}{Whichever of the two above is the model of record
#'     (adjusted if available, else unadjusted) -- used for diagnostics.}
#'   \item{ph_test}{The [survival::cox.zph()] result for `primary_model`.}
#'   \item{influence}{A list with the dfbeta/dfbetas matrices (indexed to
#'     `used_rows`, in fitted-model order), `flagged_rows` (indices into the
#'     original `data` you passed in -- so `data[fit$influence$flagged_rows, ]`
#'     always works, even if `coxph()` silently dropped rows with missing
#'     values), `used_rows` (which original rows the model was actually
#'     fitted on), and `n_used`.}
#'   \item{functional_form}{A named list of data.frames (martingale residual
#'     vs. covariate value) for each continuous covariate in `covariates`.}
#'   \item{plots}{A named list of `ggplot` objects: `ph_plot`, `influence_plot`,
#'     and one `functional_form_<var>` entry per continuous covariate.}
#'   \item{verdicts}{A data.frame summarizing every diagnostic check, its
#'     pass/flag status, and a short interpretive note.}
#'   \item{alpha}{The significance threshold used.}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 200
#' dat <- data.frame(
#'   time = rexp(n, rate = 0.05),
#'   event = rbinom(n, 1, 0.7),
#'   age = rnorm(n, 55, 10),
#'   arm = factor(sample(c("control", "treatment"), n, replace = TRUE))
#' )
#' fit <- fit_coxph_pipeline(
#'   data = dat, time_var = "time", event_var = "event",
#'   covariates = "arm", adjust_for = "age"
#' )
#' print(fit)
#' }
#'
#' @export
fit_coxph_pipeline <- function(data,
                                time_var,
                                event_var,
                                covariates,
                                adjust_for = NULL,
                                id_var = NULL,
                                frailty = FALSE,
                                alpha = 0.05,
                                influence_cutoff_sd = 1) {

  assert_columns(data, c(time_var, event_var, covariates, adjust_for, id_var), "data")

  if (frailty && is.null(id_var)) {
    stop("`frailty = TRUE` requires `id_var` to be supplied.", call. = FALSE)
  }

  n <- nrow(data)

  # --- build formulas -------------------------------------------------------
  cluster_term <- NULL
  if (!is.null(id_var)) {
    cluster_term <- if (frailty) {
      sprintf("frailty(%s)", id_var)
    } else {
      sprintf("cluster(%s)", id_var)
    }
  }

  build_formula <- function(vars) {
    rhs <- c(vars, cluster_term)
    stats::as.formula(sprintf(
      "survival::Surv(%s, %s) ~ %s",
      time_var, event_var, paste(rhs, collapse = " + ")
    ))
  }

  formula_unadjusted <- build_formula(covariates)
  model_unadjusted <- survival::coxph(formula_unadjusted, data = data, x = TRUE)

  model_adjusted <- NULL
  if (!is.null(adjust_for) && length(adjust_for) > 0) {
    formula_adjusted <- build_formula(c(covariates, adjust_for))
    model_adjusted <- survival::coxph(formula_adjusted, data = data, x = TRUE)
  }

  primary_model <- if (!is.null(model_adjusted)) model_adjusted else model_unadjusted

  # coxph() silently drops rows with NA in any model variable (na.action =
  # na.omit by default). Every diagnostic below operates on residuals/fitted
  # values that are only as long as the *fitted* model, not the original
  # `data` -- so we track which original rows survived and index from that,
  # rather than assuming a 1:1 match with `data`.
  omitted <- primary_model$na.action
  used_rows <- if (!is.null(omitted)) seq_len(n)[-as.integer(omitted)] else seq_len(n)
  n_used <- length(used_rows)

  # --- proportional hazards test ---------------------------------------------
  ph_test <- survival::cox.zph(primary_model)
  ph_table <- as.data.frame(ph_test$table)
  ph_table$term <- rownames(ph_table)
  global_p <- ph_table$p[ph_table$term == "GLOBAL"]
  covariate_rows <- ph_table[ph_table$term != "GLOBAL", ]

  ph_verdicts <- do.call(rbind, lapply(seq_len(nrow(covariate_rows)), function(i) {
    row <- covariate_rows[i, ]
    make_verdict(
      check = sprintf("proportional_hazards[%s]", row$term),
      passed = row$p > alpha,
      statistic = row$chisq,
      p_value = row$p,
      note = if (row$p > alpha) {
        "No evidence against proportional hazards for this term."
      } else {
        "Schoenfeld residuals correlate with time -- hazard ratio for this term may not be constant; consider a time-varying coefficient or stratification."
      }
    )
  }))
  ph_verdicts <- rbind(
    ph_verdicts,
    make_verdict(
      check = "proportional_hazards[GLOBAL]",
      passed = global_p > alpha,
      statistic = ph_table$chisq[ph_table$term == "GLOBAL"],
      p_value = global_p,
      note = if (global_p > alpha) {
        "Global test finds no overall violation of the proportional hazards assumption."
      } else {
        "Global test flags a violation of proportional hazards somewhere in the model."
      }
    )
  )

  # --- influence diagnostics --------------------------------------------------
  dfbeta_mat <- stats::residuals(primary_model, type = "dfbeta")
  dfbetas_mat <- stats::residuals(primary_model, type = "dfbetas")
  if (is.null(dim(dfbeta_mat))) {
    dfbeta_mat <- matrix(dfbeta_mat, ncol = 1, dimnames = list(NULL, names(stats::coef(primary_model))[1]))
    dfbetas_mat <- matrix(dfbetas_mat, ncol = 1, dimnames = list(NULL, names(stats::coef(primary_model))[1]))
  }
  cutoff <- influence_cutoff_sd * 2 / sqrt(n_used)
  flagged_rows_fitted <- unique(which(abs(dfbetas_mat) > cutoff, arr.ind = TRUE)[, 1])
  # Map back to row positions in the original `data` the caller passed in,
  # so `data[fit$influence$flagged_rows, ]` works even when rows were dropped.
  flagged_rows <- used_rows[flagged_rows_fitted]
  prop_flagged <- length(flagged_rows) / n_used

  influence_verdict <- make_verdict(
    check = "influence[dfbetas]",
    passed = prop_flagged < 0.05,
    statistic = prop_flagged,
    p_value = NA_real_,
    note = sprintf(
      "%d of %d observations used in the fit (%.1f%%) exceed the |dfbetas| > %.3f cutoff.%s %s",
      length(flagged_rows), n_used, 100 * prop_flagged, cutoff,
      if (n_used < n) sprintf(" (%d row(s) were dropped due to missing values.)", n - n_used) else "",
      if (prop_flagged < 0.05) {
        "No single observation appears to be driving the fit."
      } else {
        "A non-trivial share of observations are individually influential -- inspect flagged rows before reporting hazard ratios."
      }
    )
  )

  # --- functional form check (martingale residuals) ---------------------------
  # Scan every continuous term actually in the primary (adjusted, if present)
  # model -- not just `covariates` -- since adjustment covariates are just as
  # subject to the linear log-hazard assumption.
  all_model_vars <- c(covariates, adjust_for)
  numeric_covs <- all_model_vars[vapply(data[all_model_vars], is.numeric, logical(1))]
  functional_form <- list()
  ff_verdicts <- NULL

  if (length(numeric_covs) > 0) {
    martingale_resid <- stats::residuals(primary_model, type = "martingale")
    for (v in numeric_covs) {
      df_v <- data.frame(
        covariate_value = data[[v]][used_rows],
        martingale_residual = martingale_resid
      )
      functional_form[[v]] <- df_v
      ff_verdicts <- rbind(
        ff_verdicts,
        make_verdict(
          check = sprintf("functional_form[%s]", v),
          passed = NA,
          statistic = NA_real_,

          p_value = NA_real_,
          note = "Martingale residuals plotted against this covariate; a loess smooth that is flat suggests the linear (log-hazard) form is adequate. Requires visual review -- not auto-scored."
        )
      )
    }
    ff_verdicts$verdict <- "review"
  }

  verdicts <- rbind(ph_verdicts, influence_verdict, ff_verdicts)
  rownames(verdicts) <- NULL

  # --- plots -------------------------------------------------------------------
  plots <- list()

  plots$ph_plot <- ggplot2::ggplot(
    covariate_rows,
    ggplot2::aes(x = stats::reorder(.data$term, .data$p), y = -log10(.data$p))
  ) +
    ggplot2::geom_col(fill = "#3B6E8F") +
    ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed", color = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Proportional hazards test by term",
      subtitle = sprintf("Dashed line = alpha = %.2f (values above the line indicate a violation)", alpha),
      x = NULL, y = expression(-log[10](p))
    ) +
    theme_omicsuite()

  influence_df <- data.frame(
    index = used_rows,
    max_abs_dfbetas = apply(abs(dfbetas_mat), 1, max),
    flagged = used_rows %in% flagged_rows
  )
  plots$influence_plot <- ggplot2::ggplot(
    influence_df,
    ggplot2::aes(x = .data$index, y = .data$max_abs_dfbetas, color = .data$flagged)
  ) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = cutoff, linetype = "dashed", color = "firebrick") +
    ggplot2::scale_color_manual(values = c(`FALSE` = "grey50", `TRUE` = "firebrick")) +
    ggplot2::labs(
      title = "Influence diagnostics (max |dfbetas| across terms)",
      x = "Observation index", y = "max |dfbetas|", color = "Flagged"
    ) +
    theme_omicsuite()

  for (v in names(functional_form)) {
    plots[[sprintf("functional_form_%s", v)]] <- ggplot2::ggplot(
      functional_form[[v]],
      ggplot2::aes(x = .data$covariate_value, y = .data$martingale_residual)
    ) +
      ggplot2::geom_point(alpha = 0.4, color = "grey40") +
      ggplot2::geom_smooth(method = "loess", formula = y ~ x, color = "#3B6E8F", se = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted") +
      ggplot2::labs(
        title = sprintf("Functional form check: %s", v),
        subtitle = "A flat smooth around zero supports a linear (log-hazard) effect",
        x = v, y = "Martingale residual"
      ) +
      theme_omicsuite()
  }

  structure(
    list(
      call = match.call(),
      model_unadjusted = model_unadjusted,
      model_adjusted = model_adjusted,
      primary_model = primary_model,
      ph_test = ph_test,
      influence = list(
        dfbeta = dfbeta_mat,
        dfbetas = dfbetas_mat,
        flagged_rows = flagged_rows,
        cutoff = cutoff,
        used_rows = used_rows,
        n_used = n_used
      ),
      functional_form = functional_form,
      plots = plots,
      verdicts = verdicts,
      alpha = alpha
    ),
    class = "coxph_pipeline"
  )
}

#' @export
print.coxph_pipeline <- function(x, ...) {
  cat("<omicsuite Cox PH pipeline>\n\n")
  cat("Unadjusted model:\n")
  print(x$model_unadjusted)
  if (!is.null(x$model_adjusted)) {
    cat("\nAdjusted model:\n")
    print(x$model_adjusted)
  }
  cat("\nDiagnostic verdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' @export
summary.coxph_pipeline <- function(object, ...) {
  list(
    unadjusted = summary(object$model_unadjusted),
    adjusted = if (!is.null(object$model_adjusted)) summary(object$model_adjusted) else NULL,
    ph_test = object$ph_test,
    verdicts = object$verdicts
  )
}

#' Plot a `coxph_pipeline` object
#'
#' Displays the diagnostic plots produced by [fit_coxph_pipeline()] one at a
#' time. Use `which` to select a subset instead of stepping through all of them.
#'
#' @param x A `coxph_pipeline` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.coxph_pipeline <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
