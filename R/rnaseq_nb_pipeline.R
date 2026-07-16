#' Fit a Bayesian Hierarchical Negative Binomial Mixed Model for RNA-Seq Counts
#'
#' Fits a multilevel negative binomial model to long-format RNA-seq count
#' data (one row per gene x sample) via [brms::brm()]. The effect of primary
#' biological interest (`condition_var`) gets a gene-level random slope, so
#' genes borrow statistical strength from one another (partial pooling)
#' rather than being tested one at a time in isolation. An optional
#' subject-level random intercept accounts for repeated-measures /
#' longitudinal correlation within individuals.
#'
#' `brms` (and the Stan toolchain it depends on) is listed as a `Suggests`,
#' not a hard dependency of `omicsuite`, since it requires a working C++
#' compiler. Install it with `install.packages("brms")` before calling this
#' function.
#'
#' @param data A long-format data.frame: one row per gene x sample
#'   combination.
#' @param gene_var Character. Column identifying the gene. Gets a
#'   gene-level random intercept and random slope for `condition_var`.
#' @param count_var Character. Column of raw (unnormalized) RNA-seq counts.
#' @param condition_var Character. The fixed effect of primary biological
#'   interest (e.g. treatment arm, genotype). Given a gene-level random
#'   slope so per-gene effects are partially pooled toward the population
#'   average. Should be numeric or a two-level factor: the shrinkage
#'   summary and plot currently reflect only the first non-reference
#'   coefficient, so a factor with more than two levels will emit a warning
#'   and only show one contrast.
#' @param covariates Character vector of additional fixed effects to adjust
#'   for (population-level only -- no gene-specific slope). Default `NULL`.
#' @param subject_var Optional character. Subject/individual identifier for
#'   repeated-measures or longitudinal designs; gets a subject-level random
#'   intercept. Default `NULL` (no repeated-measures structure).
#' @param sample_var Character. Column identifying the sample -- used to
#'   compute the library-size offset when `offset_var` is not supplied.
#' @param offset_var Optional character. Column of pre-computed
#'   normalization factors (raw scale, e.g. total counts or a DESeq2/edgeR
#'   size factor) to use as the model offset. If `NULL` (the default), a
#'   simple total-count-per-sample offset is computed internally -- adequate
#'   for a first pass, but a proper size-factor normalization (e.g.
#'   median-of-ratios) is recommended for a publication-grade analysis and
#'   can be supplied here instead.
#' @param chains,iter,warmup,cores Passed to [brms::brm()]. Defaults
#'   (`chains = 4, iter = 2000, warmup = 1000, cores = 1`) are a reasonable
#'   starting point for exploration; increase `cores` to `chains` if your
#'   machine allows running them in parallel, and increase `iter` if the
#'   convergence verdict flags a low effective-sample-size ratio.
#' @param seed Optional integer seed passed to [brms::brm()]. Default `NA`
#'   (brms's own default -- no fixed seed). Note this must be `NA`, not
#'   `NULL`, since `brms::brm()` cannot coerce `NULL` to a numeric value.
#' @param cache_file Optional path. If it exists, the cached fit is loaded
#'   via [readRDS()] and no model is refit. If it does not exist, the model
#'   is fit and then saved there via [saveRDS()]. Recommended for any model
#'   that takes more than a minute or two to fit, since Stan model
#'   compilation and sampling both restart from scratch otherwise.
#' @param ... Additional arguments passed through to [brms::brm()].
#'
#' @return An object of class `"rnaseq_nb_pipeline"`, a list with elements:
#' \describe{
#'   \item{model}{The fitted `brmsfit` object.}
#'   \item{convergence}{A list with `max_rhat`, `min_neff_ratio`, and the
#'     full per-parameter summaries.}
#'   \item{shrinkage}{A data.frame with one row per gene: the partially
#'     pooled (population + gene-level deviation) estimate of the
#'     `condition_var` effect, its credible interval, and the population-level
#'     estimate for comparison.}
#'   \item{dispersion}{A list with the posterior summary of the negative
#'     binomial shape parameter.}
#'   \item{plots}{A named list of `ggplot` objects: `pp_check`,
#'     `shrinkage_plot`, and `rhat_plot`.}
#'   \item{verdicts}{A data.frame summarizing convergence, posterior
#'     predictive fit, and dispersion, in the same style as
#'     [fit_coxph_pipeline()].}
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("brms", quietly = TRUE)) {
#'   set.seed(1)
#'   n_genes <- 8
#'   n_samples <- 12
#'   dat <- expand.grid(gene = paste0("gene_", seq_len(n_genes)),
#'                       sample = paste0("sample_", seq_len(n_samples)))
#'   dat$condition <- factor(rep(rep(c("control", "treated"), each = n_samples / 2),
#'                                n_genes))
#'   dat$count <- stats::rnbinom(nrow(dat), mu = 50, size = 5)
#'
#'   fit <- fit_rnaseq_nb_pipeline(
#'     data = dat, gene_var = "gene", count_var = "count",
#'     condition_var = "condition", sample_var = "sample",
#'     chains = 1, iter = 200, warmup = 100
#'   )
#'   print(fit)
#' }
#' }
#'
#' @export
fit_rnaseq_nb_pipeline <- function(data,
                                    gene_var,
                                    count_var,
                                    condition_var,
                                    covariates = NULL,
                                    subject_var = NULL,
                                    sample_var,
                                    offset_var = NULL,
                                    chains = 4,
                                    iter = 2000,
                                    warmup = 1000,
                                    cores = 1,
                                    seed = NA,
                                    cache_file = NULL,
                                    ...) {

  if (!brms_is_available()) {
    stop(
      "fit_rnaseq_nb_pipeline() requires the 'brms' package (and a working ",
      "C++ toolchain for Stan), which is not installed. Install it with ",
      "install.packages(\"brms\") and try again.",
      call. = FALSE
    )
  }

  assert_columns(
    data,
    c(gene_var, count_var, condition_var, covariates, subject_var, sample_var, offset_var),
    "data"
  )

  if (!is.null(cache_file) && file.exists(cache_file)) {
    model <- readRDS(cache_file)
  } else {

    # --- offset -----------------------------------------------------------
    if (is.null(offset_var)) {
      lib_size <- stats::aggregate(
        stats::as.formula(sprintf("%s ~ %s", count_var, sample_var)),
        data = data, FUN = sum
      )
      names(lib_size)[2] <- ".omicsuite_lib_size"
      data <- merge(data, lib_size, by = sample_var, all.x = TRUE, sort = FALSE)
      data$.omicsuite_log_offset <- log(pmax(data$.omicsuite_lib_size, 1))
    } else {
      data$.omicsuite_log_offset <- log(pmax(data[[offset_var]], .Machine$double.eps))
    }

    # --- formula ------------------------------------------------------------
    fixed_terms <- c(condition_var, covariates)
    random_terms <- sprintf("(1 + %s | %s)", condition_var, gene_var)
    if (!is.null(subject_var)) {
      random_terms <- c(random_terms, sprintf("(1 | %s)", subject_var))
    }
    rhs <- paste(c(fixed_terms, "offset(.omicsuite_log_offset)", random_terms), collapse = " + ")
    model_formula <- stats::as.formula(sprintf("%s ~ %s", count_var, rhs))

    model <- brms::brm(
      formula = model_formula,
      data = data,
      family = brms::negbinomial(),
      chains = chains,
      iter = iter,
      warmup = warmup,
      cores = cores,
      seed = seed,
      ...
    )

    if (!is.null(cache_file)) {
      saveRDS(model, cache_file)
    }
  }

  # --- convergence ------------------------------------------------------------
  rhat_vec <- brms::rhat(model)
  neff_vec <- brms::neff_ratio(model)
  max_rhat <- max(rhat_vec, na.rm = TRUE)
  min_neff_ratio <- min(neff_vec, na.rm = TRUE)

  convergence_verdict <- make_verdict(
    check = "mcmc_convergence",
    passed = max_rhat < 1.01 && min_neff_ratio > 0.1,
    statistic = max_rhat,
    p_value = NA_real_,
    note = sprintf(
      "Max Rhat = %.4f, min effective-sample-size ratio = %.3f. %s",
      max_rhat, min_neff_ratio,
      if (max_rhat < 1.01 && min_neff_ratio > 0.1) {
        "Chains appear to have converged and mixed adequately."
      } else {
        "Convergence is questionable -- consider more iterations, more chains, or reparameterizing before trusting the posterior."
      }
    )
  )

  # --- dispersion (negative binomial shape parameter) --------------------------
  shape_summary <- brms::posterior_summary(model, variable = "shape")
  dispersion_verdict <- make_verdict(
    check = "dispersion[shape]",
    passed = NA,
    statistic = shape_summary[1, "Estimate"],
    p_value = NA_real_,
    note = sprintf(
      "Posterior mean shape = %.2f (95%% CI [%.2f, %.2f]). %s",
      shape_summary[1, "Estimate"], shape_summary[1, "Q2.5"], shape_summary[1, "Q97.5"],
      if (shape_summary[1, "Estimate"] > 30) {
        "A large shape parameter means the fitted distribution is close to Poisson -- the extra overdispersion parameter may not be doing much work here."
      } else {
        "Overdispersion beyond Poisson is clearly present and the negative binomial family is earning its keep."
      }
    )
  )
  dispersion_verdict$verdict <- "info"

  ppcheck_verdict <- make_verdict(
    check = "posterior_predictive_check",
    passed = NA,
    statistic = NA_real_,
    p_value = NA_real_,
    note = "Compare the plots$pp_check overlay of observed vs. replicated counts. Requires visual review -- not auto-scored."
  )
  ppcheck_verdict$verdict <- "review"

  verdicts <- rbind(convergence_verdict, dispersion_verdict, ppcheck_verdict)
  rownames(verdicts) <- NULL

  # --- shrinkage (gene-level partial pooling of condition_var) ------------------
  fixef_est <- brms::fixef(model)
  condition_coef_names <- grep(
    sprintf("^%s", condition_var),
    rownames(fixef_est), value = TRUE
  )
  primary_coef <- condition_coef_names[1]
  if (length(condition_coef_names) > 1) {
    warning(sprintf(
      "`%s` has more than one non-reference level (%s); the shrinkage data.frame and shrinkage_plot reflect only `%s`. Multi-level factors for condition_var aren't fully supported in this release -- consider recoding to a two-level contrast of primary interest.",
      condition_var, paste(condition_coef_names, collapse = ", "), primary_coef
    ), call. = FALSE)
  }
  population_estimate <- fixef_est[primary_coef, "Estimate"]

  ranef_arr <- brms::ranef(model)[[gene_var]]
  gene_slope_dim <- grep(sprintf("^%s", condition_var), dimnames(ranef_arr)[[3]], value = TRUE)[1]

  shrinkage <- data.frame(
    gene = dimnames(ranef_arr)[[1]],
    gene_deviation = ranef_arr[, "Estimate", gene_slope_dim],
    gene_deviation_lower = ranef_arr[, "Q2.5", gene_slope_dim],
    gene_deviation_upper = ranef_arr[, "Q97.5", gene_slope_dim],
    stringsAsFactors = FALSE
  )
  shrinkage$pooled_estimate <- population_estimate + shrinkage$gene_deviation
  shrinkage$population_estimate <- population_estimate
  shrinkage <- shrinkage[order(shrinkage$pooled_estimate), ]
  shrinkage$gene <- factor(shrinkage$gene, levels = shrinkage$gene)

  # --- plots -------------------------------------------------------------------
  plots <- list()

  plots$pp_check <- tryCatch(
    brms::pp_check(model, ndraws = 100) +
      ggplot2::labs(title = "Posterior predictive check") +
      theme_omicsuite(),
    error = function(e) NULL
  )

  plots$shrinkage_plot <- ggplot2::ggplot(
    shrinkage,
    ggplot2::aes(x = .data$pooled_estimate, y = .data$gene)
  ) +
    ggplot2::geom_vline(xintercept = population_estimate, linetype = "dashed", color = "firebrick") +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = .data$population_estimate + .data$gene_deviation_lower,
        xmax = .data$population_estimate + .data$gene_deviation_upper
      ),
      height = 0, color = "grey50"
    ) +
    ggplot2::geom_point(color = "#3B6E8F", size = 2) +
    ggplot2::labs(
      title = sprintf("Gene-level %s effect (partially pooled)", condition_var),
      subtitle = "Dashed line = population-level estimate; points shrink toward it as evidence per gene weakens",
      x = "Effect estimate (log scale)", y = NULL
    ) +
    theme_omicsuite()

  rhat_df <- data.frame(parameter = names(rhat_vec), rhat = as.numeric(rhat_vec))
  plots$rhat_plot <- ggplot2::ggplot(rhat_df, ggplot2::aes(x = .data$rhat)) +
    ggplot2::geom_histogram(bins = 30, fill = "#3B6E8F") +
    ggplot2::geom_vline(xintercept = 1.01, linetype = "dashed", color = "firebrick") +
    ggplot2::labs(
      title = "Rhat distribution across all parameters",
      subtitle = "Dashed line = 1.01 (common convergence threshold)",
      x = "Rhat", y = "Count of parameters"
    ) +
    theme_omicsuite()

  structure(
    list(
      call = match.call(),
      model = model,
      convergence = list(
        max_rhat = max_rhat,
        min_neff_ratio = min_neff_ratio,
        rhat = rhat_vec,
        neff_ratio = neff_vec
      ),
      dispersion = list(shape_summary = shape_summary),
      shrinkage = shrinkage,
      plots = plots,
      verdicts = verdicts
    ),
    class = "rnaseq_nb_pipeline"
  )
}

#' @export
print.rnaseq_nb_pipeline <- function(x, ...) {
  cat("<omicsuite RNA-seq negative binomial mixed model pipeline>\n\n")
  print(x$model)
  cat("\nDiagnostic verdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' @export
summary.rnaseq_nb_pipeline <- function(object, ...) {
  list(
    model_summary = summary(object$model),
    convergence = object$convergence,
    dispersion = object$dispersion,
    verdicts = object$verdicts
  )
}

#' Plot an `rnaseq_nb_pipeline` object
#'
#' @param x An `rnaseq_nb_pipeline` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.rnaseq_nb_pipeline <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    if (!is.null(x$plots[[nm]])) print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
