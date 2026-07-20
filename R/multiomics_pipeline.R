#' Integrate Multi-Omics Data Blocks via Regularized Generalized Canonical
#' Correlation Analysis
#'
#' Fits a regularized generalized canonical correlation analysis (RGCCA,
#' via [RGCCA::rgcca()]) across two or more omics blocks measured on the
#' same samples -- e.g. transcriptomics, proteomics, and methylation for the
#' same patients. Unlike stacking blocks into one matrix and running PCA,
#' RGCCA finds components that maximize covariation *between* blocks, so the
#' result reflects signal shared across omics layers rather than whichever
#' layer happens to have the most features.
#'
#' The loading stability check here is a case-resampling bootstrap built on
#' top of [RGCCA::rgcca()] rather than a call to `RGCCA::rgcca_bootstrap()`,
#' so the stability metric (sign-agreement rate for each block's top-loading
#' features) is transparent and doesn't depend on the exact internals of any
#' one RGCCA version.
#'
#' @param blocks A named list of matrices or data.frames, one per omics
#'   layer, samples in rows and features in columns. Row names are used to
#'   align samples across blocks; if blocks don't share identical sample
#'   sets, only the common samples are used (with a note in `verdicts`).
#' @param group Optional vector of group/phenotype labels for coloring the
#'   sample-score plot (e.g. a cancer subtype). Either a named vector keyed
#'   by sample ID, or an unnamed vector in the same order as the common
#'   samples -- if unnamed, order is assumed to match and is not verified.
#'   Not used in the RGCCA fit itself, only for plotting.
#' @param ncomp Integer. Number of components per block. Default `2`.
#' @param scheme Character. RGCCA connection scheme: `"factorial"`,
#'   `"centroid"`, or `"horst"`. Default `"factorial"`.
#' @param scale Logical. Standardize each block before fitting. Default
#'   `TRUE`.
#' @param tau Shrinkage parameter passed to [RGCCA::rgcca()]. Default
#'   `"optimal"`, which uses Schafer-Strimmer analytical shrinkage estimated
#'   separately per block.
#' @param n_boot Integer. Number of case-resampling bootstrap replicates
#'   for the loading stability check. Default `50`. Each replicate refits
#'   the full RGCCA model, so this is the slowest part of the pipeline --
#'   reduce it for exploration and increase it before reporting results.
#' @param boot_top_n Integer. Per block, only the top `boot_top_n` features
#'   by absolute component-1 loading are tracked for stability (tracking
#'   every feature in a high-dimensional block would make the bootstrap
#'   prohibitively slow without adding much information). Default `20`.
#' @param seed Optional integer seed.
#'
#' @return An object of class `"multiomics_pipeline"`, a list with elements:
#' \describe{
#'   \item{model}{The fitted `rgcca` object from [RGCCA::rgcca()].}
#'   \item{scores}{A named list of data.frames, one per block: sample scores
#'     on each component, with `sample_id` and `group` (if supplied)
#'     columns.}
#'   \item{variance_explained}{A data.frame: average variance explained
#'     (AVE) per block per component.}
#'   \item{stability}{A data.frame: per-block, per-feature bootstrap
#'     sign-agreement rate for the top-loading features on component 1.}
#'   \item{plots}{A named list of `ggplot` objects: `block_scores`,
#'     `variance_explained`, `stability`.}
#'   \item{verdicts}{A data.frame summarizing sample alignment, variance
#'     explained, and loading stability, in the same style as
#'     [fit_coxph_pipeline()].}
#' }
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' n <- 40
#' block1 <- matrix(rnorm(n * 15), nrow = n, dimnames = list(paste0("s", 1:n), NULL))
#' block2 <- matrix(rnorm(n * 10), nrow = n, dimnames = list(paste0("s", 1:n), NULL))
#' # give the blocks some shared signal
#' shared <- rnorm(n)
#' block1[, 1] <- block1[, 1] + shared
#' block2[, 1] <- block2[, 1] + shared
#'
#' fit <- integrate_multiomics(
#'   blocks = list(omics_a = block1, omics_b = block2),
#'   ncomp = 2, n_boot = 20
#' )
#' print(fit)
#' }
#'
#' @export
integrate_multiomics <- function(blocks,
                                  group = NULL,
                                  ncomp = 2,
                                  scheme = c("factorial", "centroid", "horst"),
                                  scale = TRUE,
                                  tau = "optimal",
                                  n_boot = 50,
                                  boot_top_n = 20,
                                  seed = NULL) {

  scheme <- match.arg(scheme)

  if (is.null(names(blocks)) || any(names(blocks) == "")) {
    stop("`blocks` must be a fully named list (one name per omics layer).", call. = FALSE)
  }
  if (length(blocks) < 2) {
    stop("`blocks` must contain at least two omics layers.", call. = FALSE)
  }
  blocks <- lapply(blocks, as.matrix)
  if (any(vapply(blocks, function(b) is.null(rownames(b)), logical(1)))) {
    stop("Every block must have row names giving the sample IDs, so samples can be aligned across blocks.", call. = FALSE)
  }

  if (!is.null(seed)) set.seed(seed)

  # --- align samples across blocks -----------------------------------------------
  common_samples <- Reduce(intersect, lapply(blocks, rownames))
  n_dropped <- sum(vapply(blocks, nrow, integer(1))) - length(common_samples) * length(blocks)
  blocks <- lapply(blocks, function(b) b[common_samples, , drop = FALSE])

  alignment_verdict <- make_verdict(
    check = "sample_alignment",
    passed = n_dropped == 0,
    statistic = length(common_samples),
    p_value = NA_real_,
    note = if (n_dropped == 0) {
      sprintf("All %d samples were present in every block.", length(common_samples))
    } else {
      sprintf(
        "%d sample-block row(s) outside the %d samples common to every block were dropped -- check that this matches your expectation before interpreting scores.",
        n_dropped, length(common_samples)
      )
    }
  )

  # --- group labels for plotting ---------------------------------------------------
  group_df <- NULL
  if (!is.null(group)) {
    if (!is.null(names(group))) {
      group_df <- data.frame(sample_id = common_samples, group = group[common_samples], stringsAsFactors = FALSE)
    } else if (length(group) == length(common_samples)) {
      group_df <- data.frame(sample_id = common_samples, group = group, stringsAsFactors = FALSE)
    } else {
      warning("`group` is unnamed and its length doesn't match the number of common samples -- ignoring it for plotting.", call. = FALSE)
    }
  }

  # --- fit RGCCA --------------------------------------------------------------------
  n_blocks <- length(blocks)
  connection <- matrix(1, n_blocks, n_blocks)
  diag(connection) <- 0

  fit_rgcca_once <- function(block_list) {
    RGCCA::rgcca(
      blocks = block_list,
      connection = connection,
      tau = tau,
      ncomp = ncomp,
      scheme = scheme,
      scale = scale,
      quiet = TRUE
    )
  }

  model <- tryCatch(
    fit_rgcca_once(blocks),
    error = function(e) {
      stop(
        "RGCCA::rgcca() failed to fit: ", conditionMessage(e),
        ". If this looks like an argument-name mismatch, check packageVersion(\"RGCCA\") -- ",
        "the rgcca() interface has changed across major versions.",
        call. = FALSE
      )
    }
  )

  # --- variance explained -----------------------------------------------------------
  ave_df <- do.call(rbind, lapply(names(blocks), function(bname) {
    ave_vec <- model$AVE$AVE_X[[bname]]
    data.frame(
      block = bname,
      component = seq_along(ave_vec),
      variance_explained = as.numeric(ave_vec)
    )
  }))

  ave_verdicts <- do.call(rbind, lapply(names(blocks), function(bname) {
    comp1_ave <- ave_df$variance_explained[ave_df$block == bname & ave_df$component == 1]
    make_verdict(
      check = sprintf("variance_explained[%s]", bname),
      passed = NA,
      statistic = comp1_ave,
      p_value = NA_real_,
      note = sprintf(
        "Component 1 explains %.1f%% of variance within this block. %s",
        100 * comp1_ave,
        if (comp1_ave < 0.05) {
          "This is low -- this block may be contributing mostly noise to the shared components; consider whether it belongs in the integration."
        } else {
          "This block is contributing a meaningful share of its own variance to the shared components."
        }
      )
    )
  }))
  ave_verdicts$verdict <- "info"

  # --- sample scores ------------------------------------------------------------------
  scores <- lapply(names(blocks), function(bname) {
    score_mat <- as.data.frame(model$Y[[bname]])
    names(score_mat) <- paste0("comp", seq_len(ncol(score_mat)))
    score_mat$sample_id <- common_samples
    if (!is.null(group_df)) {
      score_mat <- merge(score_mat, group_df, by = "sample_id", sort = FALSE)
    }
    score_mat
  })
  names(scores) <- names(blocks)

  # --- loading stability via manual case-resampling bootstrap ------------------------
  original_loadings <- lapply(names(blocks), function(bname) model$a[[bname]][, 1])
  names(original_loadings) <- names(blocks)

  top_features <- lapply(names(blocks), function(bname) {
    loadings <- original_loadings[[bname]]
    ord <- order(abs(loadings), decreasing = TRUE)
    names(loadings)[ord[seq_len(min(boot_top_n, length(loadings)))]]
  })
  names(top_features) <- names(blocks)

  n_samples <- length(common_samples)
  sign_match_counts <- lapply(names(blocks), function(bname) {
    stats::setNames(rep(0L, length(top_features[[bname]])), top_features[[bname]])
  })
  names(sign_match_counts) <- names(blocks)
  n_successful_boots <- 0L

  for (b in seq_len(n_boot)) {
    resample_idx <- sample(seq_len(n_samples), size = n_samples, replace = TRUE)
    boot_blocks <- lapply(blocks, function(mat) {
      resampled <- mat[resample_idx, , drop = FALSE]
      rownames(resampled) <- paste0("boot_", seq_len(n_samples))
      resampled
    })
    boot_fit <- tryCatch(fit_rgcca_once(boot_blocks), error = function(e) NULL)
    if (is.null(boot_fit)) next
    n_successful_boots <- n_successful_boots + 1L

    for (bname in names(blocks)) {
      boot_loadings <- boot_fit$a[[bname]][, 1]
      feats <- top_features[[bname]]
      for (feat in feats) {
        if (feat %in% names(boot_loadings)) {
          same_sign <- sign(boot_loadings[[feat]]) == sign(original_loadings[[bname]][[feat]])
          if (isTRUE(same_sign)) {
            sign_match_counts[[bname]][feat] <- sign_match_counts[[bname]][feat] + 1L
          }
        }
      }
    }
  }

  stability <- do.call(rbind, lapply(names(blocks), function(bname) {
    data.frame(
      block = bname,
      feature = names(sign_match_counts[[bname]]),
      original_loading = original_loadings[[bname]][names(sign_match_counts[[bname]])],
      stability = if (n_successful_boots > 0) sign_match_counts[[bname]] / n_successful_boots else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(stability) <- NULL

  stability_verdicts <- do.call(rbind, lapply(names(blocks), function(bname) {
    block_stability <- stability$stability[stability$block == bname]
    mean_stability <- mean(block_stability, na.rm = TRUE)
    make_verdict(
      check = sprintf("loading_stability[%s]", bname),
      passed = mean_stability > 0.8,
      statistic = mean_stability,
      p_value = NA_real_,
      note = sprintf(
        "Top %d loadings had a mean bootstrap sign-agreement rate of %.1f%% across %d successful resamples. %s",
        length(block_stability), 100 * mean_stability, n_successful_boots,
        if (mean_stability > 0.8) {
          "Loadings for this block appear stable to resampling."
        } else {
          "A non-trivial share of top loadings flip sign under resampling -- treat the specific feature ranking within this block cautiously, even if the block-level integration is otherwise sound."
        }
      )
    )
  }))

  verdicts <- rbind(alignment_verdict, ave_verdicts, stability_verdicts)
  rownames(verdicts) <- NULL

  # --- plots -------------------------------------------------------------------
  plots <- list()

  scores_long <- do.call(rbind, lapply(names(scores), function(bname) {
    df <- scores[[bname]]
    df$block <- bname
    df
  }))
  color_aes <- if (!is.null(group_df)) ggplot2::aes(color = .data$group) else ggplot2::aes()
  plots$block_scores <- ggplot2::ggplot(
    scores_long, ggplot2::aes(x = .data$comp1, y = .data$comp2)
  ) +
    color_aes +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::facet_wrap(~block, scales = "free") +
    ggplot2::labs(
      title = "Sample scores on the first two shared components, by block",
      x = "Component 1", y = "Component 2", color = if (!is.null(group_df)) "Group" else NULL
    ) +
    theme_omicsuite()

  plots$variance_explained <- ggplot2::ggplot(
    ave_df, ggplot2::aes(x = factor(.data$component), y = .data$variance_explained, fill = .data$block)
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::labs(
      title = "Average variance explained (AVE) per block",
      x = "Component", y = "Variance explained", fill = "Block"
    ) +
    theme_omicsuite()

  plots$stability <- ggplot2::ggplot(
    stability, ggplot2::aes(x = stats::reorder(.data$feature, .data$stability), y = .data$stability)
  ) +
    ggplot2::geom_col(fill = "#3B6E8F") +
    ggplot2::geom_hline(yintercept = 0.8, linetype = "dashed", color = "firebrick") +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~block, scales = "free_y") +
    ggplot2::labs(
      title = "Bootstrap sign-agreement rate, top loadings by block",
      subtitle = "Dashed line = 0.8 (informal stability threshold)",
      x = NULL, y = "Sign-agreement rate"
    ) +
    theme_omicsuite()

  structure(
    list(
      call = match.call(),
      model = model,
      scores = scores,
      variance_explained = ave_df,
      stability = stability,
      plots = plots,
      verdicts = verdicts
    ),
    class = "multiomics_pipeline"
  )
}

#' @export
print.multiomics_pipeline <- function(x, ...) {
  cat("<omicsuite multi-omics RGCCA integration pipeline>\n\n")
  cat(sprintf("Blocks: %s\n", paste(names(x$scores), collapse = ", ")))
  cat("\nDiagnostic verdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' @export
summary.multiomics_pipeline <- function(object, ...) {
  list(
    variance_explained = object$variance_explained,
    stability = object$stability,
    verdicts = object$verdicts
  )
}

#' Plot a `multiomics_pipeline` object
#'
#' @param x A `multiomics_pipeline` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.multiomics_pipeline <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
