#' Simulate a Stochastic SIR or SEIR Epidemic via the Gillespie Algorithm
#'
#' Runs one or more exact stochastic realizations of a compartmental epidemic
#' model (SIR or SEIR) using Gillespie's direct method (the stochastic
#' simulation algorithm), then summarizes the ensemble: a median trajectory
#' with a 5-95% envelope, the distribution of final outbreak size and peak
#' timing across realizations, the basic reproduction number implied by the
#' rate parameters, and the proportion of realizations that fade out early
#' by chance rather than becoming a sustained outbreak -- a distinction a
#' deterministic ODE model can't show you at all.
#'
#' @param model Character. Either `"SIR"` or `"SEIR"`.
#' @param initial_state Named numeric vector of starting compartment counts.
#'   Must include `S`, `I`, `R` for `"SIR"`, or `S`, `E`, `I`, `R` for
#'   `"SEIR"`.
#' @param params Named list of rate parameters. Must include `beta`
#'   (transmission rate) and `gamma` (recovery rate) for `"SIR"`; `"SEIR"`
#'   additionally requires `sigma` (the rate of progression from exposed to
#'   infectious, i.e. `1 / sigma` is the mean incubation period).
#' @param t_max Numeric. Maximum simulation time.
#' @param n_sim Integer. Number of independent stochastic realizations to
#'   run. Default `100`, enough to characterize the extinction probability
#'   and the envelope around the median trajectory without being slow.
#' @param seed Optional integer seed for reproducibility.
#' @param grid_points Integer. Number of time points used when interpolating
#'   realizations onto a common grid for the summary envelope and plots.
#'   Default `200`.
#' @param extinction_threshold Numeric in `(0, 1)`. A realization is
#'   classified as an early stochastic fadeout if its final outbreak size
#'   (total ever infected) is below this fraction of the total population.
#'   Default `0.05`.
#'
#' @return An object of class `"gillespie_epidemic"`, a list with elements:
#' \describe{
#'   \item{simulations}{A list of `n_sim` data.frames, one per realization,
#'     each with columns `time` and one column per compartment.}
#'   \item{summary}{A named list (one entry per compartment) of data.frames
#'     with `time`, `median`, `lower` (5th percentile), and `upper` (95th
#'     percentile) across realizations, on a common time grid.}
#'   \item{final_size, peak_time, peak_size}{Numeric vectors, one entry per
#'     realization: total ever infected, time of peak infectious count, and
#'     peak infectious count.}
#'   \item{r0}{The basic reproduction number implied by `params`
#'     (`beta / gamma`, the same threshold quantity for both SIR and SEIR).}
#'   \item{prop_extinct}{Proportion of realizations classified as an early
#'     stochastic fadeout.}
#'   \item{plots}{A named list of `ggplot` objects: `trajectory_plot`,
#'     `final_size_hist`, `peak_time_hist`.}
#'   \item{verdicts}{A data.frame summarizing the basic reproduction number
#'     and the extinction probability, in the same style as
#'     [fit_coxph_pipeline()].}
#' }
#'
#' @examples
#' sir_fit <- simulate_gillespie_epidemic(
#'   model = "SIR",
#'   initial_state = c(S = 999, I = 1, R = 0),
#'   params = list(beta = 0.4, gamma = 0.1),
#'   t_max = 100, n_sim = 50, seed = 1
#' )
#' print(sir_fit)
#'
#' seir_fit <- simulate_gillespie_epidemic(
#'   model = "SEIR",
#'   initial_state = c(S = 999, E = 0, I = 1, R = 0),
#'   params = list(beta = 0.4, sigma = 0.2, gamma = 0.1),
#'   t_max = 150, n_sim = 50, seed = 1
#' )
#' print(seir_fit)
#'
#' @export
simulate_gillespie_epidemic <- function(model = c("SIR", "SEIR"),
                                         initial_state,
                                         params,
                                         t_max,
                                         n_sim = 100,
                                         seed = NULL,
                                         grid_points = 200,
                                         extinction_threshold = 0.05) {

  model <- match.arg(model)
  if (!is.null(seed)) set.seed(seed)

  required_state <- if (model == "SIR") c("S", "I", "R") else c("S", "E", "I", "R")
  missing_state <- setdiff(required_state, names(initial_state))
  if (length(missing_state) > 0) {
    stop(sprintf(
      "initial_state is missing required compartment(s): %s",
      paste(missing_state, collapse = ", ")
    ), call. = FALSE)
  }

  required_params <- if (model == "SIR") c("beta", "gamma") else c("beta", "sigma", "gamma")
  missing_params <- setdiff(required_params, names(params))
  if (length(missing_params) > 0) {
    stop(sprintf(
      "params is missing required rate(s): %s",
      paste(missing_params, collapse = ", ")
    ), call. = FALSE)
  }

  simulations <- vector("list", n_sim)
  for (i in seq_len(n_sim)) {
    simulations[[i]] <- if (model == "SIR") {
      run_gillespie_sir(
        S0 = initial_state[["S"]], I0 = initial_state[["I"]], R0 = initial_state[["R"]],
        beta = params[["beta"]], gamma = params[["gamma"]], t_max = t_max
      )
    } else {
      run_gillespie_seir(
        S0 = initial_state[["S"]], E0 = initial_state[["E"]],
        I0 = initial_state[["I"]], R0 = initial_state[["R"]],
        beta = params[["beta"]], sigma = params[["sigma"]], gamma = params[["gamma"]],
        t_max = t_max
      )
    }
  }

  # --- ensemble summary on a common time grid ----------------------------------
  grid <- seq(0, t_max, length.out = grid_points)
  summary_list <- list()
  for (comp in required_state) {
    mat <- vapply(simulations, function(df) {
      stats::approx(df$time, df[[comp]], xout = grid, method = "constant", rule = 2)$y
    }, numeric(length(grid)))
    summary_list[[comp]] <- data.frame(
      time = grid,
      median = apply(mat, 1, stats::median),
      lower = apply(mat, 1, function(v) stats::quantile(v, probs = 0.05)),
      upper = apply(mat, 1, function(v) stats::quantile(v, probs = 0.95))
    )
  }

  # --- final size / peak timing across realizations ----------------------------
  n_total <- sum(initial_state)
  final_size <- vapply(simulations, function(df) n_total - df$S[nrow(df)], numeric(1))
  peak_time <- vapply(simulations, function(df) df$time[which.max(df$I)], numeric(1))
  peak_size <- vapply(simulations, function(df) max(df$I), numeric(1))

  is_extinct <- final_size < (extinction_threshold * n_total)
  prop_extinct <- mean(is_extinct)

  r0 <- params[["beta"]] / params[["gamma"]]

  # --- verdicts -----------------------------------------------------------------
  r0_verdict <- make_verdict(
    check = "basic_reproduction_number",
    passed = NA, statistic = r0, p_value = NA_real_,
    note = sprintf(
      "R0 = beta / gamma = %.2f. %s",
      r0,
      if (r0 > 1) {
        "Above the epidemic threshold -- sustained transmission is expected in a well-mixed population."
      } else {
        "At or below the epidemic threshold -- the infection is expected to die out without a major outbreak, on average."
      }
    )
  )
  r0_verdict$verdict <- "info"

  extinction_verdict <- make_verdict(
    check = "stochastic_extinction",
    passed = NA, statistic = prop_extinct, p_value = NA_real_,
    note = sprintf(
      "%.1f%% of %d realizations resulted in early stochastic fadeout (fewer than %.0f%% of the population ever infected). %s",
      100 * prop_extinct, n_sim, 100 * extinction_threshold,
      if (prop_extinct > 0.3) {
        "Extinction risk is substantial even where R0 favors an outbreak -- worth reporting alongside the deterministic R0 rather than instead of it."
      } else {
        "Most realizations produced a sustained outbreak; stochastic fadeout is a minor consideration here."
      }
    )
  )
  extinction_verdict$verdict <- "info"

  verdicts <- rbind(r0_verdict, extinction_verdict)
  rownames(verdicts) <- NULL

  # --- plots -------------------------------------------------------------------
  plots <- list()

  infectious_summary <- summary_list$I
  n_spaghetti <- min(20, n_sim)
  spaghetti <- do.call(rbind, lapply(seq_len(n_spaghetti), function(i) {
    df <- simulations[[i]]
    data.frame(time = df$time, I = df$I, realization = i)
  }))

  plots$trajectory_plot <- ggplot2::ggplot() +
    ggplot2::geom_step(
      data = spaghetti,
      ggplot2::aes(x = .data$time, y = .data$I, group = .data$realization),
      color = "grey70", alpha = 0.5, linewidth = 0.3
    ) +
    ggplot2::geom_ribbon(
      data = infectious_summary,
      ggplot2::aes(x = .data$time, ymin = .data$lower, ymax = .data$upper),
      fill = "#3B6E8F", alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = infectious_summary,
      ggplot2::aes(x = .data$time, y = .data$median),
      color = "#3B6E8F", linewidth = 1
    ) +
    ggplot2::labs(
      title = sprintf("%s: infectious compartment across %d realizations", model, n_sim),
      subtitle = sprintf("Grey lines: %d individual realizations. Band: 5th-95th percentile. Solid line: median.", n_spaghetti),
      x = "Time", y = "Infectious (I)"
    ) +
    theme_omicsuite()

  plots$final_size_hist <- ggplot2::ggplot(
    data.frame(final_size = final_size),
    ggplot2::aes(x = .data$final_size)
  ) +
    ggplot2::geom_histogram(bins = 30, fill = "#3B6E8F") +
    ggplot2::labs(
      title = "Distribution of final outbreak size",
      subtitle = sprintf("Across %d realizations; population = %d", n_sim, n_total),
      x = "Total ever infected", y = "Count of realizations"
    ) +
    theme_omicsuite()

  plots$peak_time_hist <- ggplot2::ggplot(
    data.frame(peak_time = peak_time[!is_extinct]),
    ggplot2::aes(x = .data$peak_time)
  ) +
    ggplot2::geom_histogram(bins = 30, fill = "#3B6E8F") +
    ggplot2::labs(
      title = "Distribution of peak infection timing",
      subtitle = "Non-extinct realizations only",
      x = "Time of peak infectious count", y = "Count of realizations"
    ) +
    theme_omicsuite()

  structure(
    list(
      call = match.call(),
      model = model,
      initial_state = initial_state,
      params = params,
      t_max = t_max,
      n_sim = n_sim,
      simulations = simulations,
      summary = summary_list,
      final_size = final_size,
      peak_time = peak_time,
      peak_size = peak_size,
      r0 = r0,
      prop_extinct = prop_extinct,
      is_extinct = is_extinct,
      extinction_threshold = extinction_threshold,
      plots = plots,
      verdicts = verdicts
    ),
    class = "gillespie_epidemic"
  )
}

#' @keywords internal
#' @noRd
run_gillespie_sir <- function(S0, I0, R0, beta, gamma, t_max) {
  n_total <- S0 + I0 + R0
  cap <- max(64, 4L * ceiling(S0 + I0))
  times <- numeric(cap); Sv <- numeric(cap); Iv <- numeric(cap); Rv <- numeric(cap)
  times[1] <- 0; Sv[1] <- S0; Iv[1] <- I0; Rv[1] <- R0
  n <- 1L
  t <- 0; S <- S0; I <- I0; R <- R0

  while (t < t_max && I > 0) {
    rate_inf <- beta * S * I / n_total
    rate_rec <- gamma * I
    rate_total <- rate_inf + rate_rec
    if (rate_total <= 0) break
    t_new <- t + stats::rexp(1, rate_total)
    if (t_new > t_max) break
    t <- t_new
    if (stats::runif(1) < rate_inf / rate_total) {
      S <- S - 1; I <- I + 1
    } else {
      I <- I - 1; R <- R + 1
    }
    n <- n + 1L
    if (n > cap) {
      cap <- cap * 2L
      length(times) <- cap; length(Sv) <- cap; length(Iv) <- cap; length(Rv) <- cap
    }
    times[n] <- t; Sv[n] <- S; Iv[n] <- I; Rv[n] <- R
  }

  data.frame(time = times[seq_len(n)], S = Sv[seq_len(n)], I = Iv[seq_len(n)], R = Rv[seq_len(n)])
}

#' @keywords internal
#' @noRd
run_gillespie_seir <- function(S0, E0, I0, R0, beta, sigma, gamma, t_max) {
  n_total <- S0 + E0 + I0 + R0
  cap <- max(64, 4L * ceiling(S0 + E0 + I0))
  times <- numeric(cap); Sv <- numeric(cap); Ev <- numeric(cap); Iv <- numeric(cap); Rv <- numeric(cap)
  times[1] <- 0; Sv[1] <- S0; Ev[1] <- E0; Iv[1] <- I0; Rv[1] <- R0
  n <- 1L
  t <- 0; S <- S0; E <- E0; I <- I0; R <- R0

  while (t < t_max && (I > 0 || E > 0)) {
    rate_exp <- beta * S * I / n_total
    rate_prog <- sigma * E
    rate_rec <- gamma * I
    rate_total <- rate_exp + rate_prog + rate_rec
    if (rate_total <= 0) break
    t_new <- t + stats::rexp(1, rate_total)
    if (t_new > t_max) break
    t <- t_new
    u <- stats::runif(1) * rate_total
    if (u < rate_exp) {
      S <- S - 1; E <- E + 1
    } else if (u < rate_exp + rate_prog) {
      E <- E - 1; I <- I + 1
    } else {
      I <- I - 1; R <- R + 1
    }
    n <- n + 1L
    if (n > cap) {
      cap <- cap * 2L
      length(times) <- cap; length(Sv) <- cap; length(Ev) <- cap; length(Iv) <- cap; length(Rv) <- cap
    }
    times[n] <- t; Sv[n] <- S; Ev[n] <- E; Iv[n] <- I; Rv[n] <- R
  }

  data.frame(
    time = times[seq_len(n)], S = Sv[seq_len(n)], E = Ev[seq_len(n)],
    I = Iv[seq_len(n)], R = Rv[seq_len(n)]
  )
}

#' @export
print.gillespie_epidemic <- function(x, ...) {
  cat(sprintf("<omicsuite %s Gillespie epidemic simulation>\n\n", x$model))
  cat(sprintf("Realizations: %d\n", x$n_sim))
  cat("Initial state: ")
  cat(paste(sprintf("%s=%s", names(x$initial_state), x$initial_state), collapse = ", "))
  cat("\nParameters: ")
  cat(paste(sprintf("%s=%s", names(x$params), x$params), collapse = ", "))
  cat(sprintf("\n\nR0 = %.2f\n", x$r0))
  cat(sprintf("Median final outbreak size: %.0f (of %d)\n", stats::median(x$final_size), sum(x$initial_state)))
  cat(sprintf("Proportion of realizations with early stochastic extinction: %.1f%%\n", 100 * x$prop_extinct))
  cat("\nVerdicts:\n")
  print_verdicts(x$verdicts)
  invisible(x)
}

#' @export
summary.gillespie_epidemic <- function(object, ...) {
  list(
    r0 = object$r0,
    prop_extinct = object$prop_extinct,
    final_size = summary(object$final_size),
    peak_time = summary(object$peak_time[!object$is_extinct]),
    verdicts = object$verdicts
  )
}

#' Plot a `gillespie_epidemic` object
#'
#' @param x A `gillespie_epidemic` object.
#' @param which Character vector of plot names to display. Defaults to all
#'   plots in `x$plots`. Run `names(x$plots)` to see what's available.
#' @param ... Ignored.
#' @return Invisibly returns the list of plots shown.
#' @export
plot.gillespie_epidemic <- function(x, which = names(x$plots), ...) {
  for (nm in which) {
    print(x$plots[[nm]])
  }
  invisible(x$plots[which])
}
