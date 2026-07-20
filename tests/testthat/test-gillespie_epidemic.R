test_that("simulate_gillespie_epidemic (SIR) returns the expected structure", {
  fit <- simulate_gillespie_epidemic(
    model = "SIR",
    initial_state = c(S = 199, I = 1, R = 0),
    params = list(beta = 0.5, gamma = 0.1),
    t_max = 60, n_sim = 15, seed = 1
  )

  expect_s3_class(fit, "gillespie_epidemic")
  expect_identical(fit$model, "SIR")
  expect_length(fit$simulations, 15)
  expect_true(all(c("S", "I", "R") %in% names(fit$simulations[[1]])))
  expect_true(all(c("S", "I", "R") %in% names(fit$summary)))
  expect_equal(fit$r0, 5)
  expect_true(is.data.frame(fit$verdicts))
  expect_true(all(c("basic_reproduction_number", "stochastic_extinction") %in% fit$verdicts$check))
  expect_true(all(c("trajectory_plot", "final_size_hist", "peak_time_hist") %in% names(fit$plots)))
})

test_that("simulate_gillespie_epidemic (SEIR) returns the expected structure", {
  fit <- simulate_gillespie_epidemic(
    model = "SEIR",
    initial_state = c(S = 199, E = 0, I = 1, R = 0),
    params = list(beta = 0.5, sigma = 0.3, gamma = 0.1),
    t_max = 80, n_sim = 15, seed = 2
  )

  expect_identical(fit$model, "SEIR")
  expect_true(all(c("S", "E", "I", "R") %in% names(fit$simulations[[1]])))
  expect_true(all(c("S", "E", "I", "R") %in% names(fit$summary)))
})

test_that("each realization conserves total population at every recorded time point", {
  fit <- simulate_gillespie_epidemic(
    model = "SIR",
    initial_state = c(S = 99, I = 1, R = 0),
    params = list(beta = 0.6, gamma = 0.15),
    t_max = 100, n_sim = 10, seed = 3
  )
  n_total <- sum(fit$initial_state)
  for (df in fit$simulations) {
    expect_true(all(abs(rowSums(df[, c("S", "I", "R")]) - n_total) < 1e-8))
    expect_true(all(df$S >= 0 & df$I >= 0 & df$R >= 0))
  }
})

test_that("R0 below the epidemic threshold produces predominantly small outbreaks", {
  # beta < gamma => R0 < 1; with a single initial infective in a modest
  # population, most realizations should fade out quickly rather than
  # infecting a large share of the population.
  fit <- simulate_gillespie_epidemic(
    model = "SIR",
    initial_state = c(S = 199, I = 1, R = 0),
    params = list(beta = 0.05, gamma = 0.2),
    t_max = 200, n_sim = 100, seed = 4
  )
  expect_true(fit$r0 < 1)
  expect_true(fit$prop_extinct > 0.5)
})

test_that("simulate_gillespie_epidemic errors informatively on missing state/params", {
  expect_error(
    simulate_gillespie_epidemic(
      model = "SIR", initial_state = c(S = 99, I = 1),
      params = list(beta = 0.3, gamma = 0.1), t_max = 50
    ),
    "initial_state is missing"
  )
  expect_error(
    simulate_gillespie_epidemic(
      model = "SEIR", initial_state = c(S = 99, E = 0, I = 1, R = 0),
      params = list(beta = 0.3, gamma = 0.1), t_max = 50
    ),
    "params is missing"
  )
})

test_that("print.gillespie_epidemic and plot.gillespie_epidemic run without error", {
  fit <- simulate_gillespie_epidemic(
    model = "SIR",
    initial_state = c(S = 99, I = 1, R = 0),
    params = list(beta = 0.4, gamma = 0.1),
    t_max = 60, n_sim = 10, seed = 5
  )
  expect_output(print(fit), "omicsuite SIR Gillespie")
  grDevices::pdf(NULL)
  plots_shown <- plot(fit, which = "trajectory_plot")
  grDevices::dev.off()
  expect_true("trajectory_plot" %in% names(plots_shown))
})
