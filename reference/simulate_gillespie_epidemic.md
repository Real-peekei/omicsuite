# Simulate a Stochastic SIR or SEIR Epidemic via the Gillespie Algorithm

Runs one or more exact stochastic realizations of a compartmental
epidemic model (SIR or SEIR) using Gillespie's direct method (the
stochastic simulation algorithm), then summarizes the ensemble: a median
trajectory with a 5-95% envelope, the distribution of final outbreak
size and peak timing across realizations, the basic reproduction number
implied by the rate parameters, and the proportion of realizations that
fade out early by chance rather than becoming a sustained outbreak – a
distinction a deterministic ODE model can't show you at all.

## Usage

``` r
simulate_gillespie_epidemic(
  model = c("SIR", "SEIR"),
  initial_state,
  params,
  t_max,
  n_sim = 100,
  seed = NULL,
  grid_points = 200,
  extinction_threshold = 0.05
)
```

## Arguments

- model:

  Character. Either `"SIR"` or `"SEIR"`.

- initial_state:

  Named numeric vector of starting compartment counts. Must include `S`,
  `I`, `R` for `"SIR"`, or `S`, `E`, `I`, `R` for `"SEIR"`.

- params:

  Named list of rate parameters. Must include `beta` (transmission rate)
  and `gamma` (recovery rate) for `"SIR"`; `"SEIR"` additionally
  requires `sigma` (the rate of progression from exposed to infectious,
  i.e. `1 / sigma` is the mean incubation period).

- t_max:

  Numeric. Maximum simulation time.

- n_sim:

  Integer. Number of independent stochastic realizations to run. Default
  `100`, enough to characterize the extinction probability and the
  envelope around the median trajectory without being slow.

- seed:

  Optional integer seed for reproducibility.

- grid_points:

  Integer. Number of time points used when interpolating realizations
  onto a common grid for the summary envelope and plots. Default `200`.

- extinction_threshold:

  Numeric in `(0, 1)`. A realization is classified as an early
  stochastic fadeout if its final outbreak size (total ever infected) is
  below this fraction of the total population. Default `0.05`.

## Value

An object of class `"gillespie_epidemic"`, a list with elements:

- simulations:

  A list of `n_sim` data.frames, one per realization, each with columns
  `time` and one column per compartment.

- summary:

  A named list (one entry per compartment) of data.frames with `time`,
  `median`, `lower` (5th percentile), and `upper` (95th percentile)
  across realizations, on a common time grid.

- final_size, peak_time, peak_size:

  Numeric vectors, one entry per realization: total ever infected, time
  of peak infectious count, and peak infectious count.

- r0:

  The basic reproduction number implied by `params` (`beta / gamma`, the
  same threshold quantity for both SIR and SEIR).

- prop_extinct:

  Proportion of realizations classified as an early stochastic fadeout.

- plots:

  A named list of `ggplot` objects: `trajectory_plot`,
  `final_size_hist`, `peak_time_hist`.

- verdicts:

  A data.frame with plain-language interpretation rows for the basic
  reproduction number, the stochastic extinction probability, and the
  final-size distribution among sustained outbreaks, in the same style
  as
  [`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).
  A `plot` column names which entry in `plots` each row explains.

## Examples

``` r
sir_fit <- simulate_gillespie_epidemic(
  model = "SIR",
  initial_state = c(S = 999, I = 1, R = 0),
  params = list(beta = 0.4, gamma = 0.1),
  t_max = 100, n_sim = 50, seed = 1
)
print(sir_fit)
#> <omicsuite SIR Gillespie epidemic simulation>
#> 
#> Realizations: 50
#> Initial state: S=999, I=1, R=0
#> Parameters: beta=0.4, gamma=0.1
#> 
#> R0 = 4.00
#> Median final outbreak size: 978 (of 1000)
#> Proportion of realizations with early stochastic extinction: 20.0%
#> 
#> Verdicts:
#> [INTP] interpretation[R0]           R0 = beta / gamma = 4.00. Above the epidemic threshold -- sustained transmission is expected in a well-mixed population. [see plots$trajectory_plot]
#> [INTP] interpretation[extinction_probability] 20.0% of 50 realizations resulted in early stochastic fadeout (fewer than 5% of the population ever infected). Most realizations produced a sustained outbreak; stochastic fadeout is a minor consideration here. [see plots$peak_time_hist]
#> [INTP] interpretation[final_size]   Among the 40 realization(s) that became a sustained outbreak (not an early fadeout), the median final size was 980 (98.0% of the population of 1000), with an interquartile range of 977 to 984. [see plots$final_size_hist]

seir_fit <- simulate_gillespie_epidemic(
  model = "SEIR",
  initial_state = c(S = 999, E = 0, I = 1, R = 0),
  params = list(beta = 0.4, sigma = 0.2, gamma = 0.1),
  t_max = 150, n_sim = 50, seed = 1
)
print(seir_fit)
#> <omicsuite SEIR Gillespie epidemic simulation>
#> 
#> Realizations: 50
#> Initial state: S=999, E=0, I=1, R=0
#> Parameters: beta=0.4, sigma=0.2, gamma=0.1
#> 
#> R0 = 4.00
#> Median final outbreak size: 980 (of 1000)
#> Proportion of realizations with early stochastic extinction: 18.0%
#> 
#> Verdicts:
#> [INTP] interpretation[R0]           R0 = beta / gamma = 4.00. Above the epidemic threshold -- sustained transmission is expected in a well-mixed population. [see plots$trajectory_plot]
#> [INTP] interpretation[extinction_probability] 18.0% of 50 realizations resulted in early stochastic fadeout (fewer than 5% of the population ever infected). Most realizations produced a sustained outbreak; stochastic fadeout is a minor consideration here. [see plots$peak_time_hist]
#> [INTP] interpretation[final_size]   Among the 41 realization(s) that became a sustained outbreak (not an early fadeout), the median final size was 981 (98.1% of the population of 1000), with an interquartile range of 978 to 986. [see plots$final_size_hist]
```
