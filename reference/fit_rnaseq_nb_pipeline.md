# Fit a Bayesian Hierarchical Negative Binomial Mixed Model for RNA-Seq Counts

Fits a multilevel negative binomial model to long-format RNA-seq count
data (one row per gene x sample) via
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html). The
effect of primary biological interest (`condition_var`) gets a
gene-level random slope, so genes borrow statistical strength from one
another (partial pooling) rather than being tested one at a time in
isolation. An optional subject-level random intercept accounts for
repeated-measures / longitudinal correlation within individuals.

## Usage

``` r
fit_rnaseq_nb_pipeline(
  data,
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
  ...
)
```

## Arguments

- data:

  A long-format data.frame: one row per gene x sample combination.

- gene_var:

  Character. Column identifying the gene. Gets a gene-level random
  intercept and random slope for `condition_var`.

- count_var:

  Character. Column of raw (unnormalized) RNA-seq counts.

- condition_var:

  Character. The fixed effect of primary biological interest (e.g.
  treatment arm, genotype). Given a gene-level random slope so per-gene
  effects are partially pooled toward the population average. Should be
  numeric or a two-level factor: the shrinkage summary and plot
  currently reflect only the first non-reference coefficient, so a
  factor with more than two levels will emit a warning and only show one
  contrast.

- covariates:

  Character vector of additional fixed effects to adjust for
  (population-level only – no gene-specific slope). Default `NULL`.

- subject_var:

  Optional character. Subject/individual identifier for
  repeated-measures or longitudinal designs; gets a subject-level random
  intercept. Default `NULL` (no repeated-measures structure).

- sample_var:

  Character. Column identifying the sample – used to compute the
  library-size offset when `offset_var` is not supplied.

- offset_var:

  Optional character. Column of pre-computed normalization factors (raw
  scale, e.g. total counts or a DESeq2/edgeR size factor) to use as the
  model offset. If `NULL` (the default), a simple total-count-per-sample
  offset is computed internally – adequate for a first pass, but a
  proper size-factor normalization (e.g. median-of-ratios) is
  recommended for a publication-grade analysis and can be supplied here
  instead.

- chains, iter, warmup, cores:

  Passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Defaults (`chains = 4, iter = 2000, warmup = 1000, cores = 1`) are a
  reasonable starting point for exploration; increase `cores` to
  `chains` if your machine allows running them in parallel, and increase
  `iter` if the convergence verdict flags a low effective-sample-size
  ratio.

- seed:

  Optional integer seed passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Default `NA` (brms's own default – no fixed seed). Note this must be
  `NA`, not `NULL`, since
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  cannot coerce `NULL` to a numeric value.

- cache_file:

  Optional path. If it exists, the cached fit is loaded via
  [`readRDS()`](https://rdrr.io/r/base/readRDS.html) and no model is
  refit. If it does not exist, the model is fit and then saved there via
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Recommended for
  any model that takes more than a minute or two to fit, since Stan
  model compilation and sampling both restart from scratch otherwise.

- ...:

  Additional arguments passed through to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Value

An object of class `"rnaseq_nb_pipeline"`, a list with elements:

- model:

  The fitted `brmsfit` object.

- convergence:

  A list with `max_rhat`, `min_neff_ratio`, and the full per-parameter
  summaries.

- shrinkage:

  A data.frame with one row per gene: the partially pooled (population +
  gene-level deviation) estimate of the `condition_var` effect, its
  credible interval, and the population-level estimate for comparison.

- dispersion:

  A list with the posterior summary of the negative binomial shape
  parameter.

- plots:

  A named list of `ggplot` objects: `pp_check`, `shrinkage_plot`, and
  `rhat_plot`.

- verdicts:

  A data.frame summarizing convergence, dispersion, posterior predictive
  fit, and a plain-language interpretation of the population-level
  `condition_var` effect, in the same style as
  [`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).
  A `plot` column names which entry in `plots` each row explains.

## Details

`brms` (and the Stan toolchain it depends on) is listed as a `Suggests`,
not a hard dependency of `omicsuite`, since it requires a working C++
compiler. Install it with `install.packages("brms")` before calling this
function.

## Examples

``` r
# \donttest{
if (requireNamespace("brms", quietly = TRUE)) {
  set.seed(1)
  n_genes <- 8
  n_samples <- 12
  dat <- expand.grid(gene = paste0("gene_", seq_len(n_genes)),
                      sample = paste0("sample_", seq_len(n_samples)))
  dat$condition <- factor(rep(rep(c("control", "treated"), each = n_samples / 2),
                               n_genes))
  dat$count <- stats::rnbinom(nrow(dat), mu = 50, size = 5)

  fit <- fit_rnaseq_nb_pipeline(
    data = dat, gene_var = "gene", count_var = "count",
    condition_var = "condition", sample_var = "sample",
    chains = 1, iter = 200, warmup = 100
  )
  print(fit)
}
#> Compiling Stan program...
#> Error in .fun(model_code = .x1): Boost not found; call install.packages('BH')
# }
```
