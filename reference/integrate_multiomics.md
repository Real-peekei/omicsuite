# Integrate Multi-Omics Data Blocks via Regularized Generalized Canonical Correlation Analysis

Fits a regularized generalized canonical correlation analysis (RGCCA,
via
[`RGCCA::rgcca()`](https://rgcca-factory.github.io/RGCCA/reference/rgcca.html))
across two or more omics blocks measured on the same samples – e.g.
transcriptomics, proteomics, and methylation for the same patients.
Unlike stacking blocks into one matrix and running PCA, RGCCA finds
components that maximize covariation *between* blocks, so the result
reflects signal shared across omics layers rather than whichever layer
happens to have the most features.

## Usage

``` r
integrate_multiomics(
  blocks,
  group = NULL,
  ncomp = 2,
  scheme = c("factorial", "centroid", "horst"),
  scale = TRUE,
  tau = "optimal",
  n_boot = 50,
  boot_top_n = 20,
  alpha_group_separation = 0.05,
  seed = NULL
)
```

## Arguments

- blocks:

  A named list of matrices or data.frames, one per omics layer, samples
  in rows and features in columns. Row names are used to align samples
  across blocks; if blocks don't share identical sample sets, only the
  common samples are used (with a note in `verdicts`).

- group:

  Optional vector of group/phenotype labels for coloring the
  sample-score plot (e.g. a cancer subtype). Either a named vector keyed
  by sample ID, or an unnamed vector in the same order as the common
  samples – if unnamed, order is assumed to match and is not verified.
  Not used in the RGCCA fit itself, only for plotting.

- ncomp:

  Integer. Number of components per block. Default `2`.

- scheme:

  Character. RGCCA connection scheme: `"factorial"`, `"centroid"`, or
  `"horst"`. Default `"factorial"`.

- scale:

  Logical. Standardize each block before fitting. Default `TRUE`.

- tau:

  Shrinkage parameter passed to
  [`RGCCA::rgcca()`](https://rgcca-factory.github.io/RGCCA/reference/rgcca.html).
  Default `"optimal"`, which uses Schafer-Strimmer analytical shrinkage
  estimated separately per block.

- n_boot:

  Integer. Number of case-resampling bootstrap replicates for the
  loading stability check. Default `50`. Each replicate refits the full
  RGCCA model, so this is the slowest part of the pipeline – reduce it
  for exploration and increase it before reporting results.

- boot_top_n:

  Integer. Per block, only the top `boot_top_n` features by absolute
  component-1 loading are tracked for stability (tracking every feature
  in a high-dimensional block would make the bootstrap prohibitively
  slow without adding much information). Default `20`.

- alpha_group_separation:

  Significance threshold for the group-separation interpretation (a
  Kruskal-Wallis test of each block's component 1 by `group`), only used
  when `group` is supplied. Default `0.05`.

- seed:

  Optional integer seed.

## Value

An object of class `"multiomics_pipeline"`, a list with elements:

- model:

  The fitted `rgcca` object from
  [`RGCCA::rgcca()`](https://rgcca-factory.github.io/RGCCA/reference/rgcca.html).

- scores:

  A named list of data.frames, one per block: sample scores on each
  component, with `sample_id` and `group` (if supplied) columns.

- variance_explained:

  A data.frame: average variance explained (AVE) per block per
  component.

- stability:

  A data.frame: per-block, per-feature bootstrap sign-agreement rate for
  the top-loading features on component 1.

- plots:

  A named list of `ggplot` objects: `block_scores`,
  `variance_explained`, `stability`.

- verdicts:

  A data.frame summarizing sample alignment, variance explained, and
  loading stability, plus a per-block group-separation interpretation
  when `group` is supplied, in the same style as
  [`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).
  A `plot` column names which entry in `plots` each row explains.

## Details

The loading stability check here is a case-resampling bootstrap built on
top of
[`RGCCA::rgcca()`](https://rgcca-factory.github.io/RGCCA/reference/rgcca.html)
rather than a call to
[`RGCCA::rgcca_bootstrap()`](https://rgcca-factory.github.io/RGCCA/reference/rgcca_bootstrap.html),
so the stability metric (sign-agreement rate for each block's
top-loading features) is transparent and doesn't depend on the exact
internals of any one RGCCA version.

## Examples

``` r
# \donttest{
set.seed(1)
n <- 40
block1 <- matrix(rnorm(n * 15), nrow = n, dimnames = list(paste0("s", 1:n), NULL))
block2 <- matrix(rnorm(n * 10), nrow = n, dimnames = list(paste0("s", 1:n), NULL))
# give the blocks some shared signal
shared <- rnorm(n)
block1[, 1] <- block1[, 1] + shared
block2[, 1] <- block2[, 1] + shared

fit <- integrate_multiomics(
  blocks = list(omics_a = block1, omics_b = block2),
  ncomp = 2, n_boot = 20
)
print(fit)
#> <omicsuite multi-omics RGCCA integration pipeline>
#> 
#> Blocks: omics_a, omics_b
#> 
#> Diagnostic verdicts:
#> [OK  ] sample_alignment             All 40 samples were present in every block.
#> [INFO] variance_explained[omics_a]  Component 1 explains 13.0% of variance within this block. This block is contributing a meaningful share of its own variance to the shared components. [see plots$variance_explained]
#> [INFO] variance_explained[omics_b]  Component 1 explains 19.1% of variance within this block. This block is contributing a meaningful share of its own variance to the shared components. [see plots$variance_explained]
#> [FLAG] loading_stability[omics_a]   Top 15 loadings had a mean bootstrap sign-agreement rate of 69.7% across 20 successful resamples. A non-trivial share of top loadings flip sign under resampling -- treat the specific feature ranking within this block cautiously, even if the block-level integration is otherwise sound. [see plots$stability]
#> [FLAG] loading_stability[omics_b]   Top 10 loadings had a mean bootstrap sign-agreement rate of 73.5% across 20 successful resamples. A non-trivial share of top loadings flip sign under resampling -- treat the specific feature ranking within this block cautiously, even if the block-level integration is otherwise sound. [see plots$stability]
# }
```
