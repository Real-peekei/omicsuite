# Summarize an `rnaseq_nb_pipeline` object

Summarize an `rnaseq_nb_pipeline` object

## Usage

``` r
# S3 method for class 'rnaseq_nb_pipeline'
summary(object, ...)
```

## Arguments

- object:

  An `rnaseq_nb_pipeline` object, as returned by
  [`fit_rnaseq_nb_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_rnaseq_nb_pipeline.md).

- ...:

  Ignored.

## Value

A list with elements `model_summary` (the `brmsfit` summary),
`convergence`, `dispersion`, and `verdicts`.
