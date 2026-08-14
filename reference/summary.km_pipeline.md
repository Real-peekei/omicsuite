# Summarize a `km_pipeline` object

Summarize a `km_pipeline` object

## Usage

``` r
# S3 method for class 'km_pipeline'
summary(object, ...)
```

## Arguments

- object:

  A `km_pipeline` object, as returned by
  [`fit_km_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_km_pipeline.md).

- ...:

  Ignored.

## Value

A list with elements `median_survival`, `logrank_test`, and `verdicts`.
