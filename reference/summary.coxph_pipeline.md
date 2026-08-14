# Summarize a `coxph_pipeline` object

Summarize a `coxph_pipeline` object

## Usage

``` r
# S3 method for class 'coxph_pipeline'
summary(object, ...)
```

## Arguments

- object:

  A `coxph_pipeline` object, as returned by
  [`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).

- ...:

  Ignored.

## Value

A list with elements `unadjusted` and `adjusted` (each a `summary.coxph`
object, `adjusted` being `NULL` if no adjusted model was fit),
`ph_test`, and `verdicts`.
