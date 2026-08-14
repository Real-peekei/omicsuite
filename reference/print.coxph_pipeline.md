# Print a `coxph_pipeline` object

Prints the unadjusted (and adjusted, if fit) Cox model summaries
followed by the full verdict table.

## Usage

``` r
# S3 method for class 'coxph_pipeline'
print(x, ...)
```

## Arguments

- x:

  A `coxph_pipeline` object, as returned by
  [`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md).

- ...:

  Ignored.

## Value

Invisibly returns `x`.
