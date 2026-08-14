# Plot a `coxph_pipeline` object

Displays the diagnostic plots produced by
[`fit_coxph_pipeline()`](https://real-peekei.github.io/omicsuite/reference/fit_coxph_pipeline.md)
one at a time. Use `which` to select a subset instead of stepping
through all of them.

## Usage

``` r
# S3 method for class 'coxph_pipeline'
plot(x, which = names(x$plots), ...)
```

## Arguments

- x:

  A `coxph_pipeline` object.

- which:

  Character vector of plot names to display. Defaults to all plots in
  `x$plots`. Run `names(x$plots)` to see what's available.

- ...:

  Ignored.

## Value

Invisibly returns the list of plots shown.
