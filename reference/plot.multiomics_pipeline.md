# Plot a `multiomics_pipeline` object

Plot a `multiomics_pipeline` object

## Usage

``` r
# S3 method for class 'multiomics_pipeline'
plot(x, which = names(x$plots), ...)
```

## Arguments

- x:

  A `multiomics_pipeline` object.

- which:

  Character vector of plot names to display. Defaults to all plots in
  `x$plots`. Run `names(x$plots)` to see what's available.

- ...:

  Ignored.

## Value

Invisibly returns the list of plots shown.
