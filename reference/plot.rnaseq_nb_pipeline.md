# Plot an `rnaseq_nb_pipeline` object

Plot an `rnaseq_nb_pipeline` object

## Usage

``` r
# S3 method for class 'rnaseq_nb_pipeline'
plot(x, which = names(x$plots), ...)
```

## Arguments

- x:

  An `rnaseq_nb_pipeline` object.

- which:

  Character vector of plot names to display. Defaults to all plots in
  `x$plots`. Run `names(x$plots)` to see what's available.

- ...:

  Ignored.

## Value

Invisibly returns the list of plots shown.
