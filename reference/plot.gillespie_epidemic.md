# Plot a `gillespie_epidemic` object

Plot a `gillespie_epidemic` object

## Usage

``` r
# S3 method for class 'gillespie_epidemic'
plot(x, which = names(x$plots), ...)
```

## Arguments

- x:

  A `gillespie_epidemic` object.

- which:

  Character vector of plot names to display. Defaults to all plots in
  `x$plots`. Run `names(x$plots)` to see what's available.

- ...:

  Ignored.

## Value

Invisibly returns the list of plots shown.
