# Summarize a `gillespie_epidemic` object

Summarize a `gillespie_epidemic` object

## Usage

``` r
# S3 method for class 'gillespie_epidemic'
summary(object, ...)
```

## Arguments

- object:

  A `gillespie_epidemic` object, as returned by
  [`simulate_gillespie_epidemic()`](https://real-peekei.github.io/omicsuite/reference/simulate_gillespie_epidemic.md).

- ...:

  Ignored.

## Value

A list with elements `r0`, `prop_extinct`, `final_size` (a
[`summary()`](https://rdrr.io/r/base/summary.html) of the final-size
distribution across realizations), `peak_time` (likewise, restricted to
non-extinct realizations), and `verdicts`.
