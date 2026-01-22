# Define a custom family called 'quasiproductr' to be used with the mgcv package.

The `quasiproductr` family models relationships where the link function
is the hyperbolic arctangent (`artanh`). This family defines custom
link, variance, residual, and initialization functions for use with
Generalized Additive Models (GAM).

## Usage

``` r
quasiproductr(link = "artanh")
```

## Arguments

- link:

  The link function to be used. Default is `"artanh"`.

## Value

A family object compatible with mgcv, including variance, link
functions, and residual deviance calculation.
