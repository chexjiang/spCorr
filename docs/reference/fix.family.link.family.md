# Redefine fix.family.link.family to add derivatives to a family object

This function adds the second, third, and fourth derivatives of the link
function with respect to mu to the family object, which is used for
Newton-like optimization.

## Usage

``` r
# S3 method for class 'family'
fix.family.link(fam)
```

## Arguments

- fam:

  A family object for which to add derivatives.

## Value

A modified family object with additional components for derivatives.
