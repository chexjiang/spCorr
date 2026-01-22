# Fit conditional margins to standard Gaussian distributions for a list of genes.

This function applies marginal fitting and transforms each gene's
expression to follow a standard Gaussian distribution.

## Usage

``` r
fit_marginals(
  gene_list,
  count_mat,
  cov_mat,
  formula1,
  family1,
  DT = TRUE,
  epsilon = 1e-06,
  ncores = ncores,
  seed = 123
)
```

## Arguments

- gene_list:

  A vector of gene names or indices (row names or row numbers of
  `count_mat`) to process.

- count_mat:

  A matrix of raw gene expression counts. Rows correspond to genes,
  columns to observations (cells/spots).

- cov_mat:

  A matrix or data frame of covariates used for marginal modeling (e.g.,
  spatial coordinates or experimental annotations).

- formula1:

  A formula object or string (e.g., `"~ covariate1 + covariate2"`)
  specifying the model structure for the mean.

- family1:

  A string specifying the distribution family to be used for modeling.
  Supported values include `"gaussian"`, `"poisson"`, `"nb"`, or
  `"zinb"`.

- DT:

  Logical; if `TRUE`, applies a discrete transformation suitable for
  count data. Default is `TRUE`.

- epsilon:

  A small numeric constant to avoid boundary issues (e.g., `0` or `1`
  values in uniform distribution). Default is `1e-6`.

- ncores:

  Integer specifying the number of cores to use for parallel processing
  via [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html).

- seed:

  Random seed for reproducibility. Default is `123`.

## Value

A list containing two matrices:

- marginal:

  A matrix of transformed values from each gene, transformed to follow a
  standard normal distribution.

- residual:

  A matrix of values transformed to standard uniform distribution
  (before applying the Gaussian quantile function).

## Examples

``` r
data(test_data)
# Fit standardized marginals for gene expressions
marginal_res <- fit_marginals(
  gene_list = test_data$gene_list,
  count_mat = test_data$count_mat,
  cov_mat = test_data$cov_mat,
  formula1 = "layer_annotations",
  family1 = "nb",
  DT = TRUE,
  epsilon = 1e-6,
  ncores = 2
)
```
