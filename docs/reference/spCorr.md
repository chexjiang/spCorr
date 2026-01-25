# The wrapper for the whole spCorr pipeline

This function fits conditional margins and models local correlation for
a given list of genes and gene pairs using GAM-based models. It also
performs statistical testing to identify significant patterns in gene
co-expression.

## Usage

``` r
spCorr(
  count_mat,
  gene_list,
  gene_pair_list,
  cov_mat,
  formula1 = "1",
  family1 = "nb",
  formula2 = "s(x1, x2, bs='tp', k=50)",
  family2 = quasiproductr(),
  DT = TRUE,
  global_test = "lrt",
  return_models = FALSE,
  return_coefs = FALSE,
  return_pi = FALSE,
  check_morani = FALSE,
  preconstruct_smoother = TRUE,
  ncores = 2,
  control = list(),
  epsilon = 1e-06,
  seed = 123
)
```

## Arguments

- count_mat:

  A matrix of raw gene expression counts (genes by spots/cells).

- gene_list:

  A vector of gene names or row indices for which marginals will be fit.

- gene_pair_list:

  A two-column data frame or matrix specifying gene pairs (by name or
  index).

- cov_mat:

  A data frame of covariates used in both marginal and product fitting
  (must contain `x1` and `x2` for spatial coordinates).

- formula1:

  Formula or string specifying the model structure for marginals (e.g.,
  `"~ covariate"`). Use `"1"` for intercept-only.

- family1:

  Distribution family for marginal models. Options: `"gaussian"`,
  `"poisson"`, `"nb"`, or `"zinb"`. Default: `"nb"`.

- formula2:

  Formula or string specifying the smoother for GAMs (e.g.,
  `"s(x1, x2, bs='tp', k=50)"`).

- family2:

  A GAM family object for product modeling (e.g.,
  [`quasiproductr()`](https://chexjiang.github.io/spCorr/reference/quasiproductr.md)).

- DT:

  Logical; if `TRUE`, applies a discrete transformation to marginals.
  Default is `TRUE`.

- global_test:

  Method for global testing in product models. Options: `"lrt"`
  (likelihood ratio test) or `"wald"` (Wald-style smooth term test).
  Default is `"wald"`.

- return_models:

  Logical; if `TRUE`, returns full GAM model objects for each gene pair.
  Default is `FALSE`.

- return_coefs:

  Logical; if `TRUE`, returns model coefficients and variances. Default
  is `FALSE`.

- check_morani:

  Logical; if `TRUE`, filters gene pairs using Moran's I on the product.
  Default is `FALSE`.

- preconstruct_smoother:

  Logical; if `TRUE`, replaces `bs='tp'`/`'gp'` with
  `tpcached`/`gpcached` for faster computation. Default is `TRUE`.

- ncores:

  Integer number of cores for parallel processing. Default is `2`.

- control:

  A list of control parameters passed to
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) during product
  fitting.

- epsilon:

  A small constant to avoid boundary issues in the uniform-to-Gaussian
  transformation. Default is `1e-6`.

- seed:

  Random seed for reproducibility. Default is `123`.

## Value

A named list containing:

- res_global:

  A vector of adjusted p-values (FDR) from global tests for each gene
  pair.

- res_local:

  A matrix of local fitted values (spatial correlation estimates) for
  each pair across spatial spots.

- marginals:

  A matrix of standardized marginal values (standard normal) for each
  gene.

- residuals:

  A matrix of uniform-transformed residuals for each gene.

- model_list:

  (Optional) List of fitted GAM models if `return_models = TRUE`.

- model_coef_list:

  (Optional) List of model coefficients if `return_coefs = TRUE`.

## Details

The pipeline consists of:

1.  Fitting conditional marginal distributions to individual genes.

2.  Calculating pairwise product expressions and optionally filtering
    via Moran's I.

3.  Fitting GAMs to model local spatial correlations between gene pairs.

4.  Outputting statistical results and fitted values (or models).

## See also

[`fit_marginals()`](https://chexjiang.github.io/spCorr/reference/fit_marginals.md),
[`check_products()`](https://chexjiang.github.io/spCorr/reference/check_products.md),
[`fit_products()`](https://chexjiang.github.io/spCorr/reference/fit_products.md)

## Examples

``` r
data(test_data)
result <- spCorr(
  count_mat = test_data$count_mat,
  gene_list = test_data$gene_list,
  gene_pair_list = test_data$gene_pair_list,
  cov_mat = test_data$cov_mat,
  formula1 = "layer_annotations",
  family1 = "nb",
  formula2 = "s(x1, x2, bs='tp', k=50)",
  family2 = quasiproductr(),
  DT = TRUE,
  global_test = "lrt",
  return_models = FALSE,
  return_coefs = FALSE,
  check_morani = FALSE,
  preconstruct_smoother = TRUE,
  ncores = 2,
  control = list(),
  epsilon = 1e-6,
  seed = 123
)
#> Start Marginal Fitting for 20 genes
#> Start Cross-Product Fitting for 10 gene pairs
#> Warning: all scheduled cores encountered errors in user code
```
