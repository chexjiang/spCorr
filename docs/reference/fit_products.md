# Fit product distributions to standardized gene expressions for a list of gene pairs.

This function fits product distributions to gene pairs using their
standardized expression values. It applies a Generalized Additive Model
(GAM) to each pair of genes from `gene_pair_list` to model the joint
effect while incorporating covariates.

## Usage

``` r
fit_products(
  gene_pair_list_subset,
  product_list,
  cov_mat,
  formula2,
  family2 = quasiproductr(),
  control = list(),
  ncores,
  global_test,
  critical_value = 0.05,
  return_models = FALSE,
  return_coefs = FALSE,
  return_pi = FALSE,
  preconstruct_smoother = FALSE
)
```

## Arguments

- gene_pair_list_subset:

  A data frame or matrix of gene pairs that passed filtering (e.g., via
  Moran's I test).

- product_list:

  A named list of product vectors corresponding to gene expression
  cross-products.

- cov_mat:

  A matrix or data frame containing covariates (e.g., spatial
  coordinates like `x1`, `x2`).

- formula2:

  A formula or string specifying the GAM structure for modeling spatial
  interactions (e.g., `s(x1, x2)`).

- family2:

  A distribution family to be used for model fitting. Default is
  [`quasiproductr()`](https://chexjiang.github.io/spCorr/reference/quasiproductr.md),
  which must be defined elsewhere.

- control:

  A list of control parameters passed to
  [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html).

- ncores:

  Number of cores to use for parallel processing with `mclapply`.

- global_test:

  Method for global testing in product models. Options: `"lrt"`
  (likelihood ratio test) or `"wald"` (Wald-style smooth term test).
  Default is `"wald"`.

- return_models:

  Logical; if `TRUE`, returns the full model object for each gene pair.

- return_coefs:

  Logical; if `TRUE`, returns model coefficients and variance-covariance
  matrices.

- return_pi:

  Logical; if `TRUE`, returns predicted interval for fitted correlation.

- preconstruct_smoother:

  Logical; if `TRUE`, modifies the smoother basis (e.g., `'tp'` to
  `'tpcached'`) for caching and speed optimization. Default is `FALSE`.

## Value

A list where each element corresponds to a gene pair. The contents
depend on `return_models` and `return_coefs`:

- res_global:

  Result of the global test (either p-value or list with EDF and
  p-value).

- fitted_rho:

  Fitted values from the GAM.

- model:

  (Optional) Full [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html)
  object, if `return_models = TRUE`.

- model_coef:

  (Optional) List with `beta` (coefficients) and `beta_cov`
  (variance-covariance matrix), if `return_coefs = TRUE`.

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

# Compute gene pair products and optionally subset based on Moran's I
check_result <- check_products(
  gene_pair_list = test_data$gene_pair_list,
  marginals = marginal_res$marginal,
  cov_mat = test_data$cov_mat,
  check_morani = FALSE,
  ncores = 2
)

# Fit product distributions to gene_pair_list_subset
model_list <- fit_products(
  gene_pair_list_subset = check_result$gene_pair_list_subset,
  product_list = check_result$product_list,
  cov_mat = test_data$cov_mat,
  formula2 = "s(x1, x2, bs='tp', k=50)",
  family2 = quasiproductr(),
  control = list(),
  ncores = 2,
  global_test = "wald",
  return_models = FALSE,
  return_coefs = FALSE,
  preconstruct_smoother = TRUE
)
#> Warning: all scheduled cores encountered errors in user code
```
