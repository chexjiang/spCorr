#' The wrapper for the whole spCorr pipeline
#'
#' This function fits conditional margins and models local correlation for a given list of genes and gene pairs using GAM-based models.
#' It also performs statistical testing to identify significant patterns in gene co-expression.
#'
#' The pipeline consists of:
#' \enumerate{
#'   \item Fitting conditional marginal distributions to individual genes.
#'   \item Calculating pairwise product expressions and optionally filtering via Moran's I.
#'   \item Fitting GAMs to model local spatial correlations between gene pairs.
#'   \item Outputting statistical results and fitted values (or models).
#' }
#'
#' @param count_mat A matrix of raw gene expression counts (genes × spots/cells).
#' @param gene_list A vector of gene names or row indices for which marginals will be fit.
#' @param gene_pair_list A two-column data frame or matrix specifying gene pairs (by name or index).
#' @param cov_mat A data frame of covariates used in both marginal and product fitting (must contain `x1` and `x2` for spatial coordinates).
#' @param formula1 Formula or string specifying the model structure for marginals (e.g., `"~ covariate"`). Use `"1"` for intercept-only.
#' @param family1 Distribution family for marginal models. Options: `"gaussian"`, `"poisson"`, `"nb"`, or `"zinb"`. Default: `"nb"`.
#' @param formula2 Formula or string specifying the smoother for GAMs (e.g., `"s(x1, x2, bs='tp', k=50)"`).
#' @param family2 A GAM family object for product modeling (e.g., `quasiproductr()`).
#' @param DT Logical; if `TRUE`, applies a discrete transformation to marginals. Default is `TRUE`.
#' @param global_test Method for global testing in product models. Options: `"lrt"` (likelihood ratio test) or `"wald"` (Wald-style smooth term test). Default is `"wald"`.
#' @param return_models Logical; if `TRUE`, returns full GAM model objects for each gene pair. Default is `FALSE`.
#' @param return_coefs Logical; if `TRUE`, returns model coefficients and variances. Default is `FALSE`.
#' @param check_morani Logical; if `TRUE`, filters gene pairs using Moran's I on the product. Default is `FALSE`.
#' @param preconstruct_smoother Logical; if `TRUE`, replaces `bs='tp'`/`'gp'` with `tpcached`/`gpcached` for faster computation. Default is `TRUE`.
#' @param ncores Integer number of cores for parallel processing. Default is `2`.
#' @param control A list of control parameters passed to `mgcv::gam()` during product fitting.
#' @param epsilon A small constant to avoid boundary issues in the uniform-to-Gaussian transformation. Default is `1e-6`.
#' @param seed Random seed for reproducibility. Default is `123`.
#'
#' @return A named list containing:
#' \describe{
#'   \item{res_global}{A vector of adjusted p-values (FDR) from global tests for each gene pair.}
#'   \item{res_local}{A matrix of local fitted values (spatial correlation estimates) for each pair across spatial spots.}
#'   \item{marginals}{A matrix of standardized marginal values (standard normal) for each gene.}
#'   \item{residuals}{A matrix of uniform-transformed residuals for each gene.}
#'   \item{model_list}{(Optional) List of fitted GAM models if `return_models = TRUE`.}
#'   \item{model_coef_list}{(Optional) List of model coefficients if `return_coefs = TRUE`.}
#' }
#'
#' @examples
#' data(test_data)
#' result <- spCorr(
#'   count_mat = test_data$count_mat,
#'   gene_list = test_data$gene_list,
#'   gene_pair_list = test_data$gene_pair_list,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
#'   family1 = "nb",
#'   formula2 = "s(x1, x2, bs='tp', k=50)",
#'   family2 = quasiproductr(),
#'   DT = TRUE,
#'   global_test = "LRT",
#'   return_models = FALSE,
#'   return_coefs = FALSE,
#'   check_morani = FALSE,
#'   preconstruct_smoother = TRUE,
#'   ncores = 2,
#'   control = list(),
#'   epsilon = 1e-6,
#'   seed = 123
#' )
#' @seealso [fit_marginals()], [check_products()], [fit_products()]
#' @importFrom mgcv gam
#' @importFrom parallel mclapply
#' @export
spCorr <- function(count_mat,
                   gene_list,
                   gene_pair_list,
                   cov_mat,
                   formula1 = "1",
                   family1 = "nb",
                   formula2 = "s(x1, x2, bs='tp', k=50)",
                   family2 = quasiproductr(),
                   DT = TRUE,
                   global_test = "wald",
                   return_models = FALSE,
                   return_coefs = FALSE,
                   check_morani = FALSE,
                   preconstruct_smoother = TRUE,
                   ncores = 2,
                   control = list(),
                   epsilon = 1e-6,
                   seed = 123) {
  # Set reproducibility seed
  set.seed(seed)

  # Create the caching environment
  .smoother_env <- new.env(parent = emptyenv())
  # Temporarily assign to global environment
  assign(".smoother_env", .smoother_env, envir = .GlobalEnv)

  on.exit(
    {
      # Clean up .smoother_env after spCorr finishes
      rm(".smoother_env", envir = .GlobalEnv)
    },
    add = TRUE
  )


  ## Fit conditional margins to gene_list
  message("Start Marginal Fitting for ", length(gene_list), " genes")
  marginal_res <- fit_marginals(
    gene_list = gene_list,
    count_mat = count_mat,
    cov_mat = cov_mat,
    formula1 = formula1,
    family = family1,
    DT = DT,
    ncores = ncores
  )
  marginals <- marginal_res$marginal
  residuals <- marginal_res$residual


  ## Check and subset spatially varying cross product
  message("Start Extracting Spatially Varying Gene Pairs")
  check_product_res <- check_products(
    gene_pair_list = gene_pair_list,
    marginals = marginals,
    cov_mat = cov_mat,
    check_morani = check_morani,
    ncores = ncores
  )
  product_list <- check_product_res$product_list
  gene_pair_list_subset <- check_product_res$gene_pair_list_subset



  ## Fit product distributions to gene_pair_list_subset
  message("Start Product Fitting for ", nrow(gene_pair_list), " gene pairs")
  product_res_list <- fit_products(
    gene_pair_list_subset = gene_pair_list_subset,
    product_list = product_list,
    cov_mat = cov_mat,
    formula2 = formula2,
    family2 = quasiproductr(),
    control = control,
    ncores = ncores,
    global_test = global_test,
    return_models = return_models,
    return_coefs = return_coefs,
    preconstruct_smoother = preconstruct_smoother
  )

  ## Extract the global testing result
  res_global <- do.call(rbind, lapply(product_res_list, function(x) {
    tryCatch(x$res_global, error = function(e) NULL)
  }))
  if (!is.null(res_global) && length(res_global) > 0) {
    res_global <- p.adjust(res_global, method = "fdr")
  }

  ## Extract local fitted values
  res_local <- do.call(rbind, lapply(product_res_list, function(x) {
    tryCatch(x$fitted_rho, error = function(e) NULL)
  }))
  colnames(res_local) <- row.names(cov_mat)


  ## Default is only return fitted values
  if (return_models) {
    # Model
    model_list <- lapply(product_res_list, function(x) x$model)
    return(list(
      res_global = res_global,
      res_local = res_local,
      marginals = marginals,
      residuals = residuals,
      model_list = model_list
    ))
  } else if (return_coefs) {
    # Model coef
    model_coef_list <- lapply(product_res_list, function(x) {
      tryCatch(x$model_coef, error = function(e) NULL)
    })

    return(list(
      res_global = res_global,
      marginals = marginals,
      residuals = residuals,
      model_coef_list = model_coef_list
    ))
  } else {
    # Fitted values
    return(list(
      res_global = res_global,
      res_local = res_local,
      residuals = residuals,
      marginals = marginals
    ))
  }
}
