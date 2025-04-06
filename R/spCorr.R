#' The wrapper for the whole spCorr pipeline
#'
#' This function fits conditional margins and models local correlation for a given list of genes and gene pairs using GAM-based models.
#' It also performs statistical testing to identify significant patterns in gene co-expression.
#'
#' @param count_mat A matrix of counts where rows represent genes and columns represent observations.
#' @param gene_list A vector of gene names or indices for which the conditional margins are to be fit.
#' @param gene_pair_list A data frame or matrix containing pairs of gene names (or indices) to be analyzed.
#' @param cov_mat A matrix or data frame of covariates used for fitting.
#' @param formula1 A formula object or string specifying the model for fitting the marginal distributions (e.g., `~ covariate`).
#' @param family1 The distribution family for marginal fitting. Options include `'gaussian'`, `'poisson'`, or `'nb'`. Default is `'nb'`.
#' @param formula2 A formula object or string specifying the model for fitting the product distributions (e.g., `s(x1, x2, bs='tp', k=50)`).
#' @param family2 The distribution family for product fitting. Default is `quasiproductr()`.
#' @param DT Logical. If `TRUE`, applies discrete transformation during margin fitting. Default is `TRUE`.
#' @param return_models Logical. If `TRUE`, returns the fitted model objects along with results. Default is `FALSE`.
#' @param ncores Integer. The number of cores to use for parallel processing. Default is `2`.
#' @param control A list of control parameters passed to the fitting functions.
#' @param seed Integer. Seed value for reproducibility. Default is `123`.
#' @param local_testing Logical. If `TRUE`, performs local testing for each gene pair. Default is `FALSE`.
#' @param preconstruct_smoother Logical. If `TRUE`, uses a cached smoother to speed up computations. Default is `TRUE`.
#' @return A list containing:
#' \describe{
#'   \item{test_res}{A list of results from testing the fitted product distributions.}
#'   \item{gene_expr}{A list of gene expression matrices for each gene pair.}
#'   \item{cov_mat}{The covariate matrix used in the fitting.}
#'   \item{model_list}{The fitted model objects, if `return_models = TRUE`.}
#' }
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
#'   return_models = FALSE,
#'   ncores = 2,
#'   control = list(),
#'   seed = 123,
#'   local_testing = FALSE,
#'   preconstruct_smoother = TRUE
#' )
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
                   global_test = "LRT",
                   ncores = 2,
                   control = list(),
                   epsilon = 1e-6,
                   seed = 123,
                   preconstruct_smoother = TRUE,
                   return_models = FALSE,
                   return_coefs = FALSE,
                   check_morani = FALSE) {
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
