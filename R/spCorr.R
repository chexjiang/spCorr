# Create a new environment for caching smoothers
smoother_env <- new.env()

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
                   family1 = 'nb',
                   formula2 = "s(x1, x2, bs='tp', k=50)",
                   family2 = quasiproductr(),
                   DT = TRUE,
                   return_models = FALSE,
                   ncores = 2,
                   control = list(),
                   seed = 123,
                   local_testing = FALSE,
                   preconstruct_smoother = TRUE){

  set.seed(seed)

  # Fit marginal to gene_list
  message("Start Marginal Fitting for ", length(gene_list), " genes")
  marginals <- fit_marginals(gene_list = gene_list,
                             count_mat = count_mat,
                             cov_mat = cov_mat,
                             formula1 = formula1,
                             family = family1,
                             to = 'gaussian',
                             DT = DT,
                             ncores = ncores)


  # Fit product to gene_pair_list
  message("Start Product Fitting for ", nrow(gene_pair_list), " gene pairs")
  #smoother_env <- new.env()
  model_list <- fit_products(gene_pair_list = gene_pair_list,
                             marginals = marginals,
                             cov_mat = cov_mat,
                             formula2 = formula2,
                             family2 = quasiproductr(),
                             control = control,
                             ncores = ncores,
                             preconstruct_smoother = preconstruct_smoother)

  # Extract the gene expression list
  gene_pair_expr_list <- lapply(1:nrow(gene_pair_list), function(i) {
    # Extract the expressions for the gene pair
    y1 <- count_mat[gene_pair_list[i, 1], ]
    y2 <- count_mat[gene_pair_list[i, 2], ]

    # Combine with cov_mat
    gene_pair_expr <- cbind(y1=y1, y2=y2)
    gene_pair_expr
  })
  names(gene_pair_expr_list) <- row.names(gene_pair_list)


  # Testing for gene_pair_list
  message("Start Testing for ", nrow(gene_pair_list), " gene pairs")
  test_res <- test_models(model_list, ncores, local_testing)

  message("Finished")

  if(return_models){
    return(list(test_res=test_res,
                gene_expr=gene_pair_expr_list,
                cov_mat=cov_mat,
                model_list=model_list))
  }else{
    return(list(test_res=test_res,
                gene_expr=gene_pair_expr_list,
                cov_mat=cov_mat))
  }
}

