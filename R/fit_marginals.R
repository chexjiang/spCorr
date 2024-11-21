#' Fit conditional margins to standard Gaussian distributions for a list of genes.
#'
#' This function applies marginal fitting and transforms each gene's expression to follow a standard Gaussian distribution.
#'
#' @param gene_list A vector of gene names or indices for which the conditional margins are to be fit.
#' @param count_mat A matrix of gene expression counts where rows represent genes and columns represent observations.
#' @param cov_mat A matrix or data frame of covariates used for fitting the marginal distributions.
#' @param formula1 A formula object or string specifying the model to be used for marginal fitting (e.g., `~ covariate`).
#' @param family1 The distribution family to be used for fitting. Supported options include `'gaussian'`, `'poisson'`, or `'nb'`.
#' @param to The target distribution to transform to. Options are `'uniform'` or `'gaussian'`. Default is `'gaussian'`.
#' @param DT A logical value indicating whether discrete transformation should be applied. Default is `TRUE`.
#' @param epsilon A small number to avoid boundary values during transformation. Default is `1e-6`.
#' @param ncores The number of CPU cores to use for parallel processing. Default is the number of available cores.
#' @return A matrix where each row represents the transformed values for a gene.
#' @param model_list A list of fitted GAM model objects to be tested.
#' @param ncores The number of cores to use for parallel processing.
#' @param local_testing Logical; whether to perform local significance testing for each model. Default is `FALSE`.
#' @return A list containing test results for each model in the input `model_list`.
#' Each element contains global test results (`global_p` and `edf`) and local test results.
#' @examples
#' data(test_data)
#' # Fit standardized marginals for gene expressions
#' marginals <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "~ covariate",
#'   family1 = "nb",
#'   to = "gaussian",
#'   DT = TRUE,
#'   ncores = 2
#' )
#'
#' # Fit product distributions for gene pairs
#' model_list <- fit_products(
#'   gene_pair_list = test_data$gene_pair_list,
#'   marginals = marginals,
#'   cov_mat = test_data$cov_mat,
#'   formula2 = "s(x1, x2, bs='tp', k=50)",
#'   family2 = quasiproductr(),
#'   control = list(),
#'   ncores = 2,
#'   preconstruct_smoother = TRUE
#' )
#'
#' @importFrom parallel mclapply
#' @export

fit_marginals <- function(gene_list,
                          count_mat,
                          cov_mat,
                          formula1,
                          family1,
                          to = 'gaussian',
                          DT = TRUE,
                          epsilon = 1e-6,
                          ncores = ncores) {

  # Apply covert_cond_marginal to each gene in gene_list
  result <- parallel::mclapply(gene_list, function(gene) {
    tryCatch({
      fit_marginal(gene, count_mat, cov_mat, formula1, family1, to, DT, epsilon)
    }, error = function(e) {
      message("Error with model: ", e$message)
      return(NULL)  # Return NULL if an error occurs
    })
  }, mc.cores = ncores)


  # Convert the result list to a matrix
  result_matrix <- do.call(rbind, result)
  row.names(result_matrix) <- gene_list

  return(result_matrix)
}


fit_marginal <- function(gene,
                         count_mat,
                         cov_mat,
                         formula1,
                         family1,
                         to = c('uniform', 'gaussian'),
                         DT = TRUE,
                         epsilon = 1e-6) {
  y <- as.matrix(count_mat[gene,])
  dat <- cbind(y, cov_mat)
  mgcv_formula <- stats::reformulate(formula1, response = "y")

  fit_model <- function(fam) {
    mgcv::gam(mgcv_formula, data = dat, family = fam)
  }

  get_params <- switch(family1,
                       gaussian = {
                         res <- fit_model("gaussian")
                         list(
                           mean_vec = stats::predict(res, type = "response"),
                           theta_vec = rep(sqrt(res$sig2), nrow(y))
                         )
                       },
                       poisson = {
                         res <- fit_model("poisson")
                         list(
                           mean_vec = stats::predict(res, type = "response"),
                           theta_vec = NA
                         )
                       },
                       nb = {
                         res <- fit_model("nb")
                         list(
                           mean_vec = stats::predict(res, type = "response"),
                           theta_vec = rep(res$family$getTheta(TRUE), nrow(y))
                         )
                       }
  )

  family_frame <- cbind(y, get_params$mean_vec, get_params$theta_vec, 0)

  calc_pvec <- function(x) {
    switch(family1,
           gaussian = gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3])),
           poisson = stats::ppois(x[1], lambda = x[2]),
           nb = stats::pnbinom(x[1], mu = x[2], size = x[3])
    )
  }

  pvec <- apply(family_frame, 1, calc_pvec)

  if (DT && family1 %in% c("poisson", "nb")) {
    calc_pvec2 <- function(x) {
      switch(family1,
             poisson = stats::ppois(x[1] - 1, lambda = x[2]),
             nb = stats::pnbinom(x[1] - 1, mu = x[2], size = x[3])
      ) * as.integer(x[1] > 0)
    }

    pvec2 <- apply(family_frame, 1, calc_pvec2)
    v <- stats::runif(length(pvec))
    r <- pvec * v + pvec2 * (1 - v)
  } else {
    r <- pvec
  }

  r[r < epsilon] <- epsilon
  r[1 - r < epsilon] <- 1 - epsilon

  if (to == 'uniform') {
    return(r)
  } else {
    return(stats::qnorm(r))
  }
}



