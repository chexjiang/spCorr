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
#' @examples
#' data(test_data)
#' marginals <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
#'   family1 = "nb",
#'   to = "gaussian",
#'   DT = TRUE,
#'   epsilon = 1e-6,
#'   ncores = 2
#' )
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
  
  # Extract data for the specific gene
  y <- as.matrix(count_mat[gene,])
  dat <- cbind(y, cov_mat)
  
  # Build the formula
  mgcv_formula <- stats::reformulate(formula1, response = "y")
  
  
  # Fit the model based on the family
  fit_model <- function(fam) {
    mgcv::gam(mgcv_formula, data = dat, family = fam)
  }
  
  
  # Get parameters
  params <- switch(
    family1,
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

  # Generate p-values
  pvec <- mapply(function(y, mean, theta) {
    switch(
      family1,
      gaussian = gamlss.dist::pNO(y, mu = mean, sigma = abs(theta)),
      poisson = stats::ppois(y, lambda = mean),
      nb = stats::pnbinom(y, mu = mean, size = theta)
    )
  }, y, params$mean_vec, params$theta_vec)
  
  
  
  # Apply discrete transformation if necessary
  if (DT && family1 %in% c("poisson", "nb")) {
    pvec2 <- mapply(function(y, mean, theta) {
      if (family1 == "poisson") {
        stats::ppois(y - 1, lambda = mean) * as.integer(y > 0)
      } else {
        stats::pnbinom(y - 1, mu = mean, size = theta) * as.integer(y > 0)
      }
    }, y, params$mean_vec, params$theta_vec)
    v <- stats::runif(length(pvec))
    pvec <- pvec * v + pvec2 * (1 - v)
  }

  # Avoid boundary values
  pvec[pvec < epsilon] <- epsilon
  pvec[1 - pvec < epsilon] <- 1 - epsilon
  
  
  
  # Return transformed values
  if (to == "uniform") {
    return(pvec)
  } else {
    return(stats::qnorm(pvec))
  }
}



