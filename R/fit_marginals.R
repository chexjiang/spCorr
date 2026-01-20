#' Fit conditional margins to standard Gaussian distributions for a list of genes.
#'
#' This function applies marginal fitting and transforms each gene's expression to follow a standard Gaussian distribution.
#'
#' @param gene_list A vector of gene names or indices (row names or row numbers of `count_mat`) to process.
#' @param count_mat A matrix of raw gene expression counts. Rows correspond to genes, columns to observations (cells/spots).
#' @param cov_mat A matrix or data frame of covariates used for marginal modeling (e.g., spatial coordinates or experimental annotations).
#' @param formula1 A formula object or string (e.g., `"~ covariate1 + covariate2"`) specifying the model structure for the mean.
#' @param family1 A string specifying the distribution family to be used for modeling. Supported values include `"gaussian"`, `"poisson"`, `"nb"`, or `"zinb"`.
#' @param DT Logical; if `TRUE`, applies a discrete transformation suitable for count data. Default is `TRUE`.
#' @param epsilon A small numeric constant to avoid boundary issues (e.g., `0` or `1` values in uniform distribution). Default is `1e-6`.
#' @param ncores Integer specifying the number of cores to use for parallel processing via `parallel::mclapply`.
#' @param seed Random seed for reproducibility. Default is `123`.
#'
#' @return A list containing two matrices:
#' \describe{
#'   \item{marginal}{A matrix of transformed values from each gene, transformed to follow a standard normal distribution.}
#'   \item{residual}{A matrix of values transformed to standard uniform distribution (before applying the Gaussian quantile function).}
#' }
#'
#' @examples
#' data(test_data)
#' # Fit standardized marginals for gene expressions
#' marginal_res <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
#'   family1 = "nb",
#'   DT = TRUE,
#'   epsilon = 1e-6,
#'   ncores = 2
#' )
#'
#' @importFrom parallel mclapply
#' @export

fit_marginals <- function(gene_list,
                          count_mat,
                          cov_mat,
                          formula1,
                          family1,
                          DT = TRUE,
                          epsilon = 1e-6,
                          ncores = ncores,
                          seed = 123) {
  # Apply fit_marginal to each gene in gene_list
  result <- parallel::mclapply(gene_list, function(gene) {
    tryCatch(
      {
        fit_marginal(
          gene = gene,
          count_mat = count_mat,
          cov_mat = cov_mat,
          formula1 = formula1,
          family1 = family1,
          DT = DT,
          epsilon = epsilon,
          seed = seed
        )
      },
      error = function(e) {
        message("Error with model: ", e$message)
        return(NULL) # Return NULL if an error occurs
      }
    )
  }, mc.cores = ncores)


  ## Convert the result list to a matrix
  # residuals matrix
  residual_matrix <- do.call(rbind, lapply(result, function(x) x$uniform))
  row.names(residual_matrix) <- gene_list

  # marginals matrix
  marginal_matrix <- do.call(rbind, lapply(result, function(x) x$gaussian))
  row.names(marginal_matrix) <- gene_list

  return(list(
    marginal = marginal_matrix,
    residual = residual_matrix
  ))
}


fit_marginal <- function(gene,
                         count_mat,
                         cov_mat,
                         formula1,
                         family1,
                         DT = TRUE,
                         epsilon = 1e-6,
                         seed) {
  # Extract data for the specific gene
  y <- as.matrix(count_mat[gene, ])
  dat <- cbind(y, cov_mat)

  # Build the formula
  mgcv_formula <- stats::reformulate(formula1, response = "y")
  sigma_formula <- "~1"

  # Fit the model based on the family
  fit_model <- function(family1) {
    if (family1 == "poisson") {
      fit <- mgcv::gam(mgcv_formula, data = dat, family = family1)
    } else if (family1 == "nb") {
      if (sigma_formula != "~1") {
        dat$y <- round(dat$y)
        fit <- gamlss::gamlss(
          formula = mgcv_formula,
          sigma.formula = sigma_formula,
          data = dat,
          family = gamlss.dist::NBI,
          control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.01)
        )
      } else {
        fit <- mgcv::gam(mgcv_formula, data = dat, family = family1)
      }
    } else if (family1 == "zinb") {
      dat$y <- round(dat$y)
      sigma_formula <- stats::reformulate("1", response = "y")
      fit <- suppressWarnings({
        gamlss::gamlss(
          formula = mgcv_formula,
          sigma.formula = sigma_formula,
          nu.formula = mgcv_formula, ## Here nu is the dropout probability!
          data = dat,
          family = gamlss.dist::ZINBI,
          control = gamlss::gamlss.control(trace = FALSE, c.crit = 0.01)
        )
      })

      fit
    }
  }


  # Get parameters
  get_params <- switch(family1,
    gaussian = {
      res <- fit_model("gaussian")
      list(
        mean_vec = stats::predict(res, type = "response"),
        theta_vec = rep(sqrt(res$sig2), nrow(y)),
        zero_vec = rep(0, nrow(y))
      )
    },
    poisson = {
      res <- fit_model("poisson")
      list(
        mean_vec = stats::predict(res, type = "response"),
        theta_vec = NA,
        zero_vec = rep(0, nrow(y))
      )
    },
    nb = {
      res <- fit_model("nb")
      if (class(res)[1] == "gam") {
        list(
          mean_vec = stats::predict(res, type = "response"),
          theta_vec = rep(res$family$getTheta(TRUE), nrow(y)),
          zero_vec = rep(0, nrow(y))
        )
      } else if (class(res)[1] == "gamlss") {
        list(
          mean_vec = stats::predict(res, type = "response"),
          theta_vec = 1 / stats::predict(res, type = "response", what = "sigma", data = dat),
          zero_vec = rep(0, nrow(y))
        )
      }
    },
    zinb = {
      res <- fit_model("zinb")
      list(
        mean_vec = stats::predict(res, type = "response"),
        theta_vec = stats::predict(res, type = "response", what = "sigma", data = dat),
        zero_vec = stats::predict(res, type = "response", what = "nu", data = dat)
      )
    }
  )
  # Frame
  family_frame <- cbind(y, get_params$mean_vec, get_params$theta_vec, get_params$zero_vec)

  calc_pvec <- function(x) {
    set.seed(seed)
    switch(family1,
      gaussian = gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3])),
      poisson = stats::ppois(x[1], lambda = x[2]),
      nb = stats::pnbinom(x[1], mu = x[2], size = x[3]),
      zinb = gamlss.dist::pZINBI(x[1], mu = x[2], sigma = abs(x[3]), nu = x[4])
    )
  }
  pvec <- apply(family_frame, 1, calc_pvec)


  # Apply discrete transformation if necessary
  if (DT && family1 %in% c("poisson", "nb", "zinb")) {
    calc_pvec2 <- function(x) {
      set.seed(seed)
      switch(family1,
        poisson = stats::ppois(x[1] - 1, lambda = x[2]),
        nb = stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]),
        zinb = ifelse(x[1] > 0,
          gamlss.dist::pZINBI(x[1] - 1, mu = x[2], sigma = abs(x[3]), nu = x[4]),
          0
        )
      ) * as.integer(x[1] > 0)
    }

    pvec2 <- apply(family_frame, 1, calc_pvec2)

    set.seed(seed)
    v <- stats::runif(length(pvec))
    r <- pvec * v + pvec2 * (1 - v)
  } else {
    r <- pvec
  }

  # Avoid boundary values
  r[r < epsilon] <- epsilon
  r[1 - r < epsilon] <- 1 - epsilon


  # Return transformed values
  list(uniform = r, gaussian = stats::qnorm(r))
}
