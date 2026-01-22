#' Fit product distributions to standardized gene expressions for a list of gene pairs.
#'
#' This function fits product distributions to gene pairs using their standardized expression values.
#' It applies a Generalized Additive Model (GAM) to each pair of genes from `gene_pair_list` to model
#' the joint effect while incorporating covariates.
#'
#' @param gene_pair_list_subset A data frame or matrix of gene pairs that passed filtering (e.g., via Moran's I test).
#' @param product_list A named list of product vectors corresponding to gene expression cross-products.
#' @param cov_mat A matrix or data frame containing covariates (e.g., spatial coordinates like `x1`, `x2`).
#' @param formula2 A formula or string specifying the GAM structure for modeling spatial interactions (e.g., `s(x1, x2)`).
#' @param family2 A distribution family to be used for model fitting. Default is `quasiproductr()`, which must be defined elsewhere.
#' @param control A list of control parameters passed to `mgcv::gam()`.
#' @param ncores Number of cores to use for parallel processing with `mclapply`.
#' @param global_test Method for global testing in product models. Options: `"lrt"` (likelihood ratio test) or `"wald"` (Wald-style smooth term test). Default is `"wald"`.
#' @param return_models Logical; if `TRUE`, returns the full model object for each gene pair.
#' @param return_coefs Logical; if `TRUE`, returns model coefficients and variance-covariance matrices.
#' @param return_pi Logical; if `TRUE`, returns predicted interval for fitted correlation.
#' @param preconstruct_smoother Logical; if `TRUE`, modifies the smoother basis (e.g., `'tp'` to `'tpcached'`) for caching and speed optimization. Default is `FALSE`.
#'
#' @return A list where each element corresponds to a gene pair. The contents depend on `return_models` and `return_coefs`:
#' \describe{
#'   \item{res_global}{Result of the global test (either p-value or list with EDF and p-value).}
#'   \item{fitted_rho}{Fitted values from the GAM.}
#'   \item{model}{(Optional) Full `mgcv::gam` object, if `return_models = TRUE`.}
#'   \item{model_coef}{(Optional) List with `beta` (coefficients) and `beta_cov` (variance-covariance matrix), if `return_coefs = TRUE`.}
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
#' # Compute gene pair products and optionally subset based on Moran's I
#' check_result <- check_products(
#'   gene_pair_list = test_data$gene_pair_list,
#'   marginals = marginal_res$marginal,
#'   cov_mat = test_data$cov_mat,
#'   check_morani = FALSE,
#'   ncores = 2
#' )
#'
#' # Fit product distributions to gene_pair_list_subset
#' model_list <- fit_products(
#'   gene_pair_list_subset = check_result$gene_pair_list_subset,
#'   product_list = check_result$product_list,
#'   cov_mat = test_data$cov_mat,
#'   formula2 = "s(x1, x2, bs='tp', k=50)",
#'   family2 = quasiproductr(),
#'   control = list(),
#'   ncores = 2,
#'   global_test = "wald",
#'   return_models = FALSE,
#'   return_coefs = FALSE,
#'   preconstruct_smoother = TRUE
#' )
#'
#' @importFrom parallel mclapply
#' @export

fit_products <- function(gene_pair_list_subset,
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
                         preconstruct_smoother = FALSE) {
  # Subset product list based on morani result
  product_list <- product_list[row.names(gene_pair_list_subset)]

  # Fit models for each cross product in parallel
  model_list <- parallel::mclapply(product_list, function(product) {
    fit_product(
      product = product,
      cov_mat = cov_mat,
      formula2 = formula2,
      family2 = family2,
      control = control,
      global_test = global_test,
      return_models = return_models,
      return_coefs = return_coefs,
      critical_value = critical_value,
      preconstruct_smoother = preconstruct_smoother
    )
  }, mc.cores = ncores)

  return(model_list)
}


fit_product <- function(product,
                        cov_mat,
                        formula2,
                        family2 = quasiproductr(),
                        control = list(),
                        global_test,
                        return_models = FALSE,
                        return_coefs = FALSE,
                        critical_value = 0.05,
                        return_pi = FALSE,
                        preconstruct_smoother) {
  # Extract product of expressions
  dat <- cbind(z = product, cov_mat)


  # Modify formula based on preconstruct_smoother
  if (preconstruct_smoother) {
    bs_type <- sub(".*bs='(.*?)'.*", "\\1", formula2)
    if (bs_type == "tp") {
      formula2 <- gsub("bs='tp'", "bs='tpcached'", formula2)
    } else if (bs_type == "gp") {
      formula2 <- gsub("bs='gp'", "bs='gpcached'", formula2)
    } else {
      stop("Invalid bs type in formula. Must be 'tp' or 'gp'.")
    }
  }

  # Create the formula
  mgcv_formula <- as.formula(paste("z ~", formula2))

  # Fit the GAM model
  model <- suppressWarnings({
    mgcv::gam(formula = mgcv_formula, family = family2, data = dat, control = control)
  })


  ## Global testing
  if (global_test == "lrt") {
    ## Likelihood ratio test

    # Fit the null model
    mgcv_formula0 <- as.formula(paste("z ~ 1"))
    model0 <- suppressWarnings({
      mgcv::gam(formula = mgcv_formula0, family = family2, data = dat, control = control)
    })

    res_lrt <- anova(model, model0, test = "LRT")
    global_p <- res_lrt[2, "Pr(>Chi)"]
    edf <- sum(model$edf)
    res_global <- list(global_p = global_p, edf = edf)
  } else if (global_test == "wald") {
    ## Wald test

    # Extract global p-value and edf from summary table
    global_p <- summary(model)$s.table[, "p-value"]
    edf <- sum(model$edf)
    res_global <- list(global_p = global_p, edf = edf)
  } else {
    stop("Invalid global testing method!")
  }

  # Extracting local estimation results
  fitted_rho <- model$fitted.values

  # Extracting predicted interval of fitted values
  if (return_pi) {
    p <- predict(model, type = "link", se.fit = TRUE) # include standard errors
    z_val <- qnorm(1 - critical_value / 2)
    fit_link <- p$fit
    se_link <- p$se.fit
    upper_link <- fit_link + z_val * se_link
    lower_link <- fit_link - z_val * se_link
    # Transform to response scale using the inverse link
    inv_link <- model$family$linkinv
    fit_response <- inv_link(fit_link)
    upper_response <- inv_link(upper_link)
    lower_response <- inv_link(lower_link)
    pred_interval <- data.frame(
      lower = lower_response,
      upper = upper_response
    )
  } else {
    pred_interval <- NULL
  }


  # Default is to return global testing result and local correlation estimate
  if (return_models) {
    # Include model
    return(list(
      res_global = res_global,
      fitted_rho = fitted_rho,
      pred_interval = pred_interval,
      model = model
    ))
  } else if (return_coefs) {
    # Only return model coef
    # Extracting estimation result from model
    beta <- coef(model)
    beta_cov <- vcov(model)
    model_coef <- list(beta = beta, beta_cov = beta_cov)
    return(list(
      res_global = res_global,
      model_coef = model_coef
    ))
  } else {
    # Fitted values
    return(list(
      res_global = res_global,
      fitted_rho = fitted_rho,
      pred_interval = pred_interval
    ))
  }
}
