#' Testing function for a list of models
#'
#' This function performs significance testing on a list of fitted models using parallel processing.
#' It internally calls `test_model()` for each model in the provided list.
#'
#' @param model_list A list of fitted Generalized Additive Model (GAM) objects to be tested.
#' @param ncores An integer specifying the number of cores to use for parallel processing.
#' @param local_testing Logical; if `TRUE`, performs local significance testing for each model. Default is `FALSE`.
#' @return A list containing test results for each model in the input `model_list`.
#' Each element of the list includes:
#' \describe{
#'   \item{global_test}{A list containing global p-values (`global_p`) and effective degrees of freedom (`edf`).}
#'   \item{local_test}{A data frame containing local p-values and other local test metrics, if `local_testing = TRUE`.}
#' }
#' @examples
#' data(test_data)
#'
#' # Fit standardized marginals for gene expressions
#' marginals <- fit_marginals(
#'   gene_list = test_data$gene_list,
#'   count_mat = test_data$count_mat,
#'   cov_mat = test_data$cov_mat,
#'   formula1 = "layer_annotations",
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
#' # Perform significance testing on the fitted models
#' test_res <- test_models(
#'   model_list = model_list,
#'   ncores = 2,
#'   local_testing = FALSE
#' )
#'
#' @importFrom parallel mclapply
#' @export
test_models <- function(model_list, ncores = 2, local_testing = FALSE){

  # Perform testing
  test_list <- parallel::mclapply(model_list, function(model) {
    tryCatch(
      {
        test_model(model, local_testing = local_testing)
      },
      error = function(e) {
        message("Error in model testing: ", e$message)
        list(global_test = list(global_p = NA, edf = NA), local_test = NA)
      }
    )
  }, mc.cores = ncores)
  
  return(test_list)
}



## Perform global and local testing for a GAM model
test_model <- function(model, local_testing = FALSE){

  # Global testing
  global_p <- summary(model)$s.table[, "p-value"]
  edf <- sum(model$edf)
  global_test <- list(global_p=global_p, edf=edf)

  # Local testing
  local_test <- test_spots(model, local_testing = local_testing)
  fitted_rho <- local_test$fitted_rho
  rho_z <- local_test$rho_z
  rho_p <- local_test$rho_p

  if(local_testing==TRUE){
    local_test <- data.frame(fitted_rho=fitted_rho, rho_p=rho_p, rho_z=rho_z)
  }else{
    local_test <- data.frame(fitted_rho=fitted_rho)
  }

  list(global_test=global_test, local_test=local_test)
}



## Testing for significant spots
test_spots <- function(model, local_testing = FALSE){
  
  # Ensure model is not atomic
  if (is.atomic(model)) {
    stop("Model should not be an atomic vector")
  }
  
  fitted_rho <- model$fitted.values
  
  if(local_testing){
    
    # Extract the mean and covariance matrix of beta
    beta_mean <- coef(model)
    beta_cov <- vcov(model)
    
    # Access the model frame and basis function information
    X <- predict(model, type = "lpmatrix")
    
    # Delta method for mean and variance
    rho_cov <- (1 - fitted_rho^2)^2 * X %*% beta_cov %*% t(X)
    rho_se <- sqrt(diag(rho_cov))
    
    # Compute z-scores and p-values
    rho_z <- (fitted_rho) / rho_se
    rho_p <- 2 * (1 - pnorm(abs(rho_z)))
    
    return(list(fitted_rho=fitted_rho, rho_p=rho_p, rho_z=rho_z))
    
  }else{
    return(list(fitted_rho=fitted_rho))
  }
}


