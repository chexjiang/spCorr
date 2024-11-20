## Testing for significant spots
test_spots <- function(model, local_testing = FALSE){
  
  
  if(local_testing==TRUE){
    fitted_rho <- model$fitted.values
    
    # Extract the mean and covariance matrix of beta 
    beta_mean <- coef(model)
    beta_cov <- vcov(model)
    
    # Access the model frame and basis function information
    X <- predict(model, type = "lpmatrix")
    
    #### Transform the distribution of beta to rho by delta method
    # Mean
    rho_mean <- 0
    # Variance
    rho_cov <- (1-fitted_rho^2)^2 * X %*% beta_cov %*% t(X)
    rho_se <- sqrt(diag(rho_cov))
    
    # The p-values for each rho
    rho_z <- (fitted_rho-rho_mean) / rho_se
    rho_p <- 2 * (1 - pnorm(abs(rho_z)))
    
    return(list(fitted_rho=fitted_rho, rho_p=rho_p, rho_z=rho_z))
    
    
  }else{
    fitted_rho <- model$fitted.values
    
    return(list(fitted_rho=fitted_rho))
  }
}



## Testing function for a model
test_model <- function(model, local_testing = FALSE){
  
  # Global p-value and edf
  global_p <- summary(model)$s.table[, "p-value"]
  edf <- sum(model$edf)
  global_test <- list(global_p=global_p, edf=edf)
  
  # Fitted rho and local p-values
  test_res <- test_spots(model, local_testing = local_testing)
  fitted_rho <- test_res$fitted_rho
  rho_z <- test_res$rho_z
  rho_p <- test_res$rho_p
  cov <- model$model[,c('x1','x2')]
  
  if(local_testing==TRUE){
    local_test <- data.frame(fitted_rho=fitted_rho, rho_p=rho_p, rho_z=rho_z)
  }else{
    local_test <- data.frame(fitted_rho=fitted_rho)
  }
  
  list(global_test=global_test, local_test=local_test)
}


## Testing function for a model list
test_models <- function(model_list, ncores = ncores, local_testing = local_testing){
  
  test_list <- mclapply(model_list, function(model) {
    tryCatch({
      test_model(model, local_testing = local_testing)
    }, error = function(e) {
      message("Error with model: ", e$message)
      return(NULL)  # Return NULL if an error occurs
    })
  }, mc.cores = ncores)
  
}

