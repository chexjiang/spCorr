# Extract the test results
extract_test <- function(res, name){
  
  test_res <- res$test_res
  
  if(name=="global"){
    # Extract the global p-value and edf for each pair
    global_tests <- lapply(test_res, function(x) x$global_test)
    res_global_test <- do.call(rbind, global_tests)
    res_global_test <- as.data.frame(res_global_test)
    global_pval <- unlist(res_global_test$global_p)
    edf <- unlist(res_global_test$edf)
    # FDR
    global_fdr <- p.adjust(global_pval, "fdr")
    res_global_test <- data.frame(global_fdr=global_fdr, global_pval=global_pval, edf=edf)
    
    return(res_global_test)
    
  }else if(name=="local"){
    # Extract the local p-values for each pair
    local_pval <- lapply(test_res, function(x) x$local_test$rho_p)
    local_pval <- do.call(rbind, local_pval)
    local_pval <- as.data.frame(local_pval)
    local_fdr <- apply(local_pval, 2, FUN=function(x) p.adjust(x, "fdr"))
    
    # Extract the local z-scores for each pair
    local_zval <- lapply(test_res, function(x) x$local_test$rho_p)
    local_zval <- do.call(rbind, local_zval)
    local_zval <- as.data.frame(local_zval)
    
    return(fdr=local_fdr, zval=local_zval)
  }
}


# Extract the fitted rho of each spot
extract_rho <- function(res){
  
  test_res <- res$test_res
  
  # Extract the fitted rho for each pair
  fitted_rho <- lapply(test_res, function(x) x$local_test$fitted_rho)
  fitted_rho <- do.call(rbind, fitted_rho)
  return(fitted_rho)
  
}




