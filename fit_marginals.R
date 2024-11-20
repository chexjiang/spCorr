## Fit conditional margins and convert to standard gaussian for a gene list count_genes
covert_cond_marginals <- function(gene_list,
                                  count_mat,
                                  cov_mat,
                                  formula1,
                                  family1,
                                  to = 'gaussian',
                                  DT = TRUE,
                                  epsilon = 1e-6,
                                  ncores = ncores) {
  
  # Apply covert_cond_marginal to each gene in gene_list
  result <- mclapply(gene_list, function(gene) {
    tryCatch({
      covert_cond_marginal(gene, count_mat, cov_mat, formula1, family1, to, DT, epsilon)
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


covert_cond_marginal <- function(gene,
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



### Old version
covert_cond_marginal1 <- function(count_gene,
                                  cov_mat,
                                  formula1,
                                  family1,
                                  to=c('uniform', 'gaussian'),
                                  DT=TRUE,
                                  epsilon = 1e-6){
  y <- as.matrix(count_gene)
  dat <- cbind(y, cov_mat)
  mgcv_formula <-
    stats::formula(paste0("y~", formula1))
  
  
  # Convert marginals to Unif[0, 1].
  if(family1=='gaussian'){
    res <- mgcv::gam(mgcv_formula, data=dat, family = "gaussian")
    mean_vec <- stats::predict(res, type = "response")
    theta_vec <- rep(sqrt(res$sig2), length(mean_vec)) 
  }else if(family1=='poisson'){
    res <- mgcv::gam(mgcv_formula, data=dat, family = "poisson")
    mean_vec <- stats::predict(res, type = "response")
    theta_vec <- rep(NA, length(mean_vec))
  }else if(family1=='nb'){
    res <- mgcv::gam(mgcv_formula, data=dat, family = "nb")
    mean_vec <- stats::predict(res, type = "response")
    theta <- res$family$getTheta(TRUE)
    theta_vec <- rep(theta, length(mean_vec))
  }
  ## Frame
  if (!exists("zero_vec")) {
    #zero_vec <- rep(0, length(mean_vec))
    zero_vec <- 0
  }
  family_frame <- cbind(y, mean_vec, theta_vec, zero_vec)
  
  if(family1 == "gaussian"){
    pvec <- apply(family_frame, 1, function(x) {
      gamlss.dist::pNO(x[1], mu = x[2], sigma = abs(x[3]))
    })
    pvec2 <- pvec
    r <- pvec
  }else if (family1 == "poisson") {
    pvec <- apply(family_frame, 1, function(x) {
      stats::ppois(x[1], lambda = x[2])
    })
    
    if(DT==TRUE){
      pvec2 <- apply(family_frame, 1, function(x) {
        stats::ppois(x[1] - 1, lambda = x[2]) * as.integer(x[1] > 0)
      })
      
      u1 <- pvec
      u2 <- pvec2
      
      v <- stats::runif(length(pvec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
      
    }else{
      r <- pvec
    }
    
    
  } else if(family1 == "nb"){
    pvec <- apply(family_frame, 1, function(x) {
      stats::pnbinom(x[1], mu = x[2], size = x[3])
    })
    
    if(DT==TRUE){
      pvec2 <- apply(family_frame, 1, function(x) {
        stats::pnbinom(x[1] - 1, mu = x[2], size = x[3]) * as.integer(x[1] > 0)
      })
      
      u1 <- pvec
      u2 <- pvec2
      
      v <- stats::runif(length(pvec))
      ## Random mapping
      r <- u1 * v + u2 * (1 - v)
    }else{
      r <- pvec
    }
    
  }
  
  ## Avoid Inf
  idx_adjust <- which(1 - r < epsilon)
  r[idx_adjust] <- r[idx_adjust] - epsilon
  idx_adjust <- which(r < epsilon)
  r[idx_adjust] <- r[idx_adjust] + epsilon
  
  if(to=='uniform'){
    r
  }else if(to=='gaussian'){
    stats::qnorm(
      r,
      mean = 0,
      sd = 1,
      lower.tail = TRUE,
      log.p = FALSE
    )
  }
  
}
