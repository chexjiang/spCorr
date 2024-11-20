## The quasiproductr family added to mgcv family class
quasiproductr <- function (link = "arctanh"){
  if (link == "arctanh") {
    linkfun <- function(mu) {
      mu <- pmin(pmax(mu, -0.9999), 0.9999)  # Constrain input
      if(any(mu < -0.9999 | mu > 0.9999)) {
        warning("Some mu values are out of bounds")
        print(summary(mu))
      }
      atanh(mu)
    }
    linkinv <- function(eta) {
      mu <- tanh(eta)  # This naturally constrains output to (-1, 1)
      pmin(pmax(mu, -0.9999), 0.9999)  # Extra safety
    }
    mu.eta <- function(eta) 1 / (cosh(eta)^2)
    valideta <- function(eta) TRUE
    name <- "arctanh"
    
  } else {
    linktemp <- make.link(link)
    linkfun <- linktemp$linkfun
    linkinv <- linktemp$linkinv
    mu.eta <- linktemp$mu.eta
    valideta <- linktemp$valideta
  }
  
  variance <- function(mu) 1 + mu^2
  validmu <- function(mu) all(mu >= -1 & mu <= 1)
  dev.resids <- function(y, mu, wt)
  { 
    -2*wt * (y*(atan(mu) - atan(y)) - 1/2*log((mu^2+1)/(y^2+1))) 
  }
  
  aic <- function(y, n, mu, wt, dev) NA
  
  initialize <- expression({
    if (any(is.infinite(y)))
      stop("Any infinite values not allowed for the 'quasiProductr' family")
    #n <- rep.int(1, nobs)
    mustart <- tanh(y)
  })
  structure(list(family = "quasiproductr",
                 link = name,
                 linkfun = linkfun,
                 linkinv = linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 aic = aic,
                 mu.eta = mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = valideta),
            class = "family")
}





## Sampling function for product family
rproduct <- function(n, rho){
  sample_bivariate <- function(correlation, n, mu1, mu2) {
    mu <- c(mu1, mu2)  # Mean vector for gi and gj based on mu1 and mu2
    Sigma <- matrix(c(1, correlation, correlation, 1), 2, 2)  # Covariance matrix
    sample <- mvrnorm(n, mu = mu, Sigma = Sigma)  # Sampling
    return(sample)
  }
  samples <- sample_bivariate(rho, n, 0, 0)
  # Ensure samples is treated as a matrix
  if (is.vector(samples)) samples <- matrix(samples, ncol = 2)
  z <- samples[, 1] * samples[, 2]
  return(z)
}

sample_bivariate <- function(correlation, n, mu1, mu2) {
  mu <- c(mu1, mu2)  # Mean vector for gi and gj based on mu1 and mu2
  Sigma <- matrix(c(1, correlation, correlation, 1), 2, 2)  # Covariance matrix
  sample <- mvrnorm(n, mu = mu, Sigma = Sigma)  # Sampling
  return(sample)
}

