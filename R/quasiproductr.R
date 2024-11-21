#' Define a custom family called 'quasiproductr' to be used with the mgcv package.
#'
#' The `quasiproductr` family is a custom family added to the `mgcv` package to model
#' relationships where the link function is the hyperbolic arctangent (`artanh`).
#' This family includes functions to define the link, variance, residuals, and
#' initialization logic required for fitting with Generalized Additive Models (GAM).
#'
#' @param link The link function to be used. Default is `"artanh"`.
#' @return A family object compatible with mgcv, including variance, link functions, and residual deviance calculation.
#' @examples
#' # Example of using quasiproductr with mgcv's gam function
#' library(mgcv)
#' fit <- gam(y ~ s(x), family = quasiproductr(), data = my_data)
#' @export
quasiproductr <- function (link = "artanh"){
  if (link == "artanh") {
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
    name <- "artanh"

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




