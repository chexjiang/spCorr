#' Define a custom family called 'quasiproductr' to be used with the mgcv package.
#'
#' The `quasiproductr` family models relationships where the link function is the
#' hyperbolic arctangent (`artanh`). This family defines custom link, variance,
#' residual, and initialization functions for use with Generalized Additive Models (GAM).
#'
#' @param link The link function to be used. Default is `"artanh"`.
#' @return A family object compatible with mgcv, including variance, link functions, and residual deviance calculation.
#' @export
quasiproductr <- function(link = "artanh") {
  if (link == "artanh") {
    linkfun <- function(mu) {
      mu <- pmin(pmax(mu, -0.9999), 0.9999)
      if (any(mu < -0.9999 | mu > 0.9999)) {
        warning(
          "Values of 'mu' must be within the range [-1, 1]. Adjusting the values."
        )
      }
      atanh(mu)
    }
    linkinv <- function(eta) pmin(pmax(tanh(eta), -0.9999), 0.9999)
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
  dev.resids <- function(y, mu, wt) {
    -2 * wt * (y * (atan(mu) - atan(y)) - 0.5 * log((mu^2 + 1) / (y^2 + 1)))
  }

  aic <- function(y, n, mu, wt, dev) NA

  initialize <- expression({
    if (any(is.infinite(y))) {
      stop("'quasiproductr' family does not support infinite response values.")
    }
    mustart <- tanh(y)
  })

  structure(
    list(
      family = "quasiproductr",
      link = name,
      linkfun = linkfun,
      linkinv = linkinv,
      variance = variance,
      dev.resids = dev.resids,
      aic = aic,
      mu.eta = mu.eta,
      initialize = initialize,
      validmu = validmu,
      valideta = valideta
    ),
    class = "family"
  )
}
