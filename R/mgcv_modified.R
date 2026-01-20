#### Re-define the fix.family.link.family and fix.family.var in mgcv ####

#' Redefine fix.family.link.family to add derivatives to a family object
#'
#' This function adds the second, third, and fourth derivatives of the link function
#' with respect to mu to the family object, which is used for Newton-like optimization.
#'
#' @param fam A family object for which to add derivatives.
#' @return A modified family object with additional components for derivatives.
#' @importFrom stats qnorm dnorm pnorm
#' @importFrom mgcv fix.family.link
#' @export
fix.family.link.family <- function(fam)
                                   # adds d2link the second derivative of the link function w.r.t. mu
                                   # to the family supplied, as well as a 3rd derivative function
                                   # d3link...
# All d2link and d3link functions have been checked numerically.
{
  if (!inherits(fam, "family")) stop("fam not a family object")
  if (is.null(fam$canonical)) { ## note the canonical link - saves effort in full Newton
    if (fam$family == "gaussian") {
      fam$canonical <- "identity"
    } else if (fam$family == "poisson" || fam$family == "quasipoisson") {
      fam$canonical <- "log"
    } else if (fam$family == "binomial" || fam$family == "quasibinomial") {
      fam$canonical <- "logit"
    } else if (fam$family == "Gamma") {
      fam$canonical <- "inverse"
    } else if (fam$family == "inverse.gaussian") {
      fam$canonical <- "1/mu^2"
    } else {
      fam$canonical <- "none"
    }
  }
  if (!is.null(fam$d2link) && !is.null(fam$d3link) && !is.null(fam$d4link)) {
    return(fam)
  }
  link <- fam$link
  if (length(link) > 1) {
    if (fam$family == "quasi") # then it's a power link
      {
        lambda <- log(fam$linkfun(exp(1))) ## the power, if > 0
        if (lambda <= 0) {
          fam$d2link <- function(mu) -1 / mu^2
          fam$d3link <- function(mu) 2 / mu^3
          fam$d4link <- function(mu) -6 / mu^4
        } else {
          fam$d2link <- function(mu) lambda * (lambda - 1) * mu^(lambda - 2)
          fam$d3link <- function(mu) (lambda - 2) * (lambda - 1) * lambda * mu^(lambda - 3)
          fam$d4link <- function(mu) (lambda - 3) * (lambda - 2) * (lambda - 1) * lambda * mu^(lambda - 4)
        }
      } else {
      stop("unrecognized (vector?) link")
    }
  } else if (link == "identity") {
    fam$d4link <- fam$d3link <- fam$d2link <-
      function(mu) rep.int(0, length(mu))
  } else if (link == "log") {
    fam$d2link <- function(mu) -1 / mu^2
    fam$d3link <- function(mu) 2 / mu^3
    fam$d4link <- function(mu) -6 / mu^4
  } else if (link == "inverse") {
    fam$d2link <- function(mu) 2 / mu^3
    fam$d3link <- function(mu) {
      mu <- mu * mu
      -6 / (mu * mu)
    }
    fam$d4link <- function(mu) {
      mu2 <- mu * mu
      24 / (mu2 * mu2 * mu)
    }
  } else if (link == "logit") {
    fam$d2link <- function(mu) 1 / (1 - mu)^2 - 1 / mu^2
    fam$d3link <- function(mu) 2 / (1 - mu)^3 + 2 / mu^3
    fam$d4link <- function(mu) 6 / (1 - mu)^4 - 6 / mu^4
  } else if (link == "probit") {
    fam$d2link <- function(mu) {
      # eta <- fam$linkfun(mu)
      eta <- qnorm(mu)
      # eta/fam$mu.eta(eta)^2
      eta / pmax(dnorm(eta), .Machine$double.eps)^2
    }
    fam$d3link <- function(mu) {
      # eta <-  fam$linkfun(mu)
      eta <- qnorm(mu)
      # (1 + 2*eta^2)/fam$mu.eta(eta)^3
      (1 + 2 * eta^2) / pmax(dnorm(eta), .Machine$double.eps)^3
    }
    fam$d4link <- function(mu) {
      # eta <-  fam$linkfun(mu)
      eta <- qnorm(mu)
      # (7*eta + 6*eta^3)/fam$mu.eta(eta)^4
      (7 * eta + 6 * eta^3) / pmax(dnorm(eta), .Machine$double.eps)^4
    }
  } else if (link == "cloglog") {
    ## g = log(-log(1-mu)), g' = -1/(log(1-mu)*(1-mu))
    fam$d2link <- function(mu) {
      l1m <- log1p(-mu)
      -1 / ((1 - mu)^2 * l1m) * (1 + 1 / l1m)
    }
    fam$d3link <- function(mu) {
      l1m <- log1p(-mu)
      mu3 <- (1 - mu)^3
      (-2 - 3 * l1m - 2 * l1m^2) / mu3 / l1m^3
    }
    fam$d4link <- function(mu) {
      l1m <- log1p(-mu)
      mu4 <- (1 - mu)^4
      (-12 - 11 * l1m - 6 * l1m^2 - 6 / l1m) / mu4 / l1m^3
    }
  } else if (link == "sqrt") {
    fam$d2link <- function(mu) -.25 * mu^-1.5
    fam$d3link <- function(mu) .375 * mu^-2.5
    fam$d4link <- function(mu) -0.9375 * mu^-3.5
  } else if (link == "cauchit") {
    ## uses general result that if link is a quantile function then
    ## d mu / d eta = f(eta) where f is the density. Link derivative
    ## is one over this... repeated differentiation w.r.t. mu using chain
    ## rule gives results...
    fam$d2link <- function(mu) {
      # eta <- fam$linkfun(mu)
      eta <- qcauchy(mu)
      2 * pi * pi * eta * (1 + eta * eta)
    }
    fam$d3link <- function(mu) {
      # eta <- fam$linkfun(mu)
      eta <- qcauchy(mu)
      eta2 <- eta * eta
      2 * pi * pi * pi * (1 + 3 * eta2) * (1 + eta2)
    }
    fam$d4link <- function(mu) {
      # eta <- fam$linkfun(mu)
      eta <- qcauchy(mu)
      eta2 <- eta * eta
      2 * pi^4 * (8 * eta + 12 * eta2 * eta) * (1 + eta2)
    }
  } else if (link == "1/mu^2") {
    fam$d2link <- function(mu) 6 * mu^-4
    fam$d3link <- function(mu) -24 * mu^-5
    fam$d4link <- function(mu) 120 * mu^-6
  } else if (substr(link, 1, 3) == "mu^") { ## it's a power link
    ## note that lambda <=0 gives log link so don't end up here
    lambda <- get("lambda", environment(fam$linkfun))
    fam$d2link <- function(mu) {
      (lambda * (lambda - 1)) * mu^{
        lambda - 2
      }
    }
    fam$d3link <- function(mu) {
      (lambda * (lambda - 1) * (lambda - 2)) * mu^{
        lambda - 3
      }
    }
    fam$d4link <- function(mu) {
      (lambda * (lambda - 1) * (lambda - 2) * (lambda - 3)) * mu^{
        lambda - 4
      }
    }
  } else if (link == "artanh") { ## it's a power link
    ## note that lambda <=0 gives log link so don't end up here
    fam$d2link <- function(mu) {
      mu_plus <- 1 / (1 + mu)
      mu_minus <- 1 / (1 - mu)
      (-mu_plus^2 + mu_minus^2) / (mu_plus + mu_minus)
    } ## g2g productr
    fam$d3link <- function(mu) {
      mu_plus <- 1 / (1 + mu)
      mu_minus <- 1 / (1 - mu)
      2 * (mu_plus^3 + mu_minus^3) / (mu_plus + mu_minus)
    } ## g3g productr
    fam$d4link <- function(mu) {
      mu_plus <- 1 / (1 + mu)
      mu_minus <- 1 / (1 - mu)
      6 * (-mu_plus^4 + mu_minus^4) / (mu_plus + mu_minus)
    } ## g4g productr
  } else {
    stop("link not recognised")
  }
  ## avoid giant environments being stored....
  environment(fam$d2link) <- environment(fam$d3link) <- environment(fam$d4link) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.link.family


#' @export
fix.family.var <- function(fam)
                           # adds dvar the derivative of the variance function w.r.t. mu
                           # to the family supplied, as well as d2var the 2nd derivative of
# the variance function w.r.t. the mean. (All checked numerically).
{
  if (inherits(fam, "extended.family")) {
    return(fam)
  }
  if (!inherits(fam, "family")) stop("fam not a family object")
  if (!is.null(fam$dvar) && !is.null(fam$d2var) && !is.null(fam$d3var)) {
    return(fam)
  }
  family <- fam$family
  fam$scale <- -1
  if (family == "gaussian") {
    fam$d3var <- fam$d2var <- fam$dvar <- function(mu) rep.int(0, length(mu))
  } else if (family == "poisson" || family == "quasipoisson") {
    fam$dvar <- function(mu) rep.int(1, length(mu))
    fam$d3var <- fam$d2var <- function(mu) rep.int(0, length(mu))
    if (family == "poisson") fam$scale <- 1
  } else if (family == "binomial" || family == "quasibinomial") {
    fam$dvar <- function(mu) 1 - 2 * mu
    fam$d2var <- function(mu) rep.int(-2, length(mu))
    fam$d3var <- function(mu) rep.int(0, length(mu))
    if (family == "binomial") fam$scale <- 1
  } else if (family == "Gamma") {
    fam$dvar <- function(mu) 2 * mu
    fam$d2var <- function(mu) rep.int(2, length(mu))
    fam$d3var <- function(mu) rep.int(0, length(mu))
  } else if (family == "quasi") {
    fam$dvar <- switch(fam$varfun,
      constant = function(mu) rep.int(0, length(mu)),
      "mu(1-mu)" = function(mu) 1 - 2 * mu,
      mu = function(mu) rep.int(1, length(mu)),
      "mu^2" = function(mu) 2 * mu,
      "mu^3" = function(mu) 3 * mu^2
    )
    if (is.null(fam$dvar)) stop("variance function not recognized for quasi")
    fam$d2var <- switch(fam$varfun,
      constant = function(mu) rep.int(0, length(mu)),
      "mu(1-mu)" = function(mu) rep.int(-2, length(mu)),
      mu = function(mu) rep.int(0, length(mu)),
      "mu^2" = function(mu) rep.int(2, length(mu)),
      "mu^3" = function(mu) 6 * mu
    )
    fam$d3var <- switch(fam$varfun,
      constant = function(mu) rep.int(0, length(mu)),
      "mu(1-mu)" = function(mu) rep.int(0, length(mu)),
      mu = function(mu) rep.int(0, length(mu)),
      "mu^2" = function(mu) rep.int(0, length(mu)),
      "mu^3" = function(mu) rep.int(6, length(mu))
    )
  } else if (family == "inverse.gaussian") {
    fam$dvar <- function(mu) 3 * mu^2
    fam$d2var <- function(mu) 6 * mu
    fam$d3var <- function(mu) rep.int(6, length(mu))
  } else if (family == "quasiproductr") {
    fam$dvar <- function(mu) 2 * mu
    fam$d2var <- function(mu) rep.int(2, length(mu))
    fam$d3var <- function(mu) rep.int(0, length(mu))
  } else {
    stop("family not recognised")
  }
  environment(fam$dvar) <- environment(fam$d2var) <- environment(fam$d3var) <- environment(fam$linkfun)
  return(fam)
} ## fix.family.var
