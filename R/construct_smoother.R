#### Pre-construction of the 2-dim cached smoother ####

#### tp ####
#' @importFrom mgcv Predict.matrix smooth.construct
#' @export
smooth.construct.tpcached.smooth.spec <- function(object, data, knots) {
  if (exists("tpcached_smoother", envir = .smoother_env)) {
    smooth_object <- get("tpcached_smoother", envir = .smoother_env)
    message("Extracting tpcached smoother")
  } else {
    object$bs <- "tp"
    smooth_object <- mgcv:::smooth.construct.tp.smooth.spec(object, data, knots)  # Use mgcv::: for internal function
    assign("tpcached_smoother", smooth_object, envir = .smoother_env)
    message("Constructing tpcached smoother")
  }
  return(smooth_object)
}

#' @export
# Prediction matrix function for the custom smoother
Predict.matrix.tpcached.smooth <- function(object, data) {
  # Use the Predict.matrix function of the underlying smoother
  mgcv:::Predict.matrix.tp.smooth(object, data)  # Use mgcv::: for internal function
}

#### gp ####
#' @export
smooth.construct.gpcached.smooth.spec <- function(object, data, knots) {
  if (exists("gpcached_smoother", envir = .smoother_env)) {
    smooth_object <- get("gpcached_smoother", envir = .smoother_env)
    message("Extracting gpcached smoother")
  } else {
    object$bs <- "gp"
    smooth_object <- mgcv:::smooth.construct.tp.smooth.spec(object, data, knots)  # Use mgcv::: for internal function
    assign("gpcached_smoother", smooth_object, envir = .smoother_env)
    message("Constructing gpcached smoother")
  }
  return(smooth_object)
}

#' @export
# Prediction matrix function for the custom smoother
Predict.matrix.gpcached.smooth <- function(object, data) {
  # Use the Predict.matrix function of the underlying smoother
  mgcv:::Predict.matrix.gp.smooth(object, data)  # Use mgcv::: for internal function
}
