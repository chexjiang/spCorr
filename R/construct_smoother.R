#### Pre-construction of the 2-dim cached smoother


#### tp ####
smooth.construct.tpcached.smooth.spec <- function(object, data, knots) {
  if (exists("tpcached_smoother", envir = smoother_env)) {
    smooth_object <- get("tpcached_smoother", envir = smoother_env)
    #message("Extractting tpcached smoother")
  } else {
    object$bs <- "tp"  
    smooth_object <- smooth.construct.tp.smooth.spec(object, data, knots)
    assign("tpcached_smoother", smooth_object, envir = smoother_env)
    message("Constructing tpcached smoother")
  }
  return(smooth_object)
}


# Prediction matrix function for the custom smoother
Predict.matrix.tpcached.smooth <- function(object, data) {
  # Use the Predict.matrix function of the underlying smoother
  Predict.matrix.tp.smooth(object, data)
}



#### gp ####
smooth.construct.gpcached.smooth.spec <- function(object, data, knots) {
  if (exists("gpcached_smoother", envir = smoother_env)) {
    smooth_object <- get("gpcached_smoother", envir = smoother_env)
    #message("Extracting gpcached smoother")
  } else {
    object$bs <- "gp"  
    smooth_object <- smooth.construct.tp.smooth.spec(object, data, knots)
    assign("gpcached_smoother", smooth_object, envir = smoother_env)
    message("Constructing gpcached smoother")
  }
  return(smooth_object)
}

# Prediction matrix function for the custom smoother
Predict.matrix.gpcached.smooth <- function(object, data) {
  # Use the Predict.matrix function of the underlying smoother
  Predict.matrix.gp.smooth(object, data)
}
