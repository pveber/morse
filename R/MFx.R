#' Predict the Multiplication Factor leading to \eqn{x}\% of reduction in survival
#' at a specific time.
#' 
#' Predict  the x\% Multiplication Factor.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
MFx <- function(object, ...){
  UseMethod("MFx")
}

