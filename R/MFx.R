#' Predict the Multiplication Factor leading to x\% of reduction in survival
#' at a specific time.
#' 
#' Generic method for \code{MFx}, a function denoted \eqn{MF(x,t)} for 
#' \eqn{x}\% Multiplication Factor at time \eqn{t}.
#' 
#' When the \code{object} is of class \code{survFit}, see \link[=MFx.survFit]{MFx.survFit}
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
MFx <- function(object, ...){
  UseMethod("MFx")
}

