#' Predict \eqn{X}\% Lethal Concentration at the maximum time point (default).
#' 
#' Predict median and 95\% credible interval of the x\% Lethal Concentration.
#' 
#' When class of \code{object} is \code{survFit}, see \link[=LCx.survFit]{LCx.survFit}.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
LCx <- function(object, ...){
  UseMethod("LCx")
}
