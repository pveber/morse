#' Predict \eqn{X}\% Lethal Concentration at the maximum time point (default).
#' 
#' Predict median and 95\% credible interval of the x\% Lethal Concentration.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
LCx <- function(object, ...){
  UseMethod("LCx")
}
