#' Predict Lethal Concentration for x\% of the population at time \code{time_LCx} 
#' 
#' Predict median and 95\%CIs of the Lethal Concentration for x\% of the population.
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
LCx <- function(x, ...){
  UseMethod("LCx")
}
