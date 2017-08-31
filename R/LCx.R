#' Lethal Concentration for \code{x} percent of the population
#' 
#' 
#' @param x an object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
LCx <- function(x, ...){
  UseMethod("LCx")
}
