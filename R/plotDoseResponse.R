#' Plot dose-response from raw data
#' 
#' Plots the response of the effect as a function of the concentration at a given
#' target time.
#' 
#' @param x an object used to select a method \code{plotDoseRespons}
#' @param \dots Further arguments to be passed to generic methods

#' @export
plotDoseResponse <- function(x, ...){
  UseMethod("plotDoseResponse")
}
