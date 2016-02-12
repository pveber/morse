#' Plot Dose-response from raw data
#' 
#' Plots the survival rate as a function of the concentration (for a given
#' target time).
#' 
#' @param x an object used to select a method
#' @param \dots Further arguments to be passed to generic methods

#' @export
plotDoseResponse <- function(x, ...){
  UseMethod("plotDoseResponse")
}
