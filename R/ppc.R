#' Predicted vs obseved plot function
#' 
#' The \code{ppc} function plots the predicted versus observed response values.
#' 
#' @param x An object of class \code{survFitTT} or \code{reproFitTT}
#' 
#' @export
#' @import ggplot2
#' 
ppc <- function(x) {
  # test class object
  if (! is(x, "reproFitTT") && ! is(x, "survFitTT"))
    stop("The object passed in argument [out] is not of class 'reproFitTT' or 'survFitTT' !\n")
  
}