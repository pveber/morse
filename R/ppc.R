#' Posterior predictive check plot
#' 
#' Plots posterior predictive check for \code{reproFitTT}, \code{survFitTT},
#'  \code{survFitTKTD}, \code{survFitCstExp} and \code{survFitVarExp} objects.
#' 
#' Depending on the class of the object \code{x} see their links. 
#' for class \code{reproFitTT}: \link[=ppc.reproFitTT]{ppc.reproFitTT} ;  
#' for class \code{survFitTT}: \link[=ppc.survFitTT]{ppc.survFitTT} ; 
#' for class \code{survFitTKTD}: \link[=ppc.survFitTKTD]{ppc.survFitTKTD} ;
#' for class \code{survFitCstExp}: \link[=ppc.survFitCstExp]{ppc.survFitCstExp} and
#' for class \code{survFitVarExp}: \link[=ppc.survFitVarExp]{ppc.survFitVarExp}.
#' 
#' @param x an object used to select a method \code{ppc}
#' @param \dots Further arguments to be passed to generic methods

#' @export
ppc <- function(x, ...){
  UseMethod("ppc")
}
