#' Print of \code{survFitVarExp} object
#'
#' This is the generic \code{print} S3 method for the \code{survFitVarExp} class.
#' It prints the underlying JAGS model and some information on the Bayesian
#' inference procedure.
#'
#' @param x An object of class \code{survFitVarExp}
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords print
#' 
#' @examples
#' # (1) Load the data
#' data(propiconazole_pulse_exposure)
#' 
#' # (2) Create a survData object
#' dataset <- survData(propiconazole_pulse_exposure)
#' 
#' \dontrun{
#' # (3) Run the survFit function with TK-TD model 'SD' or 'IT' 
#' out <- survFit(dataset, model_type="SD")
#' 
#' # (4) Print the survFit object
#' print(out)
#' }
#'
#' @export
print.survFitVarExp <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class survFitTKTD
  
  mcmcInfo = x$mcmcInfo
  
  # M.C.M.C. informations
  nbr.thin = mcmcInfo$nbr.thin
  mcmc_info =
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("Number of iterations per chain = ", mcmcInfo$n.iter, "\n")
  cat("Thinning interval =", mcmcInfo$thin.interval, "\n")
  cat("Number of chains =", mcmcInfo$n.chains, "\n")
  cat("Number iterations in warmup per chain =", mcmcInfo$n.warmup, "\n")
  cat("Sample size per chain =", mcmcInfo$n.iter / mcmcInfo$thin.interval , "\n")
}