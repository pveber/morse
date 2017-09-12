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
#' dat <- survData(propiconazole_pulse_exposure)
#' 
#' \dontrun{
#' # (3) Run the survFit function
#' out <- survFit(dat, model_type="SD")
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
  nbr.thin = mcmc_info$nbr.thin
  mcmc_info =
    cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Number of iterations per chain = ", mcmcInfo$nbr.iter, "\n")
  cat("Thinning interval =", mcmcInfo$nbr.thin, "\n")
  cat("Number of chains =", mcmcInfo$nbr.chains, "\n")
  cat("Number iterations in warmup per chain =", mcmcInfo$nbr.warmups, "\n")
  cat("Total time computing (specific of computer) =", mcmcInfo$total_time, "\n")
  cat("Sample size per chain =", mcmcInfo$nbr.iter / mcmcInfo$nbr.thin , "\n")
}