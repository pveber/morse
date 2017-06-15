#' Print of \code{gm_survFitTKTD} object
#'
#' This is the generic \code{print} S3 method for the \code{survFitTKTD} class.
#' It prints the underlying JAGS model and some information on the Bayesian
#' inference procedure.
#'
#' @param x An object of class \code{gm_ssurvFitTKTD}
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords print
#'
#' @export
print.gm_survFitTKTD <- function(x, ...) {
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
