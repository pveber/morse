#' Print of \code{survFit} object
#' 
#' This is the generic \code{print} S3 method for the \code{survFit} class.
#' It prints the underlying JAGS model and some information on the Bayesian 
#' inference procedure.
#' 
#' @param x An object of class \code{survFit}
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @examples
#' # (1) Load the data
#' data(propiconazole)
#' 
#' # (2) Create a survData object
#' dat <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFit function
#' out <- survFit(dat, quiet = TRUE, model_type="SD")
#' 
#' # (4) Print the survFit object
#' print(out)
#' }
#' 
#' @keywords print
#' 
#' @export
print.survFit <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class survFit
  
  summary_mcmc <- summary(x$mcmc)
  
  n.chains <- summary_mcmc$nchain
  n.thin <- summary_mcmc$thin
  end <- summary_mcmc$end
  start <- summary_mcmc$start

  
  # M.C.M.C. informations
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Iterations = ", start, ":",
      end, "\n", sep = "")
  cat("Thinning interval =", n.thin, "\n")
  cat("Number of chains =", n.chains, "\n")
  cat("Sample size per chain =",
      (end - start) / n.thin + 1, "\n")
}
