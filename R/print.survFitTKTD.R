#' Print of \code{survFitTKTD} object
#' 
#' This is the generic \code{print} S3 method for the \code{survFitTKTD} class.
#' It prints the underlying JAGS model and some information on the Bayesian 
#' inference procedure.
#' 
#' @param x An object of class \code{survFitTKTD}
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @examples
#' # (1) Load the data
#' data(propiconazole)
#' 
#' # (2) Create an object of class 'survData'
#' dat <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat, quiet = TRUE)
#' 
#' # (4) Print the survFitTKTD object
#' print(out)
#' }
#' 
#' @keywords print
#' 
#' @export
print.survFitTKTD <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class survFitTKTD
  
  # M.C.M.C. informations
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Iterations = ", x$n.iter[["start"]], ":",
      x$n.iter[["end"]], "\n", sep = "")
  cat("Thinning interval =", x$n.thin, "\n")
  cat("Number of chains =", x$n.chains, "\n")
  cat("Sample size per chain =",
      (x$n.iter[["end"]] - x$n.iter[["start"]]) / x$n.thin + 1, "\n")
}
