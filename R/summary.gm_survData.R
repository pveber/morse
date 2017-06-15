#' Summary for \code{gm_survData} objects
#'
#' The generic \code{summary} S3 method for the \code{gm_survData} class provides
#' information about the structure of the dataset and the experimental design.
#'
#' @param object an object of class \code{gm_survData}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns a list with the following fields:
#' \item{NbrepTimeConc}{nb of replicates for all concentrations and time points}
#' \item{NbsurvTimeConc}{nb of surviving ind. for all concentrations and time points}
#'
#'
#' @keywords summary
#'
#' @export
#'
summary.gm_survData <- function(object, quiet = FALSE, ...) {

  object = filter( object , !is.na(Nsurv))

  # matrix of number of survivors (sum of all replicates) by time / profile
  ans <- tapply(object$Nsurv, list(as.factor(object$profile),
                                    as.factor(object$time)), sum)

  if (! quiet) {
    cat("\nNumber of survivors (sum of replicates) per 'time' and 'profile': \n")
    print(ans)
  }

  invisible(list(NbsurvTimeProfile = ans))
}
