#' Summary of \code{survDataVarExp} object
#' 
#' The generic \code{summary} S3 method for the \code{survDataVarExp} class provides
#' information about the structure of the data set and the experimental design.
#' 
#' @param object an object of class \code{survDataVarExp}
#' @param quiet when \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function returns a list with the following information:
#' \item{OccRepTime}{Occurence of replicates for all time points}
#' \item{NbsurvTimeRep}{nb of survivors. for all replicates and time points}
#' \item{ConcTimeRep}{Concentration for all replicates and time points}
#' 
#' @examples
#' # (1) Load the data
#' data(propiconazole_pulse_exposure)
#' 
#' # (2) Create a survDataVarExp object
#' out <- survData(propiconazole_pulse_exposure)
#' 
#' # (3) Summarize the data set
#' summary(out)
#' 
#' @keywords summary
#' 
#' @export
summary.survDataVarExp <- function(object, quiet = FALSE, ...) {
  
  
  # matrix of number of replicate by time / replicate   
  ans1 <- table(object[, c("replicate", "time")])
  
  # matrix of number of survival by time / replicate
  ans2 <- tapply(object$Nsurv, list(as.factor(object$replicate),
                                    as.factor(object$time)), sum)
  
  # matrix of concentration by time / replicate
  ans3 <- tapply(object$conc, list(as.factor(object$replicate),
                                    as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nOccurence of 'replicate' for each 'time': \n")
    print(ans1)
    cat("\nNumber of survivors per 'time' and 'replicate': \n")
    print(ans2)
    cat("\nConcentrations per 'time' and 'replicate': \n")
    print(ans3)
  }
  
    invisible(list(OccRepTime = ans1,
                   NbsurvTimeRep = ans2,
                   ConcTimeRep = ans3))
}
