#' Summary for \code{survDataCstC} objects
#' 
#' The generic \code{summary} S3 method for the \code{survDataCst} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param object an object of class \code{survData}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function returns a list with the following fields:
#' \item{NbrepTimeConc}{nb of replicates for all concentrations and time points}
#' \item{NbsurvTimeConc}{nb of surviving ind. for all concentrations and time points}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#' 
#' # (3) Summarize the dataset
#' summary(cadmium1)
#' 
#' @keywords summary
#' 
#' @export
summary.survData <- function(object, quiet = FALSE, ...) {
  # matrix of number of replicate by time / conc
  ans1 <- table(object[, c("conc", "time")])
  
  # matrix of number of survival (sum of all replicate) by time / conc
  ans2 <- tapply(object$Nsurv, list(as.factor(object$conc),
                                    as.factor(object$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of replicates per time and concentration: \n")
    print(ans1)
    cat("\nNumber of survivors (sum of replicate) per time and concentration: \n")
    print(ans2)
  }
  
  invisible(list(NbrepTimeConc = ans1,
                 NbsurvTimeConc = ans2))
}

#' Summary for \code{survDataVarC} objects
#' 
#' The generic \code{summary} S3 method for the \code{survDataVarC} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param object an object of class \code{survData}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function returns a list with the following fields:
#' \item{NsurvReplicateTime}{Number of survivors per replicate and time}
#' \item{concReplicateTime}{Concentration per replicate and time}
#' 
#' @examples
#' # (1) Load the data
#' data(propiconazole_pulse_exposure)
#' 
#' # (2) Create a survData object
#' propiconazole_pulse_exposure <- survData(propiconazole_pulse_exposure)
#' 
#' # (3) Summarize the dataset
#' summary(propiconazole_pulse_exposure)
#' 
#' @keywords summary
#' 
#' @export
summary.survDataVarC <- function(object, quiet = FALSE, ...) {
  
  df_ans1 = dat_var %>%
    select(-conc) %>%
    spread(time, Nsurv)
  
  df_ans2 = dat_var %>%
    select(-Nsurv) %>%
    spread(time, conc)
  
  if (! quiet) {
    cat("\nNumber of survivors per replicate and time: \n")
    print(df_ans1)
    cat("\nConcentration per replicate and time: \n")
    print(df_ans2)
  }
  
  invisible(list(NbrepTimeConc = df_ans1,
                 NbsurvTimeConc = df_ans2))
}

