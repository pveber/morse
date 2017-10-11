#' Summary of \code{survFitTKTD} object
#'
#' This is the generic \code{summary} S3 methode for the \code{survFitTKTD} class.
#' It shows the quantiles of priors and posteriors on parameters.
#'
#' @param object an object of class \code{survFitTKTD}
#' @param quiet when \code{TRUE}, does not print
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following information:
#' \item{Qpriors}{quantiles of the model priors}
#' \item{Qposteriors}{quantiles of the model posteriors}
#'
#' @examples
#' # (1) Load the data
#' data(propiconazole)
#'
#' # (2) Create a survData object
#' dat <- survData(propiconazole)
#'
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#'
#' # (4) summarize the survFitTKTD object
#' summary(out)
#' }
#'
#' @keywords summary
#'
#' @importFrom stats qnorm qunif
#' 
#' @export
summary.survFitTKTD <- function(object, quiet = FALSE, ...) {
  
  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start
  
  # kd
  log10kd <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = object$jags.data$meanlog10kd,
                   sd = 1 / sqrt(object$jags.data$taulog10kd))
  
  kd <- 10^log10kd
  
  # ks
  log10ks <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = object$jags.data$meanlog10ks,
                   sd = 1 / sqrt(object$jags.data$taulog10ks))
  
  ks <- 10^log10ks
  
  # nec
  log10nec <- qnorm(p = c(0.5, 0.025, 0.975),
                    mean = object$jags.data$meanlog10nec,
                    sd = 1 / sqrt(object$jags.data$taulog10nec))
  
  nec <- 10^log10nec
  
  # m0
  log10m0 <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = object$jags.data$meanlog10m0,
                   sd = 1 / sqrt(object$jags.data$taulog10m0))
  
  m0 <- 10^log10m0
  
  res <- rbind(kd, ks, nec, m0)
  
  ans1 <- format(data.frame(res), scientific = TRUE, digits = 4)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")
  
  # quantiles of estimated model parameters
  ans2 <- format(object$estim.par, scientific = TRUE, digits = 4)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")
  
  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1, row.names = FALSE)
    cat("\nPosteriors of the parameters (quantiles):\n\n")
    print(ans2, row.names = FALSE)
  }
  
  invisible(list(Qpriors = ans1,
                 Qpost = ans2))
}
