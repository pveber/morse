#' Summary for survFitTKTD objects
#'
#' The summary shows the quantiles of priors and posteriors on parameters
#'
#' @param object an object of class \code{survFitTKTD}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following fields:
#' \item{Qpriors}{quantiles for the model's prior}
#' \item{Qposteriors}{quantiles for the model's posteriors}
#'
#' @seealso survFitTKTD
#'
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTT(dat, distr = "norm")
#'
#' # (4) summarize the survFitTT object
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
  
  if (object$distr == "norm") {
    # ke
    log10ke <- qnorm(p = c(0.5, 0.025, 0.975),
                     mean = object$jags.data$meanlog10ke,
                     sd = 1 / sqrt(object$jags.data$taulog10ke))
    
    ke <- 10^log10ke
    
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
    
  } else if (object$distr == "unif") {
    # ke
    log10ke <- qunif(p = c(0.5, 0.025, 0.975),
                     min = object$jags.data$kemin,
                     max = object$jags.data$kemax)
    
    ke <- 10^log10ke
    
    # ks
    log10ks <- qunif(p = c(0.5, 0.025, 0.975),
                     min = object$jags.data$ksmin,
                     max = object$jags.data$ksmax)
    
    ks <- 10^log10ks
    
    # nec
    log10nec <- qunif(p = c(0.5, 0.025, 0.975),
                      min = object$jags.data$concmin,
                      max = object$jags.data$concmax)
    
    nec <- 10^log10nec
    
    # m0
    log10m0 <-qunif(p = c(0.5, 0.025, 0.975),
                    min = object$jags.data$m0min,
                    max = object$jags.data$m0max)
    
    m0 <- 10^log10m0
  }
  
  res <- rbind(ke, ks, nec, m0)
  
  ans1 <- round(data.frame(res), digits = 3)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")
  
  # quantiles of estimated model parameters
  ans2 <- round(object$estim.par, digits = 3)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")
  
  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("The ", object$distr, " model was used !\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosterior of the parameters (quantiles):\n\n")
    print(ans2)
  }
  
  invisible(list(Qpriors = ans1,
                 Qpost = ans2))
}
