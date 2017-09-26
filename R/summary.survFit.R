#' Summary for \code{survFit} objects
#'
#' This is the generic \code{summary} S3 methode for the \code{survFit} class.
#' It shows the quantiles of priors and posteriors on parameters.
#'
#' @param object an object of class \code{survFit}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following fields:
#' \item{Qpriors}{quantiles for the model's prior}
#' \item{Qposteriors}{quantiles for the model's posteriors}
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
#' out <- survFit(dat, model_type = "SD")
#'
#' # (4) summarize the survFit object
#' summary(out)
#' }
#'
#' @keywords summary
#'
#' @importFrom stats qnorm qunif
#' 
#' @export
#' 

summary.survFit <- function(object, quiet = FALSE, ...) {
  
  param <- object$jags.data
  
  # kd
  kd_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = param$kd_meanlog10,
                   sd = param$kd_sdlog10)
  
  kd <- 10^kd_log10
  
  
  # hb
  hb_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = param$hb_meanlog10,
                   sd = param$hb_sdlog10)
  
  hb <- 10^hb_log10
  
  if(object$model_type == "SD"){
    
    # kk
    kk_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                      mean = param$kk_meanlog10,
                      sd = param$kk_sdlog10)
    
    kk <- 10^kk_log10
    
    ## z
    z_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                     mean = param$z_meanlog10,
                     sd = param$z_sdlog10)
    
    z <- 10^z_log10
    
    res <- data.frame(parameters = c("kd", "hb", "z", "kk"),
                      median = c(kd[1], hb[1], z[1], kk[1]),
                      Q2.5 = c(kd[2], hb[2], z[2], kk[2]),
                      Q97.5 = c(kd[3], hb[3], z[3], kk[3]))
    
  }
  if(object$model_type == "IT"){
    
    # alpha
    alpha_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                     mean = param$alpha_meanlog10,
                     sd = param$alpha_sdlog10)
    
    alpha <- 10^alpha_log10
    
    # beta
    beta_log10 <- qunif(p = c(0.5, 0.025, 0.975),
                      min = param$beta_minlog10,
                      max = param$beta_maxlog10)
    
    beta <- 10^beta_log10
    
    res <- data.frame(parameters = c("kd", "hb", "alpha", "beta"),
                      median = c(kd[1], hb[1], alpha[1], beta[1]),
                      Q2.5 = c(kd[2], hb[2], alpha[2], beta[2]),
                      Q97.5 = c(kd[3], hb[3], alpha[3], beta[3]))
    
  }
  
  
  ans1 <- format(data.frame(res), scientific = TRUE, digits = 4)
  
  # quantiles of estimated model parameters
  ans2 <- format(object$estim.par, scientific = TRUE, digits = 4)
  
  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("Priors of the parameters (quantiles) (select with '$Qpriors'):\n\n")
    print(ans1, row.names = FALSE)
    cat("\nPosteriors of the parameters (quantiles) (select with '$Qposteriors'):\n\n")
    print(ans2, row.names = FALSE)
  }
  
  invisible(list(Qpriors = ans1,
                 Qposteriors = ans2))
}
