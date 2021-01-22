#' Summary of \code{survFit} object
#'
#' This is the generic \code{summary} S3 method for the \code{survFit} class.
#' It shows the quantiles of priors and posteriors on parameters.
#'
#' @param object An object of class \code{survFit}.
#' @param quiet When \code{TRUE}, does not print.
#' @param EFSA_name If \code{TRUE}, the current terminology by
#'  the one used in the recent EFSA PPR Scientific Opinion (2018).
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following information:
#' \item{Qpriors}{quantiles of the model priors}
#' \item{Qposteriors}{quantiles of the model posteriors}
#' 
#' @references 
#' EFSA PPR Scientific Opinion (2018)
#' \emph{Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms}
#' \url{https://www.efsa.europa.eu/en/efsajournal/pub/5377}.
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

summary.survFit <- function(object,
                            quiet = FALSE,
                            EFSA_name = FALSE,
                            ...) {
  
  estim_parameters <- object$estim.par
  
  if(EFSA_name == TRUE){
    parameters_SD_HBon = c("kD", "hb", "zw", "bw")
    parameters_SD_HBoff = c("kD", "zw", "bw")
    parameters_IT_HBon = c("kD", "hb", "mw", "beta")
    parameters_IT_HBoff = c("kD", "mw", "beta")
    
    estim_parameters$parameters <- gsub("kd","kD", estim_parameters$parameters)
    estim_parameters$parameters <- gsub("kk","bw", estim_parameters$parameters)
    estim_parameters$parameters <- gsub("z","zw", estim_parameters$parameters)
    estim_parameters$parameters <- gsub("alpha","mw", estim_parameters$parameters)
    
  } else{
    parameters_SD_HBon = c("kd", "hb", "z", "kk")
    parameters_SD_HBoff = c("kd", "z", "kk")
    parameters_IT_HBon = c("kd", "hb", "alpha", "beta")
    parameters_IT_HBoff = c("kd", "alpha", "beta")
  }
  
  param <- object$jags.data
  if("hb" %in% estim_parameters[, "parameters"]){
    hb_value = TRUE
  } else{
    hb_value = FALSE
  }
  
  # kd
  kd_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                   mean = param$kd_meanlog10,
                   sd = param$kd_sdlog10)
  
  kd <- 10^kd_log10
  
  
  # hb
  if(hb_value == TRUE){
    hb_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                      mean = param$hb_meanlog10,
                      sd = param$hb_sdlog10)
    
    hb <- 10^hb_log10
  } 
  
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
    
    if(hb_value == TRUE){
      res <- data.frame(parameters = parameters_SD_HBon,
                        median = c(kd[1], hb[1], z[1], kk[1]),
                        Q2.5 = c(kd[2], hb[2], z[2], kk[2]),
                        Q97.5 = c(kd[3], hb[3], z[3], kk[3]))
    } else{
      res <- data.frame(parameters = parameters_SD_HBoff,
                        median = c(kd[1], z[1], kk[1]),
                        Q2.5 = c(kd[2], z[2], kk[2]),
                        Q97.5 = c(kd[3], z[3], kk[3]))
    }
    
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
    
    if(hb_value == TRUE){
      res <- data.frame(parameters = parameters_IT_HBon,
                        median = c(kd[1], hb[1], alpha[1], beta[1]),
                        Q2.5 = c(kd[2], hb[2], alpha[2], beta[2]),
                        Q97.5 = c(kd[3], hb[3], alpha[3], beta[3]))
    } else{
      res <- data.frame(parameters = parameters_IT_HBoff,
                        median = c(kd[1], alpha[1], beta[1]),
                        Q2.5 = c(kd[2], alpha[2], beta[2]),
                        Q97.5 = c(kd[3], alpha[3], beta[3]))
    }
    
  }
  
  
  ans1 <- format(data.frame(res), scientific = TRUE, digits = 4)
  
  # quantiles of estimated model parameters
  ans2 <- format(estim_parameters, scientific = TRUE, digits = 4)
  
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
