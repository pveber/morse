#' Density distribution of priors.
#' 
#' Return a \class{data.frame} with priors density distribution of parameters used in
#' \class{object}.
#' 
#' When the \code{object} is of class \code{survFit}, see \link[=priors_distribution.survFit]{priors_distribution.survFit}
#' 
#' @param object An object used to select a method
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#' 
priors_distribution <- function(object, ...){
  UseMethod("priors_distribution")
}

#' Density distribution of priors from a \class{survFit} object.
#' 
#' Return a \class{data.frame} with priors distribution of parameters used in
#' \class{object}.
#' 
#' @param object An object of class \code{survFit}.
#' @param size_sample Size of the random generation of the distribution.
#' Default is \code{1e3}.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @importFrom stats rnorm runif
#' 
#' @export
#' 

priors_distribution.survFit <- function(object, size_sample = 1e3, ...){
  
  param <- object$jags.data
  
  # kd
  kd_log10 <- rnorm(size_sample,
                    mean = param$kd_meanlog10,
                    sd = param$kd_sdlog10)
  
  kd <- 10^kd_log10
  
  
  # hb
  hb_log10 <- rnorm(size_sample,
                    mean = param$hb_meanlog10,
                    sd = param$hb_sdlog10)
  
  hb <- 10^hb_log10
  
  df = data.frame(hb = hb,
                  hb_log10 = hb_log10,
                  kd = kd,
                  kd_log10 = kd_log10)
  
  
  if(object$model_type == "SD"){
    
    # kk
    kk_log10 <- rnorm(size_sample,
                      mean = param$kk_meanlog10,
                      sd = param$kk_sdlog10)
    
    kk <- 10^kk_log10
    
    ## z
    z_log10 <- rnorm(size_sample,
                     mean = param$z_meanlog10,
                     sd = param$z_sdlog10)
    
    z <- 10^z_log10
    
    df$kk = kk
    df$kk_log10 = kk_log10
    df$z = z
    df$z_log10 = z_log10
    
  }
  if(object$model_type == "IT"){
    
    # alpha
    alpha_log10 <- rnorm(size_sample,
                         mean = param$alpha_meanlog10,
                         sd = param$alpha_sdlog10)
    
    alpha <- 10^alpha_log10
    
    # beta
    beta_log10 <- runif(size_sample,
                        min = param$beta_minlog10,
                        max = param$beta_maxlog10)
    
    beta <- 10^beta_log10

    df$alpha = alpha
    df$alpha_log10 = alpha_log10
    df$beta = beta
    df$beta_log10 = beta_log10
  }
  if(object$model_type == "PROPER"){
    
    # kk
    kk_log10 <- rnorm(size_sample,
                      mean = param$kk_meanlog10,
                      sd = param$kk_sdlog10)
    
    kk <- 10^kk_log10
    
    # alpha
    alpha_log10 <- rnorm(size_sample,
                         mean = param$alpha_meanlog10,
                         sd = param$alpha_sdlog10)
    
    alpha <- 10^alpha_log10
    
    # beta
    beta_log10 <- runif(size_sample,
                        min = param$beta_minlog10,
                        max = param$beta_maxlog10)
    
    beta <- 10^beta_log10
    
    df$kk = kk
    df$kk_log10 = kk_log10
    
    df$alpha = alpha
    df$alpha_log10 = alpha_log10
    df$beta = beta
    df$beta_log10 = beta_log10
  }

  return(df)
}
