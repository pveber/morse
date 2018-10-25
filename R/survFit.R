#' Method to fit a model for survival data using Bayesian inference
#'
#' @param data an object used to select a method 'survFit'
#' @param \dots Further arguments to be passed to generic methods

#' @export
survFit <- function(data, ...){
  UseMethod("survFit")
}


################################################################################
#
#  PRIORS
#
################################################################################

#' Create a list of scalars giving priors to use in Bayesian modelling
#'
#' @param x An object of class \code{survData}
#' @param model_type TKTD model type ('SD' or 'IT')
#' 
#' @return A list for parameterization of priors for Bayesian modeling with JAGS
#'
#' @examples 
#' 
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a survData object
#' dat <- survData(cadmium1)
#' 
#' # (3) Create priors for SD model_type
#' priors_survData(dat, model_type = "SD")
#' 
#' # (4) Create priors for IT model_type
#' priors_survData(dat, model_type = "IT")
#' 
#' @export


# priors <- function(x, ...){
#   UseMethod("priors")
# }


priors_survData <- function(x, model_type = NULL){
  
  data <- filter(x, time != 0)
  
  # Parameter calculation of concentration min and max
  conc_min <- min(data$conc[data$conc != 0], na.rm = TRUE) # to remove 0 and NA
  conc_max <- max(data$conc, na.rm = TRUE)
  
  time_min <- min(data$time)
  time_max <- max(data$time)
  
  conc_unic <- sort(unique(data$conc))
  conc_unicPrec <- dplyr::lag(conc_unic)
  conc_minDelta <- min(conc_unic - conc_unicPrec, na.rm = TRUE)
  
  ##
  ## dominant rate constant: kd
  ##
  
  kd_max <- -log(0.001) / time_min
  kd_min <- -log(0.999) / time_max
  
  ##
  ## background hazard rate
  ##
  
  hb_max <- -log(0.5) / time_min
  hb_min <- -log(0.999) / time_max
  
  ##
  ## killing rate parameter: kk
  ##
  
  kk_max <- -log(0.001) / (time_min * conc_minDelta)
  kk_min <- -log(0.999) / (time_max * (conc_max - conc_min))
  
  ##
  ## beta
  ##
  
  beta_minlog10 <- -2
  beta_maxlog10 <- 2
  
  priorsMinMax <- list(
    conc_min = conc_min,
    conc_max = conc_max,
    
    kd_min = kd_min,
    kd_max = kd_max,
    
    hb_min = hb_min,
    hb_max = hb_max )
  
  ##
  ## Construction of the list of priors
  ##
  
  priorsList <-  list(
    ##
    ## dominant rate constant: kd
    ##
    kd_meanlog10 = (log10(kd_max) + log10(kd_min)) / 2 ,
    kd_sdlog10 = (log10(kd_max) - log10(kd_min)) / 4 ,
    ##
    ## background hazard rate
    ##
    hb_meanlog10 = (log10(hb_max) + log10(hb_min)) / 2 ,
    hb_sdlog10 = (log10(hb_max) - log10(hb_min)) / 4
  )
  
  if(model_type == "IT"){
    
    ## priorsMinMax
    priorsMinMax$beta_min <- beta_minlog10
    priorsMinMax$beta_max <- beta_maxlog10
    
    ## priorsList
    ### non effect threshold: scale parameter & median of a log-logistic distribution
    priorsList$alpha_meanlog10 <- (log10(conc_max) + log10(conc_min)) / 2
    priorsList$alpha_sdlog10 <- (log10(conc_max) - log10(conc_min)) / 4
    
    ### shape parameter of a log-logistic distribution
    priorsList$beta_minlog10 <- beta_minlog10
    priorsList$beta_maxlog10 <- beta_maxlog10
    
  } else if (model_type == "SD"){
    
    ## priorsMinMax
    priorsMinMax$kk_min <- kk_min
    priorsMinMax$kk_max <- kk_max
    
    ## priorsList
    ### killing rate parameter: kk
    priorsList$kk_meanlog10 <- (log10(kk_max) + log10(kk_min)) / 2
    priorsList$kk_sdlog10 <- (log10(kk_max) - log10(kk_min)) / 4
    ### non effect threshold: z
    priorsList$z_meanlog10 <- (log10(conc_max) + log10(conc_min)) / 2
    priorsList$z_sdlog10 <- (log10(conc_max) - log10(conc_min)) / 4
  } else stop("please, provide the 'model_type': 'IT' or 'SD'")
  
  
  return(list(priorsList = priorsList,
              priorsMinMax = priorsMinMax))
}


#############################################################################
#
#    survFit_TKTD_params
#
#############################################################################
  
survFit_TKTD_params <- function(mcmc, model_type, hb_value = TRUE) {
    # create the table of posterior estimated parameters
    # for the survival analyses
    # INPUT:
    # - mcmc:  list of estimated parameters for the model with each item representing
    # a chain
    # OUTPUT:
    # - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated
    # parameters)
    
    # Retrieving parameters of the model
    res.M <- summary(mcmc)
    
    kd <- 10^res.M$quantiles["kd_log10", "50%"]
    kd_inf95 <- 10^res.M$quantiles["kd_log10", "2.5%"]
    kd_sup95 <- 10^res.M$quantiles["kd_log10", "97.5%"]
    
    if(hb_value == TRUE){
      hb <- 10^res.M$quantiles["hb_log10", "50%"]
      hb_inf95 <- 10^res.M$quantiles["hb_log10", "2.5%"]
      hb_sup95 <- 10^res.M$quantiles["hb_log10", "97.5%"]
    }
    
    if(model_type == "SD"){
      kk <- 10^res.M$quantiles["kk_log10", "50%"]
      kk_inf95 <- 10^res.M$quantiles["kk_log10", "2.5%"]
      kk_sup95 <- 10^res.M$quantiles["kk_log10", "97.5%"]
      
      z <- 10^res.M$quantiles["z_log10", "50%"]
      z_inf95 <- 10^res.M$quantiles["z_log10", "2.5%"]
      z_sup95 <- 10^res.M$quantiles["z_log10", "97.5%"]
      
      if(hb_value == TRUE){
        res <- data.frame(parameters = c("kd", "hb", "z", "kk"),
                          median = c(kd, hb, z, kk),
                          Q2.5 = c(kd_inf95, hb_inf95, z_inf95, kk_inf95),
                          Q97.5 = c(kd_sup95, hb_sup95, z_sup95, kk_sup95))
      } else{
        res <- data.frame(parameters = c("kd", "z", "kk"),
                          median = c(kd, z, kk),
                          Q2.5 = c(kd_inf95, z_inf95, kk_inf95),
                          Q97.5 = c(kd_sup95, z_sup95, kk_sup95))
      }
      
    } else if (model_type == "IT"){
      alpha <- 10^res.M$quantiles["alpha_log10", "50%"]
      alpha_inf95 <- 10^res.M$quantiles["alpha_log10", "2.5%"]
      alpha_sup95 <- 10^res.M$quantiles["alpha_log10", "97.5%"]
      
      beta <- 10^res.M$quantiles["beta_log10", "50%"]
      beta_inf95 <- 10^res.M$quantiles["beta_log10", "2.5%"]
      beta_sup95 <- 10^res.M$quantiles["beta_log10", "97.5%"]
      
      if(hb_value == TRUE){
        res <- data.frame(parameters = c("kd", "hb", "alpha", "beta"),
                          median = c(kd, hb, alpha, beta),
                          Q2.5 = c(kd_inf95, hb_inf95, alpha_inf95, beta_inf95),
                          Q97.5 = c(kd_sup95, hb_sup95, alpha_sup95, beta_sup95))
      } else{
        res <- data.frame(parameters = c("kd", "alpha", "beta"),
                          median = c(kd, alpha, beta),
                          Q2.5 = c(kd_inf95, alpha_inf95, beta_inf95),
                          Q97.5 = c(kd_sup95, alpha_sup95, beta_sup95))
      }
    } else {
      stop("please, provide the 'model_type': 'IT' or 'SD'")
    }
    
    return(res)
}
