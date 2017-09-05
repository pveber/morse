#' Fits a TKTD model for survival analysis using Bayesian inference
#'
#' This function estimates the parameters of a TKTD
#' model for survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals is modeled as a function of the pollutant's
#' concentration with a mechanistic description of toxic effects on survival over
#' time.
#'
#' Details of the model are presented in the vignette accompanying the package.
#'
#' @param data An object of class \code{survData}.
#' @param model_type can be \code{"SD"} or \code{"IT"} to choose
#'   between "Stochastic Death" or "Individual Tolerance" models
#'   (resp.). See modeling vignette for details.
#' @param n.chains Number of MCMC chains. The minimum required number
#'   of chains is 2.
#' @param quiet If \code{FALSE}, prints logs and progress bar from
#'   JAGS.
#'
#' @return The function returns an object of class \code{survFitCstExp}, which is
#' a list with the following fields:
#' \item{estim.par}{a table of the estimated parameters (medians) and 95 \%
#' credible intervals}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior
#' distributions}
#' \item{model}{a JAGS model object}
#' \item{dic}{return the Deviance Information Criterion (DIC) if \code{dic.compute} is \code{TRUE}}
#' \item{warnings}{a data.frame with warning messages}
#' \item{parameters}{a list of the parameters names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{mcmcInfo}{a data.frame with the number of iteration, chains, adaptation, warmup and the thinning interval.} 
#' \item{jags.data}{a list a the data passed to the jags model}
#'
#' @keywords estimation
#
#' @export
#' @import rjags
#'

survFit.survDataVarExp <- function(data,
                                 model_type = NULL,
                                 quiet = FALSE,
                                 extend_time = 100,
                                 nbr.chain = 3,
                                 nbr.adapt = 1000,
                                 nbr.iter = NULL,
                                 nbr.warmup = NULL,
                                 thin.interval = NULL,
                                 limit.sampling = TRUE,
                                 dic.compute = FALSE,
                                 dic.type = "pD",
                                 ...){
  
  ##
  ## Pre modelling measure and tests
  ##

  ### ensures model_type is one of "SD" and "IT"
  if(is.null(model_type) || ! (model_type %in% c("SD","IT"))) {
    stop("You need to specify a 'model_type' among 'SD' or 'IT'")
  }
  
  ### check number of sample for the diagnostic procedure
  if (nbr.chain < 2) {
    stop('2 or more parallel chains required')
  }
  
  ##
  ## Data and Priors for model
  ##
  
  globalData <- modelData(x = data, model_type = model_type, extend_time = extend_time)
  
  ### Remove the information of replicate since this is not used in JAGS, and so a warning message would be show
  
  jags.data <- globalData$modelData
  
  jags.data_fit <- jags.data
  
  jags.data_fit$replicate <- NULL
  jags.data_fit$conc <- NULL
  jags.data_fit$replicate_long <- NULL
  
  priorsData = globalData$priorsMinMax
  
  ##
  ## Define model
  ##
  
  if(model_type == "SD"){
    ### Determine sampling parameters
    parameters_sampling <- c("kd_log10", "hb_log10", "z_log10", "kk_log10")
    parameters <- c("kd_log10", "hb_log10", "z_log10", "kk_log10", "psurv", "Nsurv_ppc")
    
    jags.data_fit$time = NULL # remove jags.data_fit$time for varSD model
    
    file_to_use <- jags_TKTD_varSD
      
  } else if(model_type == "IT"){
    ### Determine sampling parameters
    parameters_sampling <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10")
    parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10", "psurv", "Nsurv_ppc")
    
    file_to_use <- jags_TKTD_varIT
  }

  model <- survLoadModel(model.program = file_to_use,
                         data = jags.data_fit,
                         n.chains = nbr.chain,
                         Nadapt = nbr.adapt,
                         quiet = quiet)

  
  ##
  ## estimate the number of iteration required for convergency of chains
  ## by using the raftery.diag
  ##
  
  if(is.null(nbr.warmup) | is.null(thin.interval) | is.null(nbr.iter)){
    
    sampling.parameters <- modelSamplingParameters(model,
                                                   parameters_sampling,
                                                   n.chains = nbr.chain, quiet = quiet)
    if (sampling.parameters$niter > 5e5)
      stop("The model needs too many iterations to provide reliable parameter estimates !")
    
    nbr.warmup = sampling.parameters$burnin
    thin.interval = sampling.parameters$thin
    nbr.iter = sampling.parameters$niter
    
  }

  ### model to check priors with the model
  update(model, nbr.warmup)
  
  if(dic.compute == TRUE){ # Deviance Information Criterion
    dic <- dic.samples(model,
                       n.iter = nbr.iter,
                       thin = thin.interval,
                       type = dic.type) 
  } else dic = NULL
  
  mcmc =  coda.samples(model,
                       variable.names = parameters,
                       n.iter = nbr.iter,
                       thin = thin.interval)
  
  ##
  ## Cheking posterior range with data from experimental design:
  ##
  
  estim.par <- survFit_TKTD_params(mcmc, model_type = model_type)
  
  warnings <- msgTableCreate()
  
  if (filter(estim.par, parameters == "kd")$Q97.5 > priorsData$kd_max){
    ## store warning in warnings table
    msg <- "The estimation of the dominant rate constant (model parameter kd) lies 
    outside the range used to define its prior distribution which indicates that this
    rate is very high and difficult to estimate from this experiment !"
    warnings <- msgTableAdd(warnings, "kd_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  
  if (filter(estim.par, parameters == "hb")$Q2.5 < priorsData$hb_min){
    ## store warning in warnings table
    msg <- "The estimation of the natural instantaneous mortality rate (model 
    parameter hb) lies outside the range used to define its prior distribution 
    which indicates that this rate is very low and so difficult to estimate 
    from this experiment !"
    warnings <- msgTableAdd(warnings, "hb_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  
  ### for SD model
  if(model_type == "SD"){
    if (filter(estim.par, parameters == "kk")$Q97.5 > priorsData$kk_max){
      ## store warning in warnings table
      msg <- "The estimation of the killing rate (model parameter k) lies
      outside the range used to define its prior distribution which indicates
      that this rate is very high and difficult to estimate from this experiment !"
      warnings <- msgTableAdd(warnings, "kk_outRange", msg)
      ## print the message
      warning(msg, call. = FALSE)
    }
    
    if (filter(estim.par, parameters == "z")$Q2.5 < priorsData$conc_min ||
        filter(estim.par, parameters == "z")$Q97.5 > priorsData$conc_max){
      ## store warning in warnings table
      msg <- "The estimation of Non Effect Concentration threshold (NEC) 
      (model parameter z) lies outside the range of tested concentration 
      and may be unreliable as the prior distribution on this parameter is
      defined from this range !"
      warnings <- msgTableAdd(warnings, "z_outRange", msg)
      ## print the message
      warning(msg, call. = FALSE)
    }
  }
  
  ### for IT model
  if(model_type == "IT"){
    
    if (filter(estim.par, parameters == "alpha")$Q2.5 < priorsData$conc_min ||
        filter(estim.par, parameters == "alpha")$Q97.5 > priorsData$conc_max){
      ## store warning in warnings table
      msg <- "The estimation of log-logistic median (model parameter alpha) 
      lies outside the range of tested concentration and may be unreliable as 
      the prior distribution on this parameter is defined from this range !"
      warnings <- msgTableAdd(warnings, "alpha_outRange", msg)
      ## print the message
      warning(msg, call. = FALSE)
    }
  }
  
  ##
  ## MCMC information
  ## 
  mcmcInfo = data.frame(nbr.iter = nbr.iter,
                        nbr.chain = nbr.chain,
                        nbr.adapt = nbr.adapt,
                        thin.interval = thin.interval,
                        nbr.warmup = nbr.warmup)
  

  ##
  ##
  ##
  transformed.data <- data.frame(
    replicate = jags.data$replicate,
    time = jags.data$time,
    conc = jags.data$conc,
    Nsurv = jags.data$Nsurv
  ) %>%
    group_by(replicate) %>%
    mutate(Ninit = max(Nsurv, na.rm = TRUE))
  
  ##
  ## OUTPUT
  ##
  
  OUT <- list(estim.par = estim.par,
              mcmc = mcmc,
              model = model,
              dic = dic,
              parameters = parameters,
              mcmcInfo = mcmcInfo,
              jags.data = jags.data,
              warnings = warnings,
              model_type = model_type,
              transformed.data = transformed.data,
              original.data = data)
  

  class(OUT) <- c("survFitVarExp","survFit")
  return(OUT)
}
