#' Fits a TKTD model for survival analysis using Bayesian inference for \code{survDataVarExp} object
#'
#' This function estimates the parameters of a TKTD ('SD' or 'IT')
#' model for survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals is modeled as a function of the chemical compound
#' concentration with a mechanistic description of the effects on survival over
#' time.
#'
#' The function \code{survFit} return the parameter estimates of Toxicokinetic-toxicodynamic (TK-TD) models
#' \code{SD} for 'Stochastic Death' or \code{IT} fo 'Individual Tolerance'.
#' TK-TD models, and particularly the General Unified Threshold model of
#' Survival (GUTS), provide a consistent process-based
#' framework to analyse both time and concentration dependent datasets.
#' In GUTS-SD, all organisms are assumed to have the same internal concentration 
#' threshold (denoted \eqn{z}), and, once exceeded, the instantaneous probability
#' to die increases linearly with the internal concentration.
#' In GUTS-IT, the threshold concentration is distributed among all the organisms, and once 
#' exceeded in one individual, this individual dies immediately.
#'
#' @param data An object of class \code{survDataVarExp}.
#' @param model_type can be \code{"SD"} or \code{"IT"} to choose
#'   between "Stochastic Death" or "Individual Tolerance" models
#'   (resp.). See modeling vignette for details.
#' @param quiet If \code{FALSE}, prints logs and progress bar from
#'   JAGS.
#' @param extend_time Number of for each replicate used for linear 
#' interpolation (comprise between time to compute and fitting accuracy)
#' @param n.chains A positive integer specifying the number of MCMC chains. The minimum required number 
#' of chains is 2.
#' @param n.adapt A positive integer specifying the number of iterations for adaptation. If \code{n.adapt} = 0
#'  then no adaptation takes place.
#' @param n.iter A positive integer specifying the number of iterations to monitor for each chain.
#' @param n.warmup A positive integer specifying the number of warmup (aka burnin) iterations per chain. 
#' @param thin.interval A positive integer specifying the period to monitor.
#' @param limit.sampling if \code{FALSE} (default is \code{TRUE}), there is no limit to the number of iterations
#' in MCMC imposed by the \code{raftery.diag} test.
#' @param dic.compute if \code{TRUE} (default is \code{FALSE}), it generates penalized deviance samples to compute
#' the Deviance Information Criterion (DIC) with the \code{rjags} package
#' @param dic.type type of penalty to use. A string identifying the type of penalty: \code{pD} or \code{popt}
#'  (see function \code{\link[rjags]{dic.samples}})
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns an object of class \code{survFitVarExp}, which is
#' a list with the following information:
#' \item{estim.par}{a table of the estimated parameters as medians and 95\%
#' credible intervals}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior
#' distribution}
#' \item{model}{a JAGS model object}
#' \item{dic}{return the Deviance Information Criterion (DIC) if \code{dic.compute} is \code{TRUE}}
#' \item{warnings}{a table with warning messages}
#' \item{parameters}{a list of parameter names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{mcmcInfo}{a table with the number of iterations, chains, adaptation, warmup and the thinning interval} 
#' \item{jags.data}{a list of the data passed to the JAGS model}
#' \item{model_type}{the type of TK-TD model used: \code{SD} or \code{IT}}
#'
#' @references Jager, T., Albert, C., Preuss, T. G. and Ashauer, R. (2011) 
#' General unified threshold model of survival-a toxicokinetic-toxicodynamic
#'  framework for ecotoxicology, \emph{Environmental Science and Technology}, 45, 2529-2540.
#' 303-314.
#' 
#' @keywords estimation
#' 
#' 
#' @examples
#'
#' # (1) Load the survival data
#' data("propiconazole_pulse_exposure")
#'
#' # (2) Create an object of class "survData"
#' dataset <- survData(propiconazole_pulse_exposure)
#'
#' \dontrun{
#' # (3) Run the survFit function with TK-TD model 'SD' or 'IT' 
#' out <- survFit(dataset , model_type = "SD")
#'
#' # (4) Summarize look the estimated parameters
#' summary(out)
#'
#' # (5) Plot the fitted curve
#' plot(out, adddata = FALSE)
#'
#' # (6) Plot the fitted curve with ggplot style and CI as spaghetti
#' plot(out, spaghetti = TRUE)
#' }
#' 
#' @export
#' @import rjags
#'

survFit.survDataVarExp <- function(data,
                                 model_type = NULL,
                                 quiet = FALSE,
                                 extend_time = 100,
                                 n.chains = 3,
                                 n.adapt = 1000,
                                 n.iter = NULL,
                                 n.warmup = NULL,
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
  if (n.chains < 2) {
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
    parameters <- c("kd_log10", "hb_log10", "z_log10", "kk_log10", "psurv", "Nsurv_ppc", "Nsurv_sim")
    
    jags.data_fit$time = NULL # remove jags.data_fit$time for varSD model
    
    file_to_use <- jags_TKTD_varSD
      
  } else if(model_type == "IT"){
    ### Determine sampling parameters
    parameters_sampling <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10")
    parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10", "psurv", "Nsurv_ppc", "Nsurv_sim")
    
    file_to_use <- jags_TKTD_varIT

  }

  model <- survLoadModel(model.program = file_to_use,
                         data = jags.data_fit,
                         n.chains = n.chains,
                         Nadapt = n.adapt,
                         quiet = quiet)

  
  ##
  ## estimate the number of iteration required for convergency of chains
  ## by using the raftery.diag
  ##
  
  if(is.null(n.warmup) | is.null(thin.interval) | is.null(n.iter)){
    
    sampling.parameters <- modelSamplingParameters(model,
                                                   parameters_sampling,
                                                   n.chains = n.chains, quiet = quiet)
    if (sampling.parameters$niter > 5e5)
      stop("The model needs too many iterations to provide reliable parameter estimates !")
    
    n.warmup = sampling.parameters$burnin
    thin.interval = sampling.parameters$thin
    n.iter = sampling.parameters$niter
    
  }

  ### model to check priors with the model
  update(model, n.warmup)
  
  if(dic.compute == TRUE){ # Deviance Information Criterion
    dic <- dic.samples(model,
                       n.iter = n.iter,
                       thin = thin.interval,
                       type = dic.type) 
  } else dic = NULL
  
  mcmc =  coda.samples(model,
                       variable.names = parameters,
                       n.iter = n.iter,
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
  mcmcInfo = data.frame(n.iter = n.iter,
                        n.chains = n.chains,
                        n.adapt = n.adapt,
                        thin.interval = thin.interval,
                        n.warmup = n.warmup)
  

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
