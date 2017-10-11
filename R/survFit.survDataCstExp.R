#' Fits a TKTD model for survival analysis using Bayesian inference for \code{survDataCstExp} object
#'
#' This function estimates the parameters of a TKTD
#' model for survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals is modeled as a function of the contaminant
#' concentration with a mechanistic description of toxic effects on survival over
#' time.
#'
#' Details of the model are presented in the vignette accompanying the package.
#'
#' @param data An object of class \code{survDataCstExp}.
#' @param model_type can be \code{"SD"} or \code{"IT"} to choose
#'   between "Stochastic Death" or "Individual Tolerance" models
#'   (resp.). See modeling vignette for details.
#' @param quiet If \code{FALSE}, prints logs and progress bar from
#'   JAGS.
#' @param nbr.chain Number of MCMC chains. The minimum required number 
#' of chains is 2.
#' @param nbr.adapt the number of iterations for adaptation. If \code{nbr.adapt} = 0
#'  then no adaptation takes place.
#' @param nbr.iter number of iterations to monitor
#' @param nbr.warmup 
#' @param thin.interval thinning interval for monitors
#' @param limit.sampling if \code{FALSE} (default is \code{TRUE}), there is no limit to the number of iterations
#' in MCMC imposed by the \code{diaftery.diag} test.
#' @param dic.compute if \code{TRUE} (default is \code{FALSE}), it generate penalized deviance samples to compute
#' the Deviance Information Criterion (DIC) with the \code{rjags} package
#' @param dic.type type of penalty to use. A string identifying the type of penalty: “pD” or “popt”
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
#' @examples
#'
#' # (1) Load the survival data
#' data(propiconazole)
#'
#' # (2) Create an object of class "survData"
#' dataset  <- survData(propiconazole)
#'
#' \dontrun{
#' # (3) Run the survFit function with TK-TD model 'SD' or 'IT' 
#' out <- survFit(dataset , model_type = "SD")
#'
#' # (4) Summary look the estimated values (parameters)
#' summary(out)
#'
#' # (5) Plot the fitted curve
#' plot(out, adddata = TRUE)
#'
#' # (6) Plot the fitted curve with ggplot style and CI as spaghetti
#' plot(out, spaghetti = TRUE , adddata = TRUE,
#'      style = "ggplot")
#' }
#'
#' @export
#' @import rjags
#' @importFrom dplyr group_by summarise filter
#'
survFit.survDataCstExp <- function(data,
                                   model_type = NULL,
                                   quiet = FALSE,
                                   nbr.chain = 3,
                                   nbr.adapt = 3000,
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
  
  globalData <- modelData(data,  model_type = model_type)
  
  jags.data <- unlist(list(globalData$dataList, globalData$priorsList), recursive = FALSE)
  
  jags.data_fit <- jags.data ; jags.data_fit$replicate = NULL
  
  priorsMinMax <- globalData$priorsMinMax

  ##
  ## Define model
  ##

  if(model_type == "SD"){
    ### Determine sampling parameters
    parameters_sampling <- c("kd_log10", "hb_log10", "kk_log10", "z_log10")
    parameters <- c("kd_log10", "hb_log10", "kk_log10", "z_log10", "psurv")

    file_to_use <- jags_TKTD_cstSD

  } else if(model_type == "IT"){
    ### Determine sampling parameters
    parameters_sampling <- c("kd_log10", "hb_log10", "alpha_log10", "beta_log10")
    
    parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10", "psurv")

    file_to_use <- jags_TKTD_cstIT
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
    
    if (sampling.parameters$niter > 2e5 & limit.sampling == TRUE){
      stop("The model needs too many iterations to provide reliable parameter estimates !")
    }
      
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
                       thin = thin.interval,
                       progress.bar = ifelse(quiet, "none", "text"))

  ##
  ## Cheking posterior range with data from experimental design:
  ##

  warnings <- msgTableCreate()

  estim.par <- survFit_TKTD_params(mcmc, model_type = model_type)

  if (filter(estim.par, parameters == "kd")$Q97.5 > priorsMinMax$kd_max){
    ##store warning in warnings table
    msg <- "The estimation of the dominant rate constant (model parameter kd)
    lies outside the range used to define its prior distribution which indicates
    that this rate is very high and difficult to estimate from this experiment !"
    warnings <- msgTableAdd(warnings, "kd_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  if (filter(estim.par, parameters == "hb")$Q2.5 < priorsMinMax$hb_min){
    ##store warning in warnings table
    msg <- "The estimation of the natural instantaneous mortality rate
    (model parameter hb) lies outside the range used to define its prior
    distribution which indicates that this rate is very low and so difficult
    to estimate from this experiment !"
    warnings <- msgTableAdd(warnings, "hb_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }

  ### for SD model
  if(model_type == "SD"){
    if (filter(estim.par, parameters == "kk")$Q97.5 > priorsMinMax$kk_max){
      ##store warning in warnings table
      msg <- "The estimation of the killing rate (model parameter kk) lies
      outside the range used to define its prior distribution which indicates
      that this rate is very high and difficult to estimate from this experiment !"
      warnings <- msgTableAdd(warnings, "kk_outRange", msg)
      ## print the message
      warning(msg, call. = FALSE)
    }

    if (filter(estim.par, parameters == "z")$Q2.5 < priorsMinMax$conc_min ||
        filter(estim.par, parameters == "z")$Q97.5 > priorsMinMax$conc_max){
      ##store warning in warnings table
      msg <- "The estimation of Non Effect Concentration threshold (NEC)
      (model parameter z) lies outside the range of tested concentration and
      may be unreliable as the prior distribution on this parameter
      is defined from this range !"
      warnings <- msgTableAdd(warnings, "z_outRange", msg)
      ## print the message
      warning(msg, call. = FALSE)
    }

  }

  ### for IT model
  if(model_type == "IT"){

    if (filter(estim.par, parameters == "alpha")$Q2.5 < priorsMinMax$conc_min ||
        filter(estim.par, parameters == "alpha")$Q97.5 > priorsMinMax$conc_max){
      ##store warning in warnings table
      msg <- "The estimation of log-logistic median (model parameter alpha) lies
      outside the range of tested concentration and may be unreliable as the prior
      distribution on this parameter is defined from this range !"
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

  class(OUT) <- c("survFitCstExp", "survFit")
  return(OUT)
}
