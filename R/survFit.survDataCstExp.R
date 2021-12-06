#' @rdname survFit
#'
#' @return The function returns an object of class \code{survFitCstExp}, which is
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
#' \item{mcmcInfo}{a table with the number of iterations, chains, adaptation, warmup and the thinning interval.} 
#' \item{jags.data}{a list of the data passed to the JAGS model}
#' \item{model_type}{the type of TKTD model used: \code{SD} or \code{IT}}
#'
#' @examples
#' 
#' # Example with time-variable exposure profile#'
#' # (1) Load the survival data
#' data(propiconazole)
#' # (2) Create an object of class "survData"
#' dataset  <- survData(propiconazole)
#' \dontrun{
#' # (3) Run the survFit function with TKTD model 'SD' or 'IT' 
#' out <- survFit(dataset , model_type = "SD")
#' # (4) Summarize look the estimated parameters
#' summary(out)
#' # (5) Plot the fitted curve
#' plot(out, adddata = TRUE)
#' # (6) Plot the fitted curve with ggplot style and CI as spaghetti
#' plot(out, spaghetti = TRUE , adddata = TRUE)
#' }
#'
#' @import rjags
#' @importFrom stats update
#' @importFrom dplyr group_by summarise filter
#'
#' @export
#'
survFit.survDataCstExp <- function(data,
                                   model_type = NULL,
                                   quiet = FALSE,
                                   n.chains = 3,
                                   n.adapt = 3000,
                                   n.iter = NULL,
                                   n.warmup = NULL,
                                   thin.interval = NULL,
                                   limit.sampling = TRUE,
                                   dic.compute = FALSE,
                                   dic.type = "pD",
                                   hb_value = TRUE,
                                   hb_valueFIXED = NA,
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
  ### warning message when hb_value = NULL
  if(hb_value==FALSE){
    warning("This is not an error message: the parameter 'hb' is fixed. This means that the correlation between
            'hb' and other parameters is ignored.")
    ## set default hb_valueFIXED
    if(is.na(hb_valueFIXED)){
      hb_valueFIXED = 0
    }
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
    if(hb_value == TRUE){
      jags.data_fit$hb_value = 1
      jags.data_fit$hb_valueFIXED = -1 # just to have it in JAGS
      parameters_sampling <- c("kd_log10", "hb_log10", "kk_log10", "z_log10")
      parameters <- c("kd_log10", "hb_log10", "kk_log10", "hb", "z_log10", "psurv", "Nsurv_ppc", "Nsurv_sim")
    } else{
      jags.data_fit$hb_value = 0
      jags.data_fit$hb_valueFIXED = hb_valueFIXED
      parameters_sampling <- c("kd_log10", "kk_log10", "z_log10")
      parameters <- c("kd_log10", "kk_log10", "z_log10", "hb", "psurv", "Nsurv_ppc", "Nsurv_sim")
    }
    file_to_use <- jags_TKTD_cstSD
    

  } else if(model_type == "IT"){
    ### Determine sampling parameters
    if(hb_value == TRUE){
      jags.data_fit$hb_value = 1
      jags.data_fit$hb_valueFIXED = -1 # just to have it in JAGS
      parameters_sampling <- c("kd_log10", "hb_log10", "alpha_log10", "beta_log10")
      parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10", "hb", "psurv", "Nsurv_ppc", "Nsurv_sim")
    } else{
      jags.data_fit$hb_value = 0
      jags.data_fit$hb_valueFIXED = hb_valueFIXED
      parameters_sampling <- c("kd_log10", "alpha_log10", "beta_log10")
      parameters <- c("kd_log10","alpha_log10", "beta_log10", "hb", "psurv", "Nsurv_ppc", "Nsurv_sim")
    }
    file_to_use <- jags_TKTD_cstIT
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
    
    if (sampling.parameters$niter > 2e5 & limit.sampling == TRUE){
      stop("The model needs too many iterations to provide reliable parameter estimates !")
    }
      
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
                       thin = thin.interval,
                       progress.bar = ifelse(quiet, "none", "text"))

  ##
  ## Cheking posterior range with data from experimental design:
  ##

  warnings <- msgTableCreate()

  estim.par <- survFit_TKTD_params(mcmc, model_type = model_type, hb_value = hb_value)

  if (filter(estim.par, parameters == "kd")$Q97.5 > priorsMinMax$kd_max){
    ##store warning in warnings table
    msg <- "The estimation of the dominant rate constant (model parameter kd)
    lies outside the range used to define its prior distribution which indicates
    that this rate is very high and difficult to estimate from this experiment !"
    warnings <- msgTableAdd(warnings, "kd_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  if(hb_value == TRUE){
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
              original.data = data,
              hb_valueFIXED = hb_valueFIXED)

  class(OUT) <- c("survFitCstExp", "survFit")
  return(OUT)
}
