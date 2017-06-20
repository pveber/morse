#' Fits a TKTD for survival analysis using Bayesian inference
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
#' @param n.chains Number of MCMC chains. The minimum required number of chains
#' is 2.
#' @param quiet If \code{FALSE}, prints logs and progress bar from JAGS.
#' 
#' @return The function returns an object of class \code{survFitCstC}, which is
#' a list with the following fields:
#' \item{estim.par}{a table of the estimated parameters (medians) and 95 \%
#' credible intervals}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior
#' distributions}
#' \item{model}{a JAGS model object}
#' \item{parameters}{a list of the parameters names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{n.iter}{a list of two indices indicating the beginning and end of
#' monitored iterations}
#' \item{n.thin}{a numerical value corresponding to the thinning interval}
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
#' dat <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFit function
#' out <- survFit(dat)
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

survFit.survDataCstC <- function(data,
                                model_type = NULL,
                                quiet = FALSE,
                                nbr.chain = 3,
                                nbr.adapt = 1000,
                                nbr.iter = NULL,
                                nbr.warmup = NULL,
                                nbr.thin = NULL){
  
  ##
  ## Pre modelling measure and tests
  ##
  
  ### time start
  time_start <- Sys.time()
  
  ### ensures model_type is one of "SD" and "IT"
  if(is.null(model_type)) {
    stop("You need to specify a 'model_type': 'SD' or 'IT'")
  }
  if(!model_type %in% c("SD","IT")) {
    stop("'model_type' available for use are 'SD' or 'IT'")
  }
  
  ### check number of sample for the diagnostic procedure
  if (nbr.chain < 2) {
    stop('2 or more parallel chains required')
  }
  
  ##
  ## Data and Priors for model
  ##
  
  globalData = modelData(data,
                         model_type = model_type)
  
  modelData_ = globalData$modelData
  modelData = modelData_ ; modelData$replicate = NULL
  
  modelData_Null_ = globalData$modelData_Null
  modelData_Null = modelData_Null_ ; modelData_Null$replicate = NULL
  
  priorsData = globalData$priorsMinMax
  
  ##
  ## Define model
  ##

  if(model_type == "SD"){
    ### Determine sampling parameters
    parameters_red <- c("kd_log10", "hb_log10", "kk_log10", "z_log10")
    parameters <- c("kd_log10", "hb_log10", "kk_log10", "z_log10", "psurv", "Nsurv_ppc", "Nsurv_sim")
    
    file_to_use <- jags_TKTD_cstSD

  } else if(model_type == "IT"){
    ### Determine sampling parameters
    parameters_red <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10")
    parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10", "psurv", "Nsurv_ppc", "Nsurv_sim")
    
    file_to_use <- jags_TKTD_cstIT
  }
  
  model_Null <- survLoadModel(model.program = file_to_use,
                              data = modelData_Null,
                              n.chains = nbr.chain,
                              Nadapt = nbr.adapt,
                              quiet = TRUE)
  
  
  model <- survLoadModel(model.program = file_to_use,
                         data = modelData,
                         n.chains = nbr.chain,
                         Nadapt = nbr.adapt,
                         quiet = quiet)
  
  
  ##
  ## estimate the number of iteration required for convergency of chains
  ## by using the raftery.diag
  ##

  if(is.null(nbr.warmup) | is.null(nbr.thin) | is.null(nbr.iter)){
  
    
    sampling.parameters <- modelSamplingParameters(model,
                                                   parameters_red,
                                                   n.chains = nbr.chain, quiet)
    if (sampling.parameters$niter > 5e5)
      stop("The model needs too many iterations to provide reliable parameter estimates !")
    
    nbr.warmup = sampling.parameters$burnin
    nbr.thin = sampling.parameters$thin
    nbr.iter = sampling.parameters$niter
    
  }
  
  ### Null model to check priors with the model
  update(model_Null, nbr.warmup)
  
  mcmc_Null =  coda.samples(model_Null,
                            variable.names = parameters_red,
                            n.iter = nbr.iter,
                            thin = nbr.thin)
  
  ### model to check priors with the model
  update(model, nbr.warmup)
  mcmc =  coda.samples(model,
                       variable.names = parameters,
                       n.iter = nbr.iter,
                       thin = nbr.thin)
  
  ##
  ## Cheking posterior range with data from experimental design:
  ##
  
  estim.par <- survTKTDPARAMS(mcmc, model_type = model_type)
  
  if (filter(estim.par, parameters == "kd")$Q97.5 > priorsData$kd_max){
    warning("The estimation of the dominant rate constant (model parameter kd) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
            call. = FALSE)
  }
  if (filter(estim.par, parameters == "hb")$Q2.5 < priorsData$hb_min){
    warning("The estimation of the natural instantaneous mortality rate (model parameter hb) lies outside the range used to define its prior distribution which indicates that this rate is very low and so difficult to estimate from this experiment !",
            call. = FALSE)
  }
  
  ### for SD model
  if(model_type == "SD"){
    if (filter(estim.par, parameters == "kk")$Q97.5 > priorsData$kk_max)
      warning("The estimation of the killing rate (model parameter k) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
              call. = FALSE)
    
    if (filter(estim.par, parameters == "z")$Q2.5 < priorsData$conc_min ||
        filter(estim.par, parameters == "z")$Q97.5 > priorsData$conc_max)
      warning("The estimation of Non Effect Concentration threshold (NEC) (model parameter z) lies outside the range of tested concentration and may be unreliable as the prior distribution on this parameter is defined from this range !",
              call. = FALSE)
    
  }
  
  ### for IT model
  if(model_type == "IT"){
    
    if (filter(estim.par, parameters == "alpha")$Q2.5 < priorsData$conc_min ||
        filter(estim.par, parameters == "alpha")$Q97.5 > priorsData$conc_max)
      warning("The estimation of log-logistic median (model parameter alpha) lies outside the range of tested concentration and may be unreliable as the prior distribution on this parameter is defined from this range !",
              call. = FALSE)
    
  }
  
  ### time end
  
  time_end = Sys.time() - time_start
  
  ### MCMC information
  mcmcInfo = data.frame(nbr.iter = nbr.iter,
                        nbr.chain = nbr.chain,
                        nbr.adapt = nbr.adapt,
                        nbr.thin=nbr.thin,
                        nbr.warmup=nbr.warmup,
                        total_time = time_end)
  
  ##
  ## OUTPUT
  ##
  
  OUT <- list(mcmc = mcmc,
              mcmc_Null = mcmc_Null,
              model = model,
              mcmcInfo = mcmcInfo,
              modelData = modelData_,
              model_type = model_type,
              estim.par = estim.par)

  class(OUT) <- c("survFitCstC","survFit")
  return(OUT)
}
