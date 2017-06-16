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
#' @return The function returns an object of class \code{survFitTKTD}, which is
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
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
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
survFitTKTD <- function(data,
                        n.chains = 3,
                        quiet = FALSE) {
  
  
  warning("'survFitTKTD' is deprecated with morse version >= 3.0.0, please use the function 'survFit'")
  
  model_type = "SD"
  
  # test class object
  if(!is(data, "survData"))
    stop("survFitTKTD: object of class survData expected")
  
  
  ##
  ## Data and Priors for model
  ##
  
  globalData <- modelData(data,
                         model_type = model_type)
  
  modelData_ <- globalData$modelData
  jags.data  <- modelData_ ; jags.data$replicate = NULL

  priorsData = globalData$priorsMinMax
  
  # Define model
  
  model <- survLoadModel(model.program = jags_TKTD_cstSD,
                         data = jags.data,
                         n.chains,
                         Nadapt = 3000,
                         quiet)
  
  # Determine sampling parameters
  parameters <- c("kd_log10", "z_log10","kk_log10", "hb_log10")
  
  sampling.parameters <- modelSamplingParameters(model,
                                                 parameters, n.chains, quiet)
  
  if (sampling.parameters$niter > 200000)
    stop("The model needs too many iterations to provide reliable parameter estimates !")
  
  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  
  mcmc <- coda.samples(model, parameters,
                       n.iter = sampling.parameters$niter,
                       thin = sampling.parameters$thin,
                       progress.bar = prog.b)
  
  # check the posterior range
  priorsData <- priors_survData(data, model_type)$priorsMinMax
  
  ##
  ## Cheking posterior range with data from experimental design:
  ##
  
  estim.par <- survTKTDPARAMS(mcmc, model_type)
  
  if (filter(estim.par, parameters == "kd")$Q97.5 > priorsData$kd_max){
    warning("The estimation of the dominant rate constant (model parameter kd) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
            call. = FALSE)
  }
  if (filter(estim.par, parameters == "hb")$Q2.5 < priorsData$hb_min){
    warning("The estimation of the natural instantaneous mortality rate (model parameter hb) lies outside the range used to define its prior distribution which indicates that this rate is very low and so difficult to estimate from this experiment !",
            call. = FALSE)
  }
  
  if (filter(estim.par, parameters == "kk")$Q97.5 > priorsData$kk_max)
    warning("The estimation of the killing rate (model parameter k) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
            call. = FALSE)
  
  if (filter(estim.par, parameters == "z")$Q2.5 < priorsData$conc_min ||
      filter(estim.par, parameters == "z")$Q97.5 > priorsData$conc_max)
    warning("The estimation of Non Effect Concentration threshold (NEC) (model parameter z) lies outside the range of tested concentration and may be unreliable as the prior distribution on this parameter is defined from this range !",
            call. = FALSE)

  #OUTPUT
  OUT <- list(estim.par = estim.par,
              mcmc = mcmc,
              model = model,
              parameters = parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = data)
  
  class(OUT) <- "survFitTKTD"
  return(OUT)
}

