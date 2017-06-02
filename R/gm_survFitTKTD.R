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
#' @param data An object of class \code{gm_modelData}.
#' @param n.chains Number of MCMC chains. The minimum required number of chains
#' is 2.
#' @return The function returns an object of class \code{gm_survFitTKTD}, which is
#' a list with the following fields:
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior
#' distributions}
#' \item{model}{a JAGS or STAN model object}
#' \item{parameters}{a list of the parameters names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{n.iter}{a list of two indices indicating the beginning and end of
#' monitored iterations}
#' \item{n.thin}{a numerical value corresponding to the thinning interval}
#'
#' @keywords estimation
#
#' @examples
#'
#'
#' @export
#' @import rjags
#'


gm_survFitTKTD <- function(data,
                        model_type = NULL,
                        quiet = FALSE,
                        variable.names = NULL,
                        nbr.chain = 3,
                        nbr.adapt = 1000,
                        n.iter = NULL,
                        nbr.warmup = NULL,
                        nbr.thin = NULL){

  ##
  ## Pre modelling measure and tests
  ##

  ### time start
  time_start <- Sys.time()

  ### class of object
  if(!is(data, "gm_survData")){
    stop("gm_survFitTKTD: object of class 'gm_modelData' expected")
  }

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

  globalData = gm_modelData(data,
                            model_type = model_type)

  modelData = globalData$modelData
  modelData_Null = globalData$modelData_Null
  priorsData = globalData$priorsMinMax

  ##
  ## Define model
  ##

  cst_conc <- globalData$cst_conc

  if(model_type == "SD"){
    ### Determine sampling parameters
    parameters <- c("kd_log10", "hb_log10", "kk_log10", "z_log10")

    if(cst_conc){
      file_to_use <- jags_TKTD_cstSD
    } else{
      file_to_use <- jags_TKTD_varSD
    }
  } else if(model_type == "IT"){
    ### Determine sampling parameters
    parameters <- c("kd_log10", "hb_log10","alpha_log10", "beta_log10")

    if(cst_conc){
      file_to_use <- jags_TKTD_cstIT
    } else{
      file_to_use <- jags_TKTD_varIT
    }
  }

  model_Null <- survLoadModel(model.program = file_to_use,
                              data = modelData_Null,
                              n.chains = nbr.chain,
                              Nadapt = nbr.adapt,
                              quiet = TRUE)

  # model_Null <- jags.model(file = file_to_use,
  #                          data = modelData_Null,
  #                          inits,
  #                          n.chains = nbr.chain,
  #                          n.adapt = nbr.adapt,
  #                          quiet = TRUE)


  model <- survLoadModel(model.program = file_to_use,
                         data = modelData,
                         n.chains = nbr.chain,
                         Nadapt = nbr.adapt,
                         quiet = TRUE)

  # model <- jags.model(file = file_to_use,
  #                     data = modelData,
  #                     inits,
  #                     n.chains = nbr.chain,
  #                     n.adapt = nbr.adapt,
  #                     quiet = TRUE)


  ##
  ## estimate the number of iteration required for convergency of chains
  ## by using the raftery.diag
  ##

  if(is.null(nbr.warmup) | is.null(nbr.thin) | is.null(nbr.iter)){

    sampling.parameters <- modelSamplingParameters(model,
                                                   parameters,
                                                   n.chains)
    if (sampling.parameters$niter > 200000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    nbr.warmup = sampling.parameters$burnin
    nbr.thin = sampling.parameters$thin
    nbr.iter = sampling.parameters$niter

  }

  ### Null model to check priors with the model
  update(model_Null, nbr.warmup)

  mcmc_Null =  coda.samples(model_Null,
                            variable.names = parameters,
                            n.iter = nbr.iter,
                            thin = nbr.thin,
                            quiet = TRUE)

  ### model to check priors with the model
  update(model, nbr.warmup)
  mcmc =  coda.samples(model,
                       variable.names = parameters,
                       n.iter = nbr.iter,
                       thin = nbr.thin)

  ##
  ## Cheking posterior range with data from experimental design:
  ##

  estim.par <- gm_survTKTDPARAMS(mcmc, model_type = model_type)

  if (estim.par["kk", "Q97.5"] > priorsData$kk_max)
    warning("The estimation of the killing rate (model parameter k) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
            call. = FALSE)

  if (estim.par["kd", "Q97.5"] > priorsData$kd_max)
    warning("The estimation of the dominant rate constant (model parameter kd) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !",
            call. = FALSE)

  if (estim.par["hb", "Q2.5"] < priorsData$hb_min)
    warning("The estimation of the natural instantaneous mortality rate (model parameter hb) lies outside the range used to define its prior distribution which indicates that this rate is very low and so difficult to estimate from this experiment !",
            call. = FALSE)

  if (estim.par["z", "Q2.5"] < priorsData$z_min ||
      estim.par["z", "Q97.5"] > priorsData$z_max)
    warning("The estimation of Non Effect Concentration threshold (NEC) (model parameter z) lies outside the range of tested concentration and may be unreliable as the prior distribution on this parameter is defined from this range !",
            call. = FALSE)

  ### time end
  time_end = Sys.time() - time_start

  ##
  ## OUTPUT
  ##

  OUT <- list(mcmc = mcmc,
              mcmc_Null = mcmc_Null,
              data = gm_survData,
              estim.par = estim.par)

  class(OUT) <- "gm_survFitTKTD"
  return(OUT)
}


jags_TKTD_varIT <-"model {
  #------------------------------------------ parameter transformation
kd_taulog10 <- 1 / kd_sdlog10^2
hb_taulog10 <- 1 / hb_sdlog10^2
alpha_taulog10 <- 1 / alpha_sdlog10^2

#-------------------------------------------------------- priors
kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

alpha_log10 ~ dnorm(alpha_meanlog10, alpha_taulog10)
beta_log10 ~ dunif(beta_minlog10, beta_maxlog10)

#------------------------------------------ parameter transformation

kd <- 10**kd_log10
hb <- 10**hb_log10

alpha <- 10**alpha_log10
beta <- 10**beta_log10

##------------------------------------ model

for( i in 1:n_dataLong){

# methode des trapezes :
diff.int[profile_ID_long[i], time_ID_long[i]] <-  (exp(kd * time_long[i]) * conc_long[i] + exp(kd * tprec_long[i])
* concprec_long[i]) / 2 * (time_long[i] - tprec_long[i])

D_int[profile_ID_long[i], time_ID_long[i]] <- kd * exp(-kd * time_long[i]) *
sum( diff.int[profile_ID_long[i], 1:time_ID_long[i]] )

}


for( j in 1:n_dataRed){

D[profile_ID[j], time_ID[j]]  <- D_int[profile_ID[j], time_ID_red[j]]

D_max[profile_ID[j], time_ID[j]] <- max(D[profile_ID[j],1:time_ID[j]])

F[j]  <- D_max[profile_ID[j], time_ID[j]]**beta / ( D_max[profile_ID[j], time_ID[j]]**beta + alpha**beta )

psurv[j] <-  exp(-hb * time[j]) * (1- F[j])


Nsurv[i] ~ ifelse(time_ID[i] > 1,
dbin(psurv[i]/psurv[i-1] , Nprec[i]),
dbin(1 , Nprec[i]) )

# Nsurv[i] ~ ifelse(time_ID[i] > 1,
#      dbin(psurv[i]/psurv[i-1] , Nprec[i]),
#      dbin(1 , Nprec[i]) )
#
# ## ---------------------- generated data
#
# Nsurv_ppc[i] ~ ifelse(time_ID[i] > 1,
#     dbin(psurv[i]/psurv[i-1] , Nprec[i]),
#     dbin(1 , Nprec[i]) )
#
# Nsurv_sim[i] ~ ifelse(time_ID[i] > 1,
#     dbin(psurv[i]/psurv[i-1] , Nsurv_sim[i-1]),
#     dbin(1 , Nprec[i]) )

}
}"



gm_survTKTDPARAMS <- function(mcmc, model_type) {
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

  hb <- 10^res.M$quantiles["hb_log10m0", "50%"]
  hb_inf95 <- 10^res.M$quantiles["hb_log10m0", "2.5%"]
  hb_sup95 <- 10^res.M$quantiles["hb_log10m0", "97.5%"]

  if(model_type == "SD"){
    kk <- 10^res.M$quantiles["kk_log10", "50%"]
    kk_inf95 <- 10^res.M$quantiles["kk_log10", "2.5%"]
    kk_sup95 <- 10^res.M$quantiles["kk_log10", "97.5%"]

    z <- 10^res.M$quantiles["z_log10", "50%"]
    z_inf95 <- 10^res.M$quantiles["z_log10", "2.5%"]
    z_sup95 <- 10^res.M$quantiles["z_log10", "97.5%"]

    res <- data.frame(parameters = c("kd", "hb", "kk", "z"),
                      median = c(kd, hb, kk, z),
                      Q2.5 = c(kd_inf95, hb_inf95, kk_inf95, z_inf95),
                      Q97.5 = c(kd_sup95, hb_sup95, kk_sup95, z_sup95))

  } else if (model_type == "IT"){
    alpha <- 10^res.M$quantiles["alpha_log10", "50%"]
    alpha_inf95 <- 10^res.M$quantiles["alpha_log10", "2.5%"]
    alpha_sup95 <- 10^res.M$quantiles["alpha_log10", "97.5%"]

    beta <- 10^res.M$quantiles["beta_log10", "50%"]
    beta_inf95 <- 10^res.M$quantiles["beta_log10NEC", "2.5%"]
    beta_sup95 <- 10^res.M$quantiles["beta_log10NEC", "97.5%"]

    res <- data.frame(parameters = c("kd", "hb", "alpha", "beta"),
                      median = c(kd, hb, alpha, beta),
                      Q2.5 = c(kd_inf95, hb_inf95, alpha_inf95, beta_inf95),
                      Q97.5 = c(kd_sup95, hb_sup95, alpha_sup95, beta_sup95))
  } else {
    stop("please, provide the 'model_type': 'IT' or 'SD'")
  }

  return(res)
}
