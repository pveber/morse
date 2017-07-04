#' @importFrom dplyr filter
survTKTDCreateJagsData <- function(data, comp) {
  # Creates the parameters to define the prior of the TKTD model
  # INPUTS
  # data : object of class survData
  # comp : if true return only min and max of prior
  # OUTPUT
  # jags.data : list of data required for the jags.model function
  
  data <- data[data$time != 0, ]
  
  # Parameter calculation of concentration min and max
  concmin <- min(data$conc[data$conc != 0])
  concmax <- max(data$conc)
  
  tmin <- min(data$time)
  tmax <- max(data$time)
  conc <- sort(unique(data$conc))
  
  deltaCmin = NULL
  for (i in 2:length(conc)) {
    deltaCmin[i - 1] <- conc[i] - conc[i - 1]
  }
  deltaCmin <- min(deltaCmin)
  
  # ks parameters
  ksmax <- -log(0.001) / (tmin * deltaCmin)
  ksmin <- -log(0.999) / (tmax * (concmax - concmin))
  
  meanlog10ks <- (log10(ksmax) + log10(ksmin)) / 2
  sdlog10ks <- (log10(ksmax) - log10(ksmin)) / 4
  taulog10ks <- 1 / sdlog10ks^2
  
  # kd parameters
  kdmax <- -log(0.001) / tmin
  kdmin <- -log(0.999) / tmax
  
  meanlog10kd <- (log10(kdmax) + log10(kdmin)) / 2
  
  sdlog10kd <- (log10(kdmax) - log10(kdmin)) / 4
  taulog10kd <- 1 / sdlog10kd^2
  
  # m0 parameters
  m0max <- -log(0.5) / tmin
  m0min <- -log(0.999) / tmax
  
  meanlog10m0 <- (log10(m0max) + log10(m0min)) / 2
  sdlog10m0 <- (log10(m0max) - log10(m0min)) / 4
  taulog10m0 <- 1/ sdlog10m0^2
  
  # nec parameters
  meanlog10nec <- (log10(concmax) + log10(concmin))/2
  sdlog10nec <- (log10(concmax) - log10(concmin)) / 4 
  taulog10nec <- 1/ sdlog10nec^2
  
  if (!comp) {
    return(list( x = data$conc, y = data$N_alive,
                 t = data$time, tprec = data$tprec,
                 Nprec = data$Nprec,
                 meanlog10ks = meanlog10ks, taulog10ks = taulog10ks,
                 meanlog10kd = meanlog10kd,
                 taulog10kd = taulog10kd,
                 meanlog10m0 = meanlog10m0,
                 taulog10m0 = taulog10m0,
                 meanlog10nec = meanlog10nec, taulog10nec = taulog10nec,
                 ndat = length(data$conc),
                 bigtime = max(data$time) + 10))
  } else {
    return(list(log10necmin = log10(concmin),
                log10necmax = log10(concmax),
                log10ksmin = log10(ksmin),
                log10ksmax = log10(ksmax),
                log10kdmin = log10(kdmin),
                log10kdmax = log10(kdmax),
                log10m0min = log10(m0min),
                log10m0max = log10(m0max)))
  }
}

modelTKTDNorm <- "model {
#########priors 
log10ks ~ dnorm(meanlog10ks, taulog10ks)
log10NEC ~ dnorm(meanlog10nec, taulog10nec)
log10kd ~ dnorm(meanlog10kd, taulog10kd)
log10m0 ~ dnorm(meanlog10m0, taulog10m0)

#####parameter transformation
ks <- 10**log10ks
NEC <- 10**log10NEC
kd <- 10**log10kd
m0 <- 10**log10m0

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1/kd * log( 1- R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/kd * ( exp(-kd * t[i]) - exp(-kd * tref[i]))), 0))
  
  y[i] ~ dbin(psurv[i] , Nprec[i]) 
}
}"

survTKTDPARAMS <- function(mcmc) {
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
  
  kd <- 10^res.M$quantiles["log10kd", "50%"]
  kdinf <- 10^res.M$quantiles["log10kd", "2.5%"]
  kdsup <- 10^res.M$quantiles["log10kd", "97.5%"]
  
  ks <- 10^res.M$quantiles["log10ks", "50%"]
  ksinf <- 10^res.M$quantiles["log10ks", "2.5%"]
  kssup <- 10^res.M$quantiles["log10ks", "97.5%"]
  nec <- 10^res.M$quantiles["log10NEC", "50%"]
  necinf <- 10^res.M$quantiles["log10NEC", "2.5%"]
  necsup <- 10^res.M$quantiles["log10NEC", "97.5%"]
  
  m0 <- 10^res.M$quantiles["log10m0", "50%"]
  m0inf <- 10^res.M$quantiles["log10m0", "2.5%"]
  m0sup <- 10^res.M$quantiles["log10m0", "97.5%"]
  
  # Definition of the parameter storage and storage data
  
  rownames <- c("kd", "ks", "nec", "m0")
  params <- c(kd, ks, nec, m0)
  CIinf <- c(kdinf, ksinf, necinf, m0inf)
  CIsup <- c(kdsup, kssup, necsup, m0sup)
  
  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)
  
  return(res)
}

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
#' \item{warnings}{a data.frame with warning messages}
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
  # test class object
  if(!is(data, "survData"))
    stop("survFitTKTD: object of class survData expected")
  
  # data transformation
  data <- summarise(group_by(data, conc, time), N_alive = sum(Nsurv))
  
  n <- nrow(data)
  data$tprec <- NA
  data$Nprec <- NA
  data$N_init <- NA
  for (i in 1:n)
  {
    if (data$time[i] != 0)
    {
      data$tprec[i] <- data$time[i - 1]
      data$Nprec[i] <- data$N_alive[i - 1]
      data$N_init[i] <- data$N_alive[data$conc == data$conc[i] & data$time == 0]
    }
  }
  
  # control
  datasurv0 <- subset(data, time == min(data$time[data$time != 0]))
  datasurv0$time <- 0
  datasurv0$N_alive <- datasurv0$N_init
  data[is.na(data$tprec),
       c("tprec", "Nprec", "N_init")] <- datasurv0[, c("tprec", "Nprec", "N_init")]
  
  jags.data <- survTKTDCreateJagsData(data, FALSE)
  
  # Define model
  
  model <- survLoadModel(model.program = modelTKTDNorm,
                         data = jags.data, n.chains,
                         Nadapt = 3000, quiet)
  
  # Determine sampling parameters
  parameters <- c("log10kd", "log10NEC","log10ks", "log10m0")
  
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
  
  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survTKTDPARAMS(mcmc)
  
  # check the posterior range
  priorBonds <- survTKTDCreateJagsData(data, TRUE)
  
  warnings <- warningTableCreate() 
  
  if (log10(estim.par["ks", "Q97.5"]) > priorBonds$log10ksmax){
    ##store warning in warnings table
    msg <- "The estimation of the killing rate (model parameter ks) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !"
    warnings <- warningTableAdd(warnings, "ks_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  
  if (log10(estim.par["kd", "Q97.5"]) > priorBonds$log10kdmax){
    ##store warning in warnings table
    msg <- "The estimation of the dominant rate constant (model parameter kd) lies outside the range used to define its prior distribution which indicates that this rate is very high and difficult to estimate from this experiment !"
    warnings <- warningTableAdd(warnings, "kd_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }

  
  if (log10(estim.par["m0", "Q2.5"]) < priorBonds$log10m0min){
    ##store warning in warnings table
    msg <- "The estimation of the natural instantaneous mortality rate (model parameter m0) lies outside the range used to define its prior distribution which indicates that this rate is very low and so difficult to estimate from this experiment !"
    warnings <- warningTableAdd(warnings, "hb_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
    
  }
  
  if (log10(estim.par["nec", "Q2.5"]) < priorBonds$log10necmin ||
      log10(estim.par["nec", "Q97.5"]) > priorBonds$log10necmax){
    ##store warning in warnings table
    msg <- "The NEC estimation (model parameter nec) lies outside the range of tested concentration and may be unreliable as the prior distribution on this parameter is defined from this range !"
    warnings <- warningTableAdd(warnings, "nec_outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }
  
  #OUTPUT
  OUT <- list(estim.par = estim.par,
              mcmc = mcmc,
              warnings = warnings,
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

