#' @importFrom dplyr filter
survTKTDCreateJagsData <- function(data, distr) {
  # Creates the parameters to define the prior of the TKTD model
  # INPUTS
  # data: object of class survData
  # OUTPUT
  # jags.data : list of data required for the jags.model function
  
  data <- data[data$time != 0, ]
  
  # Parameter calculation of concentration min and max
  concmin <- min(data$conc[data$conc != 0])
  concmax <- max(data$conc)
  
  tmin <- min(filter(data, time != 0)$time)
  tmax <- max(data$time)
  
  deltaCmin = NULL
  for (i in 2:length(unique(data$conc))) {
    deltaCmin[i - 1] <- unique(data$conc)[i] - unique(data$conc)[i - 1]
  }
  deltaCmin <- min(deltaCmin)
  
  # ks parameters
  ksmax <- -log(0.01) / (tmin * deltaCmin)
  ksmin <- -log(0.99) / (tmax * (concmax - concmin))
  meanlog10ks <- (log10(ksmax) + log10(ksmin)) / 2
  sdlog10ks <- (log10(ksmax) - log10(ksmin)) / 4
  taulog10ks <- 1 / sdlog10ks^2
  
  # ke parameters
  kemax <- -log(0.01) / tmin
  kemin <- -log(0.99) / tmax
  meanlog10ke <- (log10(kemax) + log10(kemin)) / 2
  
  sdlog10ke <- (log10(kemax) - log10(kemin)) / 4
  taulog10ke <- 1 / sdlog10ke^2
  
  # m0 parameters
  m0max <- -log(0.5) / tmin
  m0min <- -log(0.99) / tmax
  meanlog10m0 <- (log10(m0max) + log10(m0min)) / 2
  
  sdlog10m0 <- (log10(m0max) - log10(m0min)) / 4
  taulog10m0 <- 1/ sdlog10m0^2
  
  # nec parameters
  meanlog10nec <- (log10(concmax) + log10(concmin))/2
  sdlog10nec <- (log10(concmax) - log10(concmin)) / 4 
  taulog10nec <- 1/ sdlog10nec^2
  
  if (distr == "unif") {
    return(list( 'x' = data$conc, 'y' = data$N_alive,
                 't' = data$time, 'tprec' = data$tprec,
                 'Nprec' = data$Nprec,
                 concmin = concmin, concmax = concmax,
                 ksmin = ksmin, ksmax = ksmax,
                 kemin = kemin, kemax = kemax,
                 m0min = m0min, m0max = m0max,
                 ndat = length(data$conc),
                 bigtime = max(data$time) + 10))
  } else if (distr == "norm") {
    return(list( 'x' = data$conc, 'y' = data$N_alive,
                 't' = data$time, 'tprec' = data$tprec,
                 'Nprec' = data$Nprec,
                 meanlog10ks = meanlog10ks, taulog10ks = taulog10ks,
                 meanlog10ke = meanlog10ke, taulog10ke = taulog10ke,
                 meanlog10m0 = meanlog10m0, taulog10m0 = taulog10m0,
                 meanlog10nec = meanlog10nec, taulog10nec = taulog10nec,
                 ndat = length(data$conc),
                 bigtime = max(data$time) + 10))
  }
}

modelTKTDUnif <- "model {
#########priors 
lks ~ dunif(log(ksmin) / log(10), log(ksmax) / log(10))
lNEC ~ dunif(log(concmin) / log(10) - 1, log(concmax) / log(10) + 1)
lke ~ dunif(log(kemin) / log(10), log(kemax) / log(10))
lm0 ~ dunif(log(m0min) / log(10), log(m0max) / log(10))

#####parameter transformation
ks <- 10**lks
NEC <- 10**lNEC
ke <- 10**lke
m0 <- 10**lm0

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1 / ke * log(1 - R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), 0))
  
  y[i] ~ dbin(psurv[i] , Nprec[i]) 
}
}"

modelTKTDNorm <- "model {
#########priors 
lks ~ dnorm(meanlog10ks, taulog10ks)
lNEC ~ dnorm(meanlog10nec, taulog10nec)
lke ~ dnorm(meanlog10ke, taulog10ke)
lm0 ~ dnorm(meanlog10m0, taulog10m0)

#####parameter transformation
ks <- 10**lks
NEC <- 10**lNEC
ke <- 10**lke
m0 <- 10**lm0

##########Computation of the likelihood
for (i in 1:ndat)
{
  tNEC[i] <- ifelse(x[i] > NEC, -1/ke * log( 1- R[i]), bigtime)
  R[i] <- ifelse(x[i] > NEC, NEC/xcor[i], 0.1)
  xcor[i] <- ifelse(x[i] > 0, x[i], 10)
  tref[i] <- max(tprec[i], tNEC[i])
  
  psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tNEC[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/ke * ( exp(-ke * t[i]) - exp(-ke * tref[i]))), 0))
  
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
  
  ke <- 10^res.M$quantiles["lke", "50%"]
  keinf <- 10^res.M$quantiles["lke", "2.5%"]
  kesup <- 10^res.M$quantiles["lke", "97.5%"]
  ks <- 10^res.M$quantiles["lks", "50%"]
  ksinf <- 10^res.M$quantiles["lks", "2.5%"]
  kssup <- 10^res.M$quantiles["lks", "97.5%"]
  nec <- 10^res.M$quantiles["lNEC", "50%"]
  necinf <- 10^res.M$quantiles["lNEC", "2.5%"]
  necsup <- 10^res.M$quantiles["lNEC", "97.5%"]
  m0 <- 10^res.M$quantiles["lm0", "50%"]
  m0inf <- 10^res.M$quantiles["lm0", "2.5%"]
  m0sup <- 10^res.M$quantiles["lm0", "97.5%"]
  
  # Definition of the parameter storage and storage data
  # If Poisson Model
  
  rownames <- c("ke", "ks", "nec", "m0")
  params <- c(ke, ks, nec, m0)
  CIinf <- c(keinf, ksinf, necinf, m0inf)
  CIsup <- c(kesup, kssup, necsup, m0sup)
  
  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)
  
  return(res)
}

#' @import coda
coda.samples.dic <- function (model, variable.names = NULL, n.iter, thin = 1, ...)
{
  load.module('dic') # necessary for pD and deviance monitor
  
  start <- model$iter() + thin
  varnames = c(variable.names, c('deviance', 'pD'))
  out <- jags.samples(model, varnames, n.iter, thin,
                      type = "trace", ...)
  deviance <- out$deviance
  pD <- out$pD
  out$deviance <- NULL
  out$pD <- NULL
  ans <- vector("list", rjags:::nchain(model))
  for (ch in 1:rjags:::nchain(model)) {
    ans.ch <- vector("list", length(out))
    vnames.ch <- NULL
    for (i in seq(along = out)) {
      varname <- names(out)[[i]]
      d <- dim(out[[i]])
      if (length(d) < 3) {
        stop("Invalid dimensions for sampled output")
      }
      vardim <- d[1:(length(d) - 2)]
      nvar <- prod(vardim)
      niter <- d[length(d) - 1]
      nchain <- d[length(d)]
      values <- as.vector(out[[i]])
      var.i <- matrix(NA, nrow = niter, ncol = nvar)
      for (j in 1:nvar) {
        var.i[, j] <- values[j + (0:(niter - 1)) * nvar +
                               (ch - 1) * niter * nvar]
      }
      vnames.ch <- c(vnames.ch, rjags:::coda.names(varname, vardim))
      ans.ch[[i]] <- var.i
    }
    ans.ch <- do.call("cbind", ans.ch)
    colnames(ans.ch) <- vnames.ch
    ans[[ch]] <- mcmc(ans.ch, start = start, thin = thin)
  }
  
  dic <- list(deviance = mean(as.vector(deviance)), penalty = mean(as.vector(pD)), type = 'pD')
  class(dic) <- "dic"
  return(list(samples = mcmc.list(ans),
              dic = round(sum(sapply(dic$deviance, mean) + sapply(dic$penalty, mean)))))
}

#' Fit a Bayesian TKTD model for survival analysis among time
#' 
#' The \code{survFitTKTD} function estimates the parameters of a TKTD model
#' for survival analysis using Bayesian inference.
#' 
#' The function returns
#' parameter estimates of the TKTD model.
#' 
#' @param data An object of class \code{survData}.
#' @param distr If \code{"norm"} gives normal distributions for prior distribution
#' parameters, else if \code{"unif"} gives uniform distributions for prior
#' distribution parameters.
#' @param n.chains Number of MCMC chains. The minimum required number of chains
#' is 2.
#' @param quiet If \code{TRUE}, make silent all prints and progress bars of
#' JAGS compilation.
#' 
#' @return The function returns an object of class \code{survFitTKTD}. A list
#' of 10 objects:
#' \item{DIC}{DIC value of the selected model.}
#' \item{estim.par}{A table of the estimated parameters as medians and 95 \%
#' credible intervals.}
#' \item{mcmc}{An object of class \code{mcmc.list} with the posterior
#' distributions.}
#' \item{model}{A JAGS model object.}
#' \item{parameters}{A list of the parameters names used in the model.}
#' \item{n.chains}{An integer value corresponding to the number of chains used
#' for the MCMC computation.}
#' \item{n.iter}{A list of two numerical value corresponding to the beginning
#' and the end of monitored iterations.}
#' \item{n.thin}{A numerical value corresponding to the thinning interval.}
#' \item{jags.data}{A list of data used by the internal \code{\link[rjags]{jags.model}}
#' function. This object is intended for the case when the user wishes to use
#' the \code{\link[rjags]{rjags}} package instead of the automatied estimation
#' function.}
#' \item{transformed.data}{The \code{survData} object.
#' See \code{\link{survData}} for details.}
#' 
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @references Plummer, M. (2013) JAGS Version 3.4.0 user manual.
#' \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}
#'
#' Spiegelhalter, D., N. Best, B. Carlin, and A. van der Linde (2002) Bayesian
#' measures of model complexity and fit (with discussion).  \emph{Journal of
#' the Royal Statistical Society}, Series B 64, 583-639.
#'
#' @keywords estimation
#
#' @examples
#' 
#' # (1) Load the survival data
#' data(cadmium1)
#' 
#' # (2) Create an object of class "survData"
#' dat <- survData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#' }
#' 
#' @export
#' @import rjags
#' @importFrom dplyr group_by summarise
#' 
survFitTKTD <- function(data,
                        distr = "norm",
                        n.chains = 3,
                        quiet = FALSE) {
  
  # data transformation
  # agregate by sum of replicate
  data <- cbind(aggregate(Nsurv ~ time + conc, data, sum),
                replicate = 1)
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
  
  jags.data <- survTKTDCreateJagsData(data, distr)
  
  # Define model
  if (distr == "norm") {
    model <- survLoadModel(model.program = modelTKTDNorm,
                           data = jags.data, n.chains,
                           Nadapt = 3000, quiet)
  } else if (distr == "unif") {
    model <- survLoadModel(model.program = modelTKTDUnif,
                           data = jags.data, n.chains,
                           Nadapt = 3000, quiet)
  }
  
  # Determine sampling parameters
  parameters <- c("lke", "lNEC","lks", "lm0")
  sampling.parameters <- modelSamplingParameters(model,
                                                 parameters, n.chains, quiet)
  
#   if (sampling.parameters$niter > 100000)
#     stop("The model needs too many iterations to provide reliable parameter estimates !")
  
  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  
  mcmc <- coda.samples.dic(model, parameters,
                           n.iter = sampling.parameters$niter,
                           thin = sampling.parameters$thin,
                           progress.bar = prog.b)
  
  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survTKTDPARAMS(mcmc$samples)
  
  #OUTPUT
  OUT <- list(DIC = mcmc$dic,
              estim.par = estim.par,
              mcmc = mcmc$samples,
              model = model,
              parameters = parameters,
              n.chains = summary(mcmc$samples)$nchain,
              n.iter = list(start = summary(mcmc$samples)$start,
                            end = summary(mcmc$samples)$end),
              n.thin = summary(mcmc$samples)$thin,
              jags.data = jags.data,
              transformed.data = data)
  
  class(OUT) <- "survFitTKTD"
  return(OUT)
}

