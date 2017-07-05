#' Fits a Bayesian exposure-response model for target-time reproduction analysis
#'
#' This function estimates the parameters of an exposure-response model for
#' target-time reproduction analysis using Bayesian inference.
#' In this model the response is the cumulated reproduction outputs of a
#' population in a given time period in presence of mortality.
#'
#' Because some individuals may die during the observation period, the
#' reproduction rate alone is not sufficient to account for the observed number
#' of offspring. In addition, we need the time individuals have stayed alive
#' during the experiment. The \code{reproFitTT} function estimates the number
#' of individual-days in an experiment between its start and the target time.
#' This covariable is then used to estimate a relation between the toxicant
#' concentration and the reproduction rate \emph{per individual-day}.
#'
#' The \code{reproFitTT} function fits two models, one where inter-individual
#' variability is neglected ("Poisson" model) and one where it is taken into
#' account ("gamma-Poisson" model). When setting \code{stoc.part} to
#' \code{"bestfit"}, a model comparison procedure is used to choose between
#' them. More details are presented in the vignette accompanying the package.
#'
#' @param data an object of class \code{reproData}
#' @param stoc.part stochastic part of the model. Possible values are \code{"bestfit"},
#' \code{"poisson"} and \code{"gammapoisson"}
#' @param target.time defines the observation period. By default the last time point
#' @param ecx desired values of \eqn{x} (in percent) for which to compute
#' \eqn{EC_{x}}{ECx}
#' @param n.chains number of MCMC chains. The minimum required number of chains is 2
#' @param quiet if \code{TRUE}, does not print messages and progress bars from JAGS
#'
#'
#' @return The function returns an object of class \code{reproFitTT} which is a list
#' of the following objects:
#' \item{DIC}{DIC value of the selected model}
#' \item{estim.ECx}{a table of the estimated 5, 10, 20 and 50 \% effective
#' concentrations (by default) and their 95 \% credible intervals}
#' \item{estim.par}{a table of the estimated parameters as medians and 95 \%
#' credible intervals}
#' \item{mcmc}{an object of class \code{mcmc.list} with the posterior distributions}
#' \item{model}{a JAGS model object}
#' \item{warnings}{a data.frame with warning messages}
#' \item{model.label}{a character string, \code{"P"} if the poisson model is used,
#' \code{"GP"} if the gamma-poisson is used}
#' \item{parameters}{a list of the parameters names used in the model}
#' \item{n.chains}{an integer value corresponding to the number of chains used
#' for the MCMC computation}
#' \item{n.iter}{a list of two indices indicating the beginning and
#' the end of monitored iterations}
#' \item{n.thin}{a numerical value corresponding to the thinning interval}
#' \item{jags.data}{a list a the data passed to the jags model}
#' \item{transformed.data}{the \code{survData} object passed to the function}
#' \item{dataTT}{the dataset with which one the parameters are estimated}
#'
#' @keywords estimation
#'
#' @examples
#'
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create an object of class "reproData"
#' dat <- reproData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the reproFitTT function with the log-logistic gamma-poisson model
#' out <- reproFitTT(dat, stoc.part = "gammapoisson",
#'                   ecx = c(5, 10, 15, 20, 30, 50, 80), quiet = TRUE)
#' }
#'
#' @import rjags
#'
#' @export
reproFitTT <- function(data,
                       stoc.part = "bestfit",
                       target.time = NULL,
                       ecx = c(5, 10, 20, 50),
                       n.chains = 3,
                       quiet = FALSE) {
  # test class object
  if (! is(data, "reproData"))
    stop("reproFitTT: object of class reproData expected")

  # stocastic verification
  stoc.partpossible <- c("poisson", "gammapoisson", "bestfit")

  if (!any(stoc.partpossible == stoc.part))
    stop("Invalid value for argument [stoc.part]")

  # check 0 Nreprocumul
  if (all(data$Nreprocumul == 0))
    stop("Nreprocumul contains only 0 values !")

  # parameters
  parameters <- list(poisson = c("d", "log10b", "log10e"),
                     gammapoisson = c("d", "log10b","log10e", "log10omega"))

  # select Data at target.time
  dataTT <- selectDataTT(data, target.time)

  # create priors parameters
  jags.data <- reproCreateJagsData(stoc.part, dataTT)

  # Poisson model only
  if (stoc.part == "poisson") {
    # Define model
    poisson.model <- reproLoadPoissonModel(model.program = llm.poisson.model.text,
                                           data = jags.data,
                                           n.chains, quiet)

    # Determine sampling parameters
    poisson.sampling.parameters <- modelSamplingParameters(poisson.model,
                                                           parameters$poisson,
                                                           n.chains, quiet)

    if (poisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    poisson.DIC <- calcDIC(poisson.model, poisson.sampling.parameters, quiet)

    # list of objet for the coda.sample function
    coda.arg <- list(model = poisson.model,
                     model.label = "P",
                     niter = poisson.sampling.parameters$niter,
                     thin = poisson.sampling.parameters$thin,
                     nburnin = poisson.sampling.parameters$burnin,
                     parameters = parameters$poisson,
                     DIC = poisson.DIC)
  }

  # Gamma-poisson model only
  if (stoc.part == "gammapoisson") {
    # Define model
    gammapoisson.model <- reproLoadGammapoissonModel(model.program = llm.gammapoisson.model.text,
                                                     data = jags.data,
                                                     n.chains, quiet)

    # Determine sampling parameters
    gammapoisson.sampling.parameters <- modelSamplingParameters(gammapoisson.model,
                                                                parameters$gammapoisson,
                                                                n.chains, quiet)

    if (gammapoisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    gammapoisson.DIC <- calcDIC(gammapoisson.model,
                                gammapoisson.sampling.parameters, quiet)

    # list of objet for the coda.sample function
    coda.arg <- list(model = gammapoisson.model,
                     model.label = "GP",
                     niter = gammapoisson.sampling.parameters$niter,
                     thin = gammapoisson.sampling.parameters$thin,
                     nburnin = gammapoisson.sampling.parameters$burnin,
                     parameters = parameters$gammapoisson,
                     DIC = gammapoisson.DIC)
  }

  # Model Selection by the DIC
  if (stoc.part == "bestfit") {
    # Define models
    poisson.model <- reproLoadPoissonModel(model.program = llm.poisson.model.text,
                                           data = jags.data,
                                           n.chains, quiet)

    gammapoisson.model <- reproLoadGammapoissonModel(model.program = llm.gammapoisson.model.text,
                                                     data = jags.data,
                                                     n.chains, quiet)
    # Determine sampling parameters
    poisson.sampling.parameters <- modelSamplingParameters(poisson.model,
                                                           parameters$poisson,
                                                           n.chains, quiet)

    gammapoisson.sampling.parameters <- modelSamplingParameters(gammapoisson.model,
                                                                parameters$gammapoisson,
                                                                n.chains, quiet)

    if (poisson.sampling.parameters$niter > 100000 && gammapoisson.sampling.parameters$niter > 100000)
      stop("The model needs too many iterations to provide reliable parameter estimates !")

    # calcul DIC
    poisson.DIC <- calcDIC(poisson.model, poisson.sampling.parameters, quiet)
    gammapoisson.DIC <- calcDIC(gammapoisson.model,
                                gammapoisson.sampling.parameters, quiet)

    if (gammapoisson.sampling.parameters$niter > 100000) {
      # list of object for the coda.sample function
      coda.arg <- list(model = poisson.model,
                       model.label = "P",
                       niter = poisson.sampling.parameters$niter,
                       thin = poisson.sampling.parameters$thin,
                       nburnin = poisson.sampling.parameters$burnin,
                       parameters = parameters$poisson,
                       DIC = poisson.DIC)
    }

    if (poisson.sampling.parameters$niter > 100000) {
      # list of object for the coda.sample function
      coda.arg <- list(model = gammapoisson.model,
                       model.label = "GP",
                       niter = gammapoisson.sampling.parameters$niter,
                       thin = gammapoisson.sampling.parameters$thin,
                       nburnin = gammapoisson.sampling.parameters$burnin,
                       parameters = parameters$gammapoisson,
                       DIC = gammapoisson.DIC)
    }
    if (poisson.sampling.parameters$niter <= 100000 && gammapoisson.sampling.parameters$niter <= 100000) {
      if (poisson.DIC <= (gammapoisson.DIC + 10)) {
        # list of objet for the coda.sample function
        coda.arg <- list(model = poisson.model,
                         model.label = "P",
                         niter = poisson.sampling.parameters$niter,
                         thin = poisson.sampling.parameters$thin,
                         nburnin = poisson.sampling.parameters$burnin,
                         parameters = parameters$poisson,
                         DIC = poisson.DIC)
      } else {
        # list of objet for the coda.sample function
        coda.arg <- list(model = gammapoisson.model,
                         model.label = "GP",
                         niter = gammapoisson.sampling.parameters$niter,
                         thin = gammapoisson.sampling.parameters$thin,
                         nburnin = gammapoisson.sampling.parameters$burnin,
                         parameters = parameters$gammapoisson,
                         DIC = gammapoisson.DIC)
      }
    }
  }

  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  mcmc <- coda.samples(coda.arg$model,
                       coda.arg$parameters,
                       n.iter = coda.arg$niter,
                       thin = coda.arg$thin,
                       progress.bar = prog.b)

  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- reproPARAMS(mcmc, coda.arg$model.label)

  # ECx calculation  estimated ECx and their CIs 95%
  # vector of ECX
  estim.ECx <- estimXCX(mcmc, ecx, "EC")

  # check if the maximum measured concentration is in the EC50's range of
  # 95% percentile

  warnings <- msgTableCreate()

  EC50 <- log10(estim.par["e", "median"])
  if (!(min(log10(data$conc)) < EC50 & EC50 < max(log10(data$conc)))){
    ##store warning in warnings table
    msg <- "The EC50 estimation (model parameter e) lies outside the range of
    tested concentration and may be unreliable as the prior distribution on
    this parameter is defined from this range !"
    warnings <- msgTableAdd(warnings, "EC50outRange", msg)
    ## print the message
    warning(msg, call. = FALSE)
  }


  # output
  OUT <- list(DIC = coda.arg$DIC,
              estim.ECx = estim.ECx,
              estim.par = estim.par,
              det.part = "loglogistic",
              mcmc = mcmc,
              warnings = warnings,
              model = coda.arg$model,
              model.label = coda.arg$model.label,
              parameters = coda.arg$parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = data,
              dataTT = dataTT)

  class(OUT) <- "reproFitTT"
  return(OUT)
}


#' @importFrom stats sd
reproCreateJagsData <- function(stoc.part, data) {
  # create the parameters to define the prior of the log-logistic model
  # for reproduction data analysis
  # INPUTS
  # stoc.part: model name
  # data: object of class reproData
  # OUTPUT
  # jags.data : list data require for the jags.model function


  # separate control data to the other
  # tab0: data at conc = 0
  tab0 <- data[data$conc == min(data$conc), ]
  # tab: data at conc != 0
  tab <- data[data$conc != min(data$conc), ]

  Nindtime <- tab$Nindtime
  NreprocumulIndtime0 <- tab0$Nreprocumul / tab0$Nindtime # cumulated number of
  # offspring / number of
  # individual-days
  conc <- tab$conc
  Ncumul <- tab$Nreprocumul
  n <- nrow(tab) # number of observation != from the control

  # Parameter calculation of concentration min and max
  concmin <- min(sort(unique(conc))[-1])
  concmax <- max(conc)

  # create priors parameters for the log logistic model

  # Params to define log10e
  meanlog10e <- (log10(concmin) + log10(concmax)) / 2
  sdlog10e <- (log10(concmax) - log10(concmin)) / 4
  taulog10e <- 1 / sdlog10e^2

  # Params to define d
  meand <- mean(NreprocumulIndtime0)
  SEd <- sd(NreprocumulIndtime0) / sqrt(length(unique(tab0$replicate)))
  taud <- 1 / (SEd)^2

  # Params to define b
  log10bmin <- -2
  log10bmax <- 2

  # list of data use by jags
  jags.data <- list(meanlog10e = meanlog10e,
                    taulog10e = taulog10e,
                    meand = meand,
                    taud = taud,
                    log10bmin = log10bmin,
                    log10bmax = log10bmax,
                    n = n,
                    xconc = conc,
                    Nindtime = Nindtime,
                    Ncumul = Ncumul)

  # Params to define overdispersion rate
  if (stoc.part == "bestfit" || stoc.part == "gammapoisson") {
    log10omegamin <- -4
    log10omegamax <- 4

    # list of data use by jags
    jags.data <- c(jags.data,
                   log10omegamin = log10omegamin,
                   log10omegamax = log10omegamax)
  }
  return(jags.data)
}

reproLoadPoissonModel <- function(model.program,
                                  data,
                                  n.chains,
                                  quiet = quiet) {
                                    # sub function to load jags poisson model
                                    reproLoadModel(model.program, F, data, n.chains, quiet = quiet)
                                  }

reproLoadGammapoissonModel <- function(model.program,
                                       data,
                                       n.chains,
                                       quiet = quiet) {
                                         # sub function to load jags gamma poisson model
                                         reproLoadModel(model.program, T, data, n.chains, quiet = quiet)
}

#' @import rjags
reproLoadModel <- function(model.program,
                           lr.bound.keep,
                           data,
                           n.chains,
                           Nadapt = 3000,
                           quiet = quiet) {
  # create the JAGS model object and called by reproLoadPoissonModel
  # and reproLoadGammapoissonModel
  # INPUTS:
  # - model.program: character string containing a jags model description
  # - lr.bound.keep: boolean value to use omega parameter or not
  # - data: list of data created by reproCreateJagsData
  # - nchains: Number of chains desired
  # - Nadapt: length of the adaptation phase
  # - quiet: silent option
  # OUTPUT:
  # - JAGS model

  # delisting of lr.bound because not used in the function
  if (!lr.bound.keep) {
    data[c("meanlog10omega", "taulog10omega", "log10omegamin",
           "log10omegamax")] <- NULL
    }

    # load model text in a temporary file
    model.file <- tempfile() # temporary file address
    fileC <- file(model.file) # open connection
    writeLines(model.program, fileC) # write text in temporary file
    close(fileC) # close connection to temporary file
    # creation of the jags model
    model <- jags.model(file = model.file, data = data, n.chains = n.chains,
                        n.adapt = Nadapt, quiet = quiet)
    unlink(model.file)
    return(model)
}

reproPARAMS <- function(mcmc, MODEL = "P") {
  # create the table of posterior estimated parameters
  # for the reproduction analyses
  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chain
  # - MODEL: a position flag model with P: poisson model and GP: gammapoisson
  # model
  # OUTPUT:
  # - data frame with 3 columns (values, CIinf, CIsup) and 3-4rows (the estimated
  # parameters)

  # Retrieving parameters of the model
  res.M <- summary(mcmc)

  b <- 10^res.M$quantiles["log10b", "50%"]
  d <- res.M$quantiles["d", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  binf <- 10^res.M$quantiles["log10b", "2.5%"]
  dinf <- res.M$quantiles["d", "2.5%"]
  einf <- 10^res.M$quantiles["log10e", "2.5%"]
  bsup <- 10^res.M$quantiles["log10b", "97.5%"]
  dsup <- res.M$quantiles["d", "97.5%"]
  esup <- 10^res.M$quantiles["log10e", "97.5%"]

  # Definition of the parameter storage and storage data

  # If Poisson Model
  if (MODEL == "P") {
    rownames <- c("b", "d", "e")
    params <- c(b, d, e)
    CIinf <- c(binf, dinf, einf)
    CIsup <- c(bsup, dsup, esup)
  }
  # If Gamma Poisson Model
  if (MODEL == "GP") {
    # Calculation of the parameter omega
    omega <- 10^res.M$quantiles["log10omega", "50%"]
    omegainf <- 10^res.M$quantiles["log10omega", "2.5%"]
    omegasup <- 10^res.M$quantiles["log10omega", "97.5%"]
    # Definition of the parameter storage and storage data
    rownames <- c("b", "d", "e", "omega")
    params <- c(b, d, e, omega)
    CIinf <- c(binf, dinf, einf, omegainf)
    CIsup <- c(bsup, dsup, esup, omegasup)
  }

  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)

  return(res)
}

llm.poisson.model.text <- "\nmodel # Loglogistic Poisson model\n{\n#\nfor (j in 1:n) # loop on replicates\n{\n# Explicit writting of a Poisson law for each replicate\n# mean is given by the theoretical curve\nytheo[j] <- d / (1 + pow(xconc[j]/e, b))\nnbtheo[j] <- ytheo[j]*Nindtime[j]\nNcumul[j] ~ dpois(nbtheo[j])\n}\n# Prior distributions\nd ~ dnorm(meand, taud)T(0,)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10,log10b)\ne <- pow(10,log10e)\n}\n"

llm.gammapoisson.model.text <- "\nmodel # Loglogisitc Gamma poisson model\n{\n#\nfor (j in 1:n) # loop on replicates\n{\n# Explicit writting of a gamma-Poisson law for each replicate\n# the mean is given by a gamma law centered on the theoretical curve\nrate[j] <- d / (1 + pow(xconc[j]/e, b)) / omega\np[j] <- 1 / (Nindtime[j] * omega + 1)\nNcumul[j] ~ dnegbin(p[j], rate[j])\n}\n# Prior distributions\nd ~ dnorm(meand, taud)T(0,)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\nlog10omega ~ dunif(log10omegamin, log10omegamax)\n\nomega <- pow(10,log10omega)\nb <- pow(10,log10b)\ne <- pow(10,log10e)\n}\n"
