survCreateJagsData <- function(det.part, data) {
  # Creates the parameters to define the prior of the log-logistic binomial model
  # INPUTS
  # det.part: model name
  # data: object of class survData
  # OUTPUT
  # jags.data : list of data required for the jags.model function

  # Parameter calculation of concentration min and max
  concmin <- min(sort(unique(data$conc))[-1])
  concmax <- max(data$conc)

  # Create prior parameters for the log logistic model

  # Params to define e
  meanlog10e <- (log10(concmin) + log10(concmax)) / 2
  sdlog10e <- (log10(concmax) - log10(concmin)) / 4
  taulog10e <- 1 / sdlog10e^2

  # Params to define b
  log10bmin <- -2
  log10bmax <- 2

  # list of data use by jags
  jags.data <- list(meanlog10e = meanlog10e,
                    Ninit = data$Ninit,
                    Nsurv = data$Nsurv,
                    taulog10e = taulog10e,
                    log10bmin = log10bmin,
                    log10bmax = log10bmax,
                    n = length(data$conc),
                    xconc = data$conc)

  # list of data use by jags
  if (det.part == "loglogisticbinom_3") {
    jags.data <- c(jags.data,
                   dmin = 0,
                   dmax = 1)
  }
  return(jags.data)
}

#' @import rjags
survLoadModel <- function(model.program,
                          data,
                          n.chains,
                          Nadapt,
                          quiet = quiet) {
  # create the JAGS model object
  # INPUTS:
  # - model.program: character string containing a jags model description
  # - data: list of data created by survCreateJagsData
  # - nchains: Number of chains desired
  # - Nadapt: length of the adaptation phase
  # - quiet: silent option
  # OUTPUT:
  # - JAGS model

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

survPARAMS <- function(mcmc, det.part) {
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

  if (det.part ==  "loglogisticbinom_3") {
    d <- res.M$quantiles["d", "50%"]
    dinf <- res.M$quantiles["d", "2.5%"]
    dsup <- res.M$quantiles["d", "97.5%"]
  }
  # for loglogisticbinom_2 and 3
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  binf <- 10^res.M$quantiles["log10b", "2.5%"]
  einf <- 10^res.M$quantiles["log10e", "2.5%"]
  bsup <- 10^res.M$quantiles["log10b", "97.5%"]
  esup <- 10^res.M$quantiles["log10e", "97.5%"]

  # Definition of the parameter storage and storage data
  # If Poisson Model

  if (det.part == "loglogisticbinom_3") {
    # if mortality in control
    rownames <- c("b", "d", "e")
    params <- c(b, d, e)
    CIinf <- c(binf, dinf, einf)
    CIsup <- c(bsup, dsup, esup)
  } else {
    # if no mortality in control
    # Definition of the parameter storage and storage data
    rownames <- c("b", "e")
    params <- c(b, e)
    CIinf <- c(binf, einf)
    CIsup <- c(bsup, esup)
  }

  res <- data.frame(median = params, Q2.5 = CIinf, Q97.5 = CIsup,
                    row.names = rownames)

  return(res)
}

llbinom3.model.text <- "\nmodel # Loglogistic binomial model with 3 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- d/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nd ~ dunif(dmin, dmax)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"

llbinom2.model.text <- "\nmodel # Loglogistic binomial model with 2 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- 1/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"

#' Print of \code{survFitTT} object
#' 
#' This is the generic \code{print} S3 method for the \code{survFitTT} class.
#' It prints the underlying JAGS model and some information on the Bayesian 
#' inference procedure.
#' 
#' @param x An object of class \code{survFitTT}
#' 
#' @seealso \code{\link{survFitTT}}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the survFitTT function with the log-logistic
#' # binomial model
#' out <- survFitTT(dat, lcx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
#' 
#' # (4) Print the survFitTT object
#' out
#' }
#' 
#' @keywords print
#' 
#' @export
#' 
print.survFitTT <- function(x, ...) {
  # print the model text and the Jags Computing information
  # for an object of class survFitTT
  
  # M.C.M.C. informations
  cat("Model:\n")
  print(x$model)
  cat("\nComputing information:\n\n")
  cat("\n", "Iterations = ", x$n.iter[["start"]], ":",
      x$n.iter[["end"]], "\n", sep = "")
  cat("Thinning interval =", x$n.thin, "\n")
  cat("Number of chains =", x$n.chains, "\n")
  cat("Sample size per chain =",
      (x$n.iter[["end"]] - x$n.iter[["start"]]) / x$n.thin + 1, "\n")
}

#' Summary of survFitTT object
#' 
#' The summary shows the quantiles of priors and posteriors on parameters
#' and the quantiles of estimated LCx.
#' 
#' @param object An object of class survFitTT
#' @param quiet If \code{TRUE}, make silent all prints
#' 
#' @seealso survFitTT
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a survData object
#' cadmium1 <- survData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the survFitTT function with the log-logistic
#' # binomial model
#' out <- survFitTT(dat, lcx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
#' 
#' # (4) summarize the survFitTT object
#' summary(out)
#' }
#' 
#' @keywords summary
#' 
#' @export
#' 
summary.survFitTT <- function(object, quiet = FALSE) {
  
  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start
  
  # b
  log10b <- qunif(p = c(0.5, 0.025, 0.975),
                     min = object$jags.data$log10bmin,
                     max = object$jags.data$log10bmax)

  b <- 10^log10b
  
  # e
  log10e <- qnorm(p = c(0.5, 0.025, 0.975),
                  mean = object$jags.data$meanlog10e,
                  sd = 1 / sqrt(object$jags.data$taulog10e))

  e <- 10^log10e
 
  # d
  if (object$det.part == "loglogisticbinom_3") {

    d <- qunif(p = c(0.5, 0.025, 0.975),
               min = object$jags.data$dmin,
               max = object$jags.data$dmax)
    
    res <- rbind(b, d, e)
  } else {
    res <- rbind(b, e)
  }
  
  ans1 <- round(data.frame(res), digits = 3)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")
  
  # quantiles of estimated model parameters
  ans2 <- round(object$estim.par, digits = 3)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")
  
  # estimated ECx and their CIs 95%
  ans3 <- round(object$estim.LCx, digits = 3)
  colnames(ans3) <- c("50%", "2.5%", "97.5%")
  
  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("The ", object$det.part, " model with a binomial stochastic part was used !\n\n")
    cat("Quantiles of priors on parameters: \n\n")
    print(ans1)
    cat("\nQuantiles of posteriors on parameters: \n\n")
    print(ans2)
    cat("\nQuantiles of the estimated LCx: \n\n")
    print(ans3)
  }
  
  invisible(list(Qpriors = ans1,
                 Qpost = ans2,
                 QECx = ans3))
}

#' Fit a Bayesian exposure-response model for endpoint survival analysis
#'
#' The \code{survFitTT} function estimates the parameters of an exposure-response
#' model for endpoint survival analysis using Bayesian inference. In this model,
#' the survival rate of individuals after some time (called target time) is modeled
#' as a function of the pollutant's concentration, called deterministic part. The
#' actual number of
#' surviving individuals is then modeled as a stochastic function of the survival
#' rate, called the stochastic part. Details of the model are presented in the
#' vignette accompanying the package.
#'
#' The function returns
#' parameter estimates of the exposure-response model and estimates of the so-called
#' LCx, that is the concentration of pollutant required to obtain an 1 - x survival
#' rate.
#'
#'
#'
# \describe{
# FIXME
# Credible limits: For 100 values of concentrations regularly spread within
# the range of tested concentrations the joint posterior distribution of
# parameters is used to simulate 5000 values of \eqn{f_{ij}}, the number of
# offspring per individual-day for various replicates. For each concentration,
# 2.5, 50 and 97.5 percentiles of simulated values are calculated, from which
# there is a point estimate and a 95 \% credible interval (Delignette-Muller
# et al., 2014).
#
# DIC: The Deviance Information Criterium (DIC) as defined by Spiegelhalter et
# al. (2002) is provided by the \code{dic.samples} function. The DIC is a
# goodness-of-fit criterion penalized by the complexity of the model
# (Delignette-Muller et al., 2014).
#
# Raftery and Lewis's diagnostic: The \code{raftery.diag} is a run length
# control diagnostic based on a criterion that calculates the appropriate
# number of iterations required to accurately estimate the parameter
# quantiles. The Raftery and Lewis's diagnostic value used in the
# \code{surFitTT} function is the \code{resmatrix} object. See the
# \code{\link[coda]{raftery.diag}} help for more details.
#
#' @param data An object of class \code{survData}.
#' @param target.time The chosen endpoint to evaluate the effect of a given
#' concentration of pollutant. By default the last time point available for
#' all concentrations.
#' @param lcx Desired values of \eqn{x} (in percent) for which to compute
#' \eqn{LC_{x}}{LCx} (default is c(5, 10, 20, 50))
#' @param n.chains Number of MCMC chains. The minimum required number of chains
#' is 2.
#' @param quiet If \code{TRUE}, make silent all prints and progress bars of
#' JAGS compilation.
#' @param x An object of class \code{survFitTT}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Concentrations}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Response}.
#' @param main A main title for the plot.
#' @param fitcol A single color to plot the fitted curve, by default
#' \code{red}.
#' @param fitlty A single line type to plot the fitted curve, by default
#' \code{1}.
#' @param fitlwd A single numeric which controls the width of the fitted curve,
#' by default \code{1}.
#' @param ci If \code{TRUE}, the 95 \% confidente interval on observed data are
#' plotted.
#' @param conf.level : binom.test confidence level
#' @param cicol A single color to plot the 95 \% confidente interval, by default
#' \code{red}.
#' @param cilty A single line type to plot 95 \% confidente interval, by default
#' \code{1}.
#' @param cilwd A single numeric which controls the width of the 95 \% confidente
#' interval, by default \code{2}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param log.scale If \code{TRUE}, a log-scale is used on the \eqn{X}-axis.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns an object of class \code{survFitTT}, which is list containing
#' the following objects:
#' \item{DIC}{DIC value of the selected model.}
#' \item{estim.LCx}{A table of the estimated 5, 10, 20 and 50 \% lethal
#' concentrations and their 95 \% credible intervals.}
#' \item{estim.par}{A table of the estimated parameters as medians and 95 \%
#' credible intervals.}
#' \item{det.part}{The name of the deterministic part of the used model.}
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
#'
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#'
# FIXME
# @seealso \code{\link[rjags]{rjags}}, \code{\link[rjags]{coda.samples}},
# \code{\link[rjags]{dic.samples}}, 
# \code{\link{survData}}, \code{\link[coda]{raftery.diag}} and
# \code{\link[ggplot2]{ggplot}}
#
# @references Plummer, M. (2013) JAGS Version 3.4.0 user manual.
# \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}
#
# Spiegelhalter, D., N. Best, B. Carlin, and A. van der Linde (2002) Bayesian
# measures of model complexity and fit (with discussion).  \emph{Journal of
# the Royal Statistical Society}, Series B 64, 583-639.
#
#' @keywords estimation
#
#' @examples
#'
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create an object of class "survData"
#' dat <- survData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the survFitTT function with the log-logistic
#' # binomial model
#' out <- survFitTT(dat, lcx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
#'
#' # (4) Summary look the estimated values (LCx and parameters)
#' out$estim.LCx
#' out$estim.par
#'
#' # (5) Plot the fitted curve
#' plot(out, log.scale = TRUE, ci = TRUE)
#'
#' # (6) Plot the fitted curve with ggplot style
#' require(ggplot2)
#' plot(out, xlab = expression("Concentration in" ~ mu~g.L^{-1}),
#' fitcol = "blue", ci = TRUE, cicol = "blue",  style = "ggplot")
#'
#' # (7) Add a specific legend with generic type
#' plot(out, addlegend = FALSE)
#' legend("left", legend = c("Without mortality", "With mortality"),
#' pch = c(19,1))
#' }
#'
#' @export
#'
#' @import rjags
#' @importFrom dplyr filter
#'
survFitTT <- function(data,
                      target.time = NULL,
                      lcx,
                      n.chains = 3,
                      quiet = FALSE) {
  # test class object
  if(! is(data, "survData"))
    stop("survFitTT: object of class survData expected")

  # select Data at target.time
  dataTT <- selectDataTT(data, target.time)

  # Choose model by testing mortality in the control
  control <- filter(dataTT, conc == 0)
  det.part <-
    if (any(control$Nsurv < control$Ninit)) "loglogisticbinom_3"
  else "loglogisticbinom_2"

  # select model text
  if (det.part == "loglogisticbinom_2") {
    model.text <- llbinom2.model.text
  }
  if (det.part == "loglogisticbinom_3") {
    model.text <- llbinom3.model.text
  }

  # parameters
  parameters <- if (det.part == "loglogisticbinom_2") {
    c("log10b", "log10e")
  } else {
    if (det.part == "loglogisticbinom_3") {
      c("log10b", "d", "log10e")}
  }

  # create priors parameters
  jags.data <- survCreateJagsData(det.part, dataTT)

  # Define model
  model <- survLoadModel(model.program = model.text,
                         data = jags.data, n.chains,
                         Nadapt = 3000, quiet)

  # Determine sampling parameters
  sampling.parameters <- modelSamplingParameters(model,
                                                 parameters, n.chains, quiet)

  if (sampling.parameters$niter > 100000)
    stop("The model needs too many iterations to provide reliable parameter estimates !")

  # calcul DIC
  modelDIC <- calcDIC(model, sampling.parameters, quiet)

  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")

  mcmc <- coda.samples(model, parameters,
                       n.iter = sampling.parameters$niter,
                       thin = sampling.parameters$thin,
                       progress.bar = prog.b)

  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survPARAMS(mcmc, det.part)

  # LCx calculation  estimated LCx and their CIs 95%
  # vector of LCX
  if (missing(lcx)) {
    lcx <- c(5, 10, 20, 50)
  }
  estim.LCx <- estimXCX(mcmc, lcx, "LC")

  # check if estimated LC50 lies in the tested concentration range
  if (50 %in% lcx) {
    LC50 <- log10(estim.LCx["LC50", "median"])
    if (!(min(log10(data$conc)) < LC50 & LC50 < max(log10(data$conc))))
      warning("The LC50 estimation lies outsides the range of tested concentration and may be unreliable !")
  }

  # output
  OUT <- list(DIC = modelDIC,
              estim.LCx = estim.LCx,
              estim.par = estim.par,
              det.part = det.part,
              mcmc = mcmc,
              model = model,
              parameters = parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = data,
              dataTT = dataTT)

  class(OUT) <- "survFitTT"
  return(OUT)
}
