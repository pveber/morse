# Ugly hack to get rid of spurious notes in package check, caused by uses
# of dplyr::{rename, filter}. R is such a sad language.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Nsurv","conc","Ninit"))


# Generates a character string vector from a data.frame using its replicate,
# conc and time columns. The result can be used as identifiers for the rows
# of the data.set.
#
#' @importFrom stringr str_c
#'
idCreate <- function(data) {
  str_c(data[, "replicate"],
        data[, "conc"],
        data[, "time"],
        sep = "_")
}

#' @importFrom dplyr filter
selectDataTt <- function(data, target.time) {
  # INPUT
  # - data: An object of class reproData or survData
  # - target.time: the time we want to consider as the last time for the analysis.
  # OUTPUT
  # - subset the dataframe for target.time used by function with target.time arg
  
  # target.time default
  if ("time" %in% colnames(data)) { # one time dataset check
    if (is.null(target.time)) {
      target.time <- max(data$time)
    }
    
    # correct target time
    if (!any(data$time == target.time))
      stop("target.time is not one of the possible time !")
    
    datatt <- filter(data, time == target.time)
  } else {
    datatt <- cbind(data, time = 1)
  }
  
  return(datatt)
}

#' @import rjags
#' @importFrom coda raftery.diag
modelSamplingParameters <- function(model, parameters, n.chains, quiet = quiet) {
  # estimate the number of iteration required for the estimation
  # by using the raftery.diag
  # INPUTS:
  # - model: jags model from loading function
  # - parameters: parameters from loading function
  # - nchains: Number of chains desired
  # - quiet: silent option
  # OUTPUTS:
  # - niter: number of iteration (mcmc)
  # - thin: thining rate parameter
  # - burnin: number of iteration burned
  
  # number of iteration for the pilote run required by raftery.diag
  # default value: 3746
  niter.init <- 5000
  prog.b <- ifelse(quiet == TRUE, "none", "text") # plot progress bar option
  mcmc <- coda.samples(model, parameters, n.iter = niter.init, thin = 1,
                       progress.bar = prog.b)
  RL <- raftery.diag(mcmc)
  
  # check raftery.diag result (number of sample for the diagnostic procedure)
  if (n.chains < 2) stop('2 or more parallel chains required !')
  
  # extract raftery diagnostic results
  resmatrix <- RL[[1]]$resmatrix
  for (i in 2: length(RL)) {
    resmatrix <- rbind(resmatrix, RL[[i]]$resmatrix)
  }
  
  # creation of sampling parameters
  thin <- round(max(resmatrix[, "I"]) + 0.5) # autocorrelation
  niter <- max(resmatrix[, "Nmin"]) * thin # number of iteration
  burnin <- max(resmatrix[, "M"]) # burnin period
  
  return(list(niter = niter, thin = thin, burnin = burnin))
}

#' @import rjags
calcDIC <- function(m.M, sampling.parameters, quiet = quiet) {
  # calculate the dic for a jags model
  # INPUTS
  # - m.M:  jags model object
  # - niter: number of iterations for the sampling
  # - thin
  # OUTPUT:
  # - numeric value of the DIC
  
  prog.b <- ifelse(quiet == TRUE, "none", "text") # plot progress bar option
  
  # estimation of DIC
  dic <- dic.samples(m.M, n.iter = sampling.parameters$niter,
                     thin = sampling.parameters$thin, progress.bar = prog.b)
  
  # return penalised DIC
  return(round(sum(sapply(dic$deviance, mean) + sapply(dic$penalty, mean))))
}

estimXCX <- function(mcmc, xcx, varx) {
  # create the table of estimated values of LCx or ECx
  # for the survival analyses
  
  # INPUT:
  # - mcmc:  list of estimated parameters for the model with each item representing
  # a chains
  # - xcx: vector of values of LCx or ECX
  # - varx: character string for lcx or ecx
  # OUTPUT:
  # - data frame with the estimated ECx and their CIs 95% (3 columns (values,
  # CIinf, CIsup) and length(x) rows)
  
  # Retrieving estimated parameters of the model
  mctot <- do.call("rbind", mcmc)
  b <- 10^mctot[, "log10b"]
  e <- 10^mctot[, "log10e"]
  
  # Calculation XCx median and quantiles
  XCx <- sapply(xcx, function(x) {e * ((100 / (100 - x)) - 1)^(1 / b)})
  
  q50 <- apply(XCx, 2, function(XCx) {quantile(XCx, probs = 0.5)})
  qinf95 <- apply(XCx, 2, function(XCx) {quantile(XCx, probs = 0.025)})
  qsup95 <- apply(XCx, 2, function(XCx) {quantile(XCx, probs = 0.975)})
  
  # defining names
  XCname <- sapply(xcx, function(x) {paste(varx, x, sep = '')})
  colnames(XCx) <- XCname
  
  # create the dataframe with ECx median and quantiles
  res <- data.frame(median = q50, Q2.5 = qinf95, Q97.5 = qsup95,
                    row.names = XCname)
  
  return(res)
}
