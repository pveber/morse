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

#' @import dclone
survLoadModelPar <- function(cl,
                             model.program,
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
  model <- parJagsModel(cl, name = "res", file = model.file, data = data,
                        n.chains = n.chains, n.adapt = Nadapt, quiet = quiet)
  unlink(model.file)
  return(model)
}

#' @import dclone
#' @importFrom coda raftery.diag
modelSamplingParametersPar <- function(cl, model, parameters, n.chains,
                                       quiet = quiet) {
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
  mcmc <- parCodaSamples(cl,"res", parameters, n.iter = niter.init, thin = 1,
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