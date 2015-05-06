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


#' 
#' @name survFitTt
#' 
#' @export
#' @import dclone
#' @import rjags
#' @importFrom dplyr filter
#' @importFrom parallel makePSOCKcluster stopCluster
#' 
survParfitTt <- function(x,
                         det.part = "loglogisticbinom_2",
                         target.time = NULL,
                         lcx,
                         n.chains = 3,
                         quiet = FALSE) {
  
  requireNamespace("snow")
  
  # test class object
  if(! is(x,"survData"))
    stop("survFitTt: object of class survData expected")
  
  # test determinist part
  if (!any(det.part == "loglogisticbinom_2" || det.part == "loglogisticbinom_3"))
    stop("The [det.part] argument is not a possible deterministic part !")
  
  # open cluster
  cl <- makePSOCKcluster(n.chains)
  
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
      c("d","log10b", "log10e")}
  }
  
  # select Data at target.time
  dataTt <- selectDataTt(x, target.time)
  
  # create priors parameters
  jags.data <- survCreateJagsData(det.part, dataTt)
  
  # Test mortality in the control
  if (any(filter(dataTt,
                 conc == 0)$Nsurv < filter(dataTt,
                                           conc == 0)$Ninit) && det.part == "loglogisticbinom_2")
    stop("Beware! There is mortality in the control. A model with three parameters must be chosen.")
  
  if (!any(filter(dataTt,
                  conc == 0)$Nsurv < filter(dataTt,
                                            conc == 0)$Ninit) && det.part == "loglogisticbinom_3")
    stop("Beware! There is no mortality in the control. A model with two parameters must be chosen.")
  
  #	Model computing
  
  # Define model
  model <- survLoadModelPar(cl, model.program = model.text,
                            data = jags.data, n.chains,
                            Nadapt = 3000, quiet)
  
  # Determine sampling parameters
  sampling.parameters <- modelSamplingParametersPar(cl, model,
                                                    parameters, n.chains, quiet)
  
  if(sampling.parameters$niter > 100000)
    stop("The model needs too many iterations to provide reliable parameter estimates !")
  
  # calcul DIC
  modeldic2 <- tempfile() # temporary file address
  fileC <- file(modeldic2) # open connection
  writeLines(model.text, fileC) # write text in temporary file
  close(fileC) # close connection to temporary file
  modeldic3 <- jags.model(file = modeldic2,
                          data = jags.data,
                          n.chains = 2,
                          n.adapt = 3000, quiet = quiet)
  modelDIC <- calcDIC(modeldic3, sampling.parameters, quiet)
  
  # Sampling
  prog.b <- ifelse(quiet == TRUE, "none", "text")
  
  mcmc <- parCodaSamples(cl, "res", parameters,
                         n.iter = sampling.parameters$niter,
                         thin = sampling.parameters$thin,
                         progress.bar = prog.b)
  
  # close cluster connection
  stopCluster(cl)
  
  # summarize estime.par et CIs
  # calculate from the estimated parameters
  estim.par <- survPARAMS(mcmc, det.part)
  
  # LCx calculation  estimated LCx and their CIs 95%
  # vector of LCX
  if (missing(lcx)) {
    lcx <- c(5, 10, 20, 50)
  }
  estim.LCx <- estimXCX(mcmc, lcx, "LC")
  
  # check if the maximum measured concentration is in the LC50's range of
  # 95% percentile
  if (50 %in% lcx) {
    if (!(min(log10(x$conc)) < log10(estim.LCx["LC50", "median"]) &
          log10(estim.LCx["LC50", "median"]) < max(log10(x$conc))))
      warning("The LC50 estimation lies outsides the range of tested concentration and may be reliable !")
  }
  
  # output
  OUT <- list(DIC = modelDIC,
              estim.LCx = estim.LCx,
              estim.par = estim.par,
              det.part = det.part,
              mcmc = mcmc,
              model = modeldic3,
              parameters = parameters,
              n.chains = summary(mcmc)$nchain,
              n.iter = list(start = summary(mcmc)$start,
                            end = summary(mcmc)$end),
              n.thin = summary(mcmc)$thin,
              jags.data = jags.data,
              transformed.data = x,
              dataTt = dataTt)
  
  class(OUT) <- "survFitTt"
  return(OUT)
}
