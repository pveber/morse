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

survCreateJagsData <- function(det.part, data) {
  # create the parameters to define the prior of the log-logistic binomial model
  # INPUTS
  # det.part: model name
  # data: object of class survData
  # OUTPUT
  # jags.data : list data require for the jags.model function
  
  # Parameter calculation of concentration min and max
  concmin <- min(sort(unique(data$conc))[-1])
  concmax <- max(data$conc)
  
  # create priors parameters for the log logistic model
  
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
    dmin = 0
    dmax = 1
    jags.data <- c(jags.data,
                   dmin = dmin,
                   dmax = dmax)
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

survFullPlotGeneric <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate
  # for generic type
  # INPUTS
  # - data: raw dataframe with 4 columns with:
  #   - replicate: replicate indentification
  #   - conc: tested concentrations
  #   - time: time of the observation
  #   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
  # OUTPUT:
  # - Plot

    .convert <- function(x) {
      # conversion of a replicate name in a number coding for color
      # INPUT
      # - x: name of replicate
      # OUTPUT
      #  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
      # - a position of replicate in vector to display color
      replicate <- unique(data$replicate)
      mat <- matrix(ncol = 2, nrow = length(replicate))
      mat[, 1] <- as.character(replicate)
      mat[, 2] <- as.character(seq(1, length(replicate)))
      return(as.integer(mat[mat[, 1] == as.character(x)][2]))
    }

    .NbPlot <- function(conc) {
      # definition of the number of subplots
      # INPUT
      # - conc: vector of tested concentrations
      # OUTPUT
      # - vector defining the number of columns and rows in par(mfrow)

      nbconc <- length(conc)
      PlotPar <- c(c(2, 2), c(2, 3), c(2, 4), c(3, 3), c(2, 5), c(3, 4), c(3, 5),
                   c(4, 4))
      NbPlotTheo <- matrix(ncol = 2, nrow = 8)
      NbPlotTheo[, 1] <- c(1, 3, 5, 7, 9, 11, 13, 15)
      NbPlotTheo[, 2] <- c(4, 6, 8, 9, 10, 12, 15, 16)
      if (nbconc < 15) {
        PlotTrue <- NbPlotTheo[NbPlotTheo[, 2] - nbconc > 0, ][1, 1]
      } else {
        PlotTrue = 15
      }
      return(c(PlotPar[PlotTrue], PlotPar[PlotTrue + 1]))
    }


    # creation of a vector of colors
    colors <- rainbow(length(unique(data$replicate)))
    pchs <- as.numeric(unique(data$replicate))
    # split of the graphical window in subplots
    par(mfrow = .NbPlot(unique(data$conc)))

    by(data, data$conc, function(x) {
      # bakground
      plot(x$time, rep(0, length(x$time)),
           xlab = xlab,
           ylab = ylab,
           ylim = c(0,max(x$Nsurv)),
           type = "n",
           col = 'white',
           yaxt = 'n')

      # axis
      axis(side = 2, at = pretty(c(0,max(x$Nsurv))))

      # lines and points
      by(x, x$replicate, function(y) {
        lines(y$time, y$Nsurv,
              type = "l",
              col = colors[.convert(unique(y$replicate))])
        points(y$time, y$Nsurv,
               pch = pchs[.convert(unique(y$replicate))],
               col = colors[.convert(unique(y$replicate))])
      })

      # title
      title(paste("Conc: ", unique(x$conc), sep = ""))
    })

    if (addlegend) {
      # creation of an empty plot to display legend
      plot(0, 0,
           xlab = "",
           ylab = "",
           xlim = c(0,5),
           ylim = c(0,5),
           type = "n",
           xaxt = "n",
           yaxt = "n")

      # Display legend
      title.legend <- "Replicate"
      mat <- matrix(nrow = length(unique(data$replicate)), ncol = 2)
      mat[, 1] <- rep(title.legend, length(unique(data$replicate)))
      mat[, 2] <- unique(as.character(data$replicate))
      name <- apply(mat, 1, function(x) {paste(x[1], x[2], sep = ": ")})

      legend("top", name,
             lty = rep(1, length(unique(data$replicate))),
             pch = pchs,
             col = colors,
             bty = "n",
             cex = 1)
    }
    par(mfrow = c(1, 1))
}

#' @importFrom lattice lattice.options xyplot
survFullPlotL <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate for lattice type
  # INPUTS
  # - data: raw dataframe with 4 columns with:
  #   - replicate: replicate indentification
  #   - conc: tested concentrations
  #   - time: time of the observation
  #   - Nsurv: number of alive individuals at time "time" and at concentration
  #     "conc"
  # OUTPUT:
  # - Plot

  # change order of reading concentrations
  lattice.options(default.args = list(as.table = TRUE))

  if (addlegend) {
    xyplot(Nsurv ~ time | factor(conc),
           data = data,
           group = replicate,
           index.cond = list(c((round(length(unique(data$conc)) / 2) + 1):length(unique(data$conc)),
                               1:(round(length(unique(data$conc)) / 2)))),
           type = "b",
           pch = 16,
           xlab = xlab,
           ylab = ylab,
           auto.key = list(space = "right",
                           title = "Replicate",
                           cex.title = 1,
                           lines = TRUE,
                           points = FALSE))
  } else {
    xyplot(Nsurv ~ time | factor(conc),
           data = data,
           group = replicate,
           index.cond = list(c((round(length(unique(data$conc)) / 2) + 1):length(unique(data$conc)),
                               1:(round(length(unique(data$conc))/2)))),
           type = "b",
           pch = 16,
           xlab = xlab,
           ylab = ylab)
  }
}

#' @import ggplot2
survFullPlotGG <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate for ggplot type
  # INPUTS
  # - data: raw dataframe with 4 columns with:
  #   - replicate: replicate indentification
  #   - conc: tested concentrations
  #   - time: time of the observation
  #   - Nsurv: number of alive individuals at time "time" and at concentration "conc"
  # OUTPUT:
  # - Plot

  time = NULL
  Nsurv = NULL
  title.legend <- "Replicate"

  # create ggplot object Nsurv / time / replicate / conc
  fg <- ggplot(data, aes(time, Nsurv, colour = factor(replicate))) +
    geom_point() +
    geom_line() +
    labs(x = xlab, y = ylab) +
    facet_wrap(~conc, nrow = 2) +
    scale_x_continuous(breaks = unique(data$time)) +
    ylim(0, max(data$Nsurv))

  # legend option
  if (addlegend){
    fd <- fg + scale_colour_hue(title.legend) # the default legend
  } else {
    fd <- fg + theme(legend.position = "none") # remove legend
  }
  return(fd)
}

llbinom3.model.text <- "\nmodel # Loglogistic binomial model with 3 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- d/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nd ~ dunif(dmin, dmax)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"

llbinom2.model.text <- "\nmodel # Loglogistic binomial model with 2 parameters\n\t\t{\t\nfor (i in 1:n)\n{\np[i] <- 1/ (1 + (xconc[i]/e)^b)\nNsurv[i]~ dbin(p[i], Ninit[i])\n}\n\n# specification of priors (may be changed if needed)\nlog10b ~ dunif(log10bmin, log10bmax)\nlog10e ~ dnorm(meanlog10e, taulog10e)\n\nb <- pow(10, log10b)\ne <- pow(10, log10e)\n}\n"

