#' @importFrom stringr str_c
idCreate <- function(data, notime = FALSE) {
  # INPUT
  # data: a raw data frame
  # notime: if TRUE there is no column time in the dataset
  # OUTPUT
  # vector of triplet replicate_conc_time or couple replicate_conc
  if (!notime) {
    return(str_c(data[, "replicate"], data[, "conc"], data[, "time"],
                 sep = "_"))
  } else {
    return(str_c(data[, "replicate"], data[, "conc"], sep = "_"))
  }
}

#' @importFrom dplyr right_join %>% rename
#'
reproTransformData <- function(data) {

  # calculate Ninit, Nreprocumul and Nindtime
  # INPUTS
  # - data: a dataframe of class reproData with:
  #   - ID: ID replicate - conc - time
  #   - replicate: replicate indentification
  #   - conc: tested concentration
  #   - time: observation time
  #   - Nsurv: number of alive individuals at concentration "conc" and at time
  #     "time"
  #   - Nrepro: number of collected offspring at concentration "conc" and at time
  #     "time"
  #  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
  #   By default, it is the last time of original data.
  # OUTPUT: a dataframe with 9 columns corresponding to survival, Nindtime and
  # repro data at the time given in the input of the function
  #   - ID: ID replicate - conc - time
  #   - replicate: replicate indentification
  #   - conc: tested concentration
  #   - time: observation time
  #   - Ninit: number of individuals at the beginning of the bioassay
  #   - Nsurv: number of alive individuals at the time specified in the input of
  #     the function
  #   - Nrepro: number of offspring at this time
  #   - Nreprocumul: cumulative number of offspring at this time
  #   - Nindtime: number of individual-days from the beginning of the bioassay
  #     to this time

  # split dataset by time to calculate Ninit
  temp <- split(data, data$time)

  # add column Ninit
  tabletime0 <- data[data$time == 0, ] # control dataset
  tabletime0[, "Ninit"] <- tabletime0[, "Nsurv"]

  reproCalculNinit <- function(t1, t2) {
    # calcul the correct number Ninit
    # for each replicate conc and time
    # INPUTS
    # - t1: list of splited dataframe
    # - t2: tabletime0
    # OUTPUTS
    # - list of splited dataframe with the news column Ninit

    . = NULL
    ID.x = NULL
    time.x = NULL
    Nsurv.x = NULL
    Nrepro.x = NULL

    right_join(t1, t2,
               by = c("replicate", "conc"))[, c("ID.x", "replicate", "conc", "time.x",
                                                "Nsurv.x", "Nrepro.x", "Ninit")] %>% rename(., time = time.x) %>% rename(.,
                                                                                                                         Nsurv = Nsurv.x) %>% rename(., Nrepro = Nrepro.x) %>% rename(., ID = ID.x)
  }

  res <- lapply(temp, function(x) reproCalculNinit(x, tabletime0)) # Ninit

  data <- do.call(rbind, res) # return to a dataframe
  rownames(data) <- 1:dim(data)[1]
  data$time <- as.numeric(data$time) # change type of time

  T <- unique(data$time) # times of obs without repetitions
  finalnbr <- match(max(data$time), T) # index of the time at which we want the
  #  results in vector T
  if (finalnbr == 1) {
    stop("!!!! It isn't possible to use the first observation time as the last observation time !!!!")
  }

  # calculation of cumulative number of repro and individual-days

  tableTi = list()

  for (i in 2:finalnbr) {
    # original dataset at T[i-1]
    dataTim1 <- subset(data, data$time == T[i-1])
    # original dataset at T[i]
    dataTi <- subset(data, data$time == T[i])
    # check if data have been properly captured
    if (any(dataTim1$replicate != dataTi$replicate) || any(dataTim1$conc != dataTi$conc))
      warning("!!!! BE CAREFUL concentrations and/or replicates are not identical at each time !!!!")

    # build the output date at time T[i]
    if (i == 2) {
      # (number of survivors at T[i]) * T[i]
      # + (number of dead organisms between T[i] and T[i-1]) * (T[i]+T[i-1])/2
      NindtimeTi <- (dataTi$Nsurv*dataTi$time +
                       (dataTim1$Nsurv - dataTi$Nsurv)*(dataTi$time + dataTim1$time)/2)
      NreprocumulTi <- dataTim1$Nrepro + dataTi$Nrepro
    } else {
      # (number of survivors at T[i]) * T[i]
      # + (number of dead organisms between t and tm1) * (t+tm1)/2
      # + number of individual-time at T[i-1] (accounting for dead organisms
      # at time before T[i-1])
      # - number of individual-time at T[i-1] that are staying alive at T[i-1]
      NindtimeTi <- ( dataTi$Nsurv*dataTi$time +
                        (dataTim1$Nsurv - dataTi$Nsurv)*(dataTi$time + dataTim1$time)/2 +
                        tableTim1$Nindtime - dataTim1$Nsurv * dataTim1$time )

      NreprocumulTi <- tableTim1$Nreprocumul + dataTi$Nrepro
    }

    tableTi[[i]] <- data.frame(ID = dataTi$ID,
                               replicate = dataTi$replicate,
                               conc = dataTi$conc,
                               time = dataTi$time,
                               Ninit = dataTi$Ninit,
                               Nsurv = dataTi$Nsurv,
                               Nrepro = dataTi$Nrepro,
                               Nreprocumul = NreprocumulTi, # cumulative number of offspring
                               Nindtime = NindtimeTi)

    # tabelTi stored as tableTim1 for next iteration
    tableTim1 <- tableTi[[i]]
  }

  tablefinale <- do.call("rbind", tableTi)
  return(tablefinale)
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

logTransXaxisFit <- function(log.scale, concentrations) {
  # remove concentrations 0 in log.scale if needed
  # INPUT :
  # - log.scale: booleen
  # - concentrations: vector of concentrations
  # OUTPUT
  # - X: new vector of concentration of length 100
  
  if (log.scale && any(concentrations == 0)) { # log.scale yes
    # 0 in the concentrations
      X <- seq(min(unique(concentrations)[-1]), max(concentrations),
               length = 100)
    } else { # no 0 in the concentrations
      X <- seq(min(concentrations), max(concentrations), length = 100)
    }
  return(X)
}

survLlbinomFit <- function(res.M, x, X) {
  # create the parameters of the fitted curve for loglogistic model
  # INPUT :
  # - res.M: mcmc summary
  # - x: repro.fit object
  # - X: vector of concentrations (xaxis)
  # OUTPUT :
  # - fNsurvtheo
  
  # unlog parameters
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  
  if (x$det.part == "loglogisticbinom_3") {
    d <- res.M$quantiles["d", "50%"]
    fNsurvtheo <- d / (1 + (X / e)^b) # mean curve equation 3 parameters
  } else {
    fNsurvtheo <- 1 / (1 + (X / e)^b) # mean curve equation 2 parameters
  }
  return(fNsurvtheo)
}

survLlbinomCi <- function(x, X) {
  # create the parameters for credible interval for the log logistic model
  # INPUT:
  # - x : object of class repro.fit
  # - X : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit
  
  mctot <- do.call("rbind", x$mcmc)
  k <- nrow(mctot)
  
  # parameters
  if (x$det.part == "loglogisticbinom_3"){
    d2 <- mctot[,"d"]
  }
  log10b2 <- mctot[,"log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[,"log10e"]
  e2 <- 10^log10e2
  
  # quantiles
  qinf95 = NULL
  qsup95 = NULL
  
  for (i in 1:length(X)) {
    if (x$det.part == "loglogisticbinom_3") {
      theomean <- d2/(1 + (X[i] / e2)^(b2))
    } else {
      theomean <- 1/(1 + (X[i] / e2)^(b2))
    }
    
    # IC 95%
    qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
    qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
  }
  # values for CI
  ci <- list(qinf95 = qinf95,
             qsup95 = qsup95)
  return(ci)
}

logTransConcFit <- function(log.scale, X, x, concentrations) {
  # transform vector conccentration in log.scale if needed
  # X axis log.scale value
  
  # select non 0 values
  sel <- if (log.scale) {
    X > 0
  } else {
    rep(TRUE, length(X))
  }
  
  if (log.scale) {
    X[sel] <- log(X[sel])
  } else {
    X[sel] <- X[sel]
  }
  
  # log transform concentrations values
  sel2 <- if (log.scale) {
    x$dataTt$conc > 0
  } else {
    rep(TRUE, length(x$dataTt$conc))
  }
  if (log.scale) {
    concentrations[sel2] <- log(concentrations[sel2])
  } else {
    concentrations[sel2] <- concentrations[sel2]
  }
  
  return(list(sel = sel,
              X = X,
              sel2 = sel2,
              concentrations = concentrations))
}

legendGgplotFit <- function(a.gplot) {
  # create multi legend for plot.reproFitTt and plot.survFitTt
  # to have differente legend for points and mean and CI curve
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}
