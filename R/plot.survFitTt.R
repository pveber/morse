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

survLlbinomCi <- function(x) {
  # create confidente interval on observed data for the log logistic
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # OUTPUT:
  # - ci : confidente interval
  x <- cbind(aggregate(Nsurv ~ time + conc, x$dataTT, sum),
             Ninit = aggregate(Ninit ~ time + conc, x$dataTT, sum)$Ninit)
  
  ci <- apply(x, 1, function(x) {
    binom.test(x["Nsurv"], x["Ninit"])$conf.int
    })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$conc
  
  return(ci)
}

survFitPlotGenericNoCi <- function(concentrations, response, x, X,
                                   fNsurvtheo, sel2, sel, mortality, 
                                   xlab, ylab, fitcol, fitlty, fitlwd,
                                   main, addlegend, legend.position, legend.title,
                                   legend.name.no, legend.name.yes, ...)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style without credible interval
  plot(concentrations[sel2], response[sel2],
       xlab = xlab,
       ylab = ylab,
       main = main,
       pch = mortality,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1.2),
       ...)
  
  # axis
  axis(side = 2, at = c(0, 1))
  axis(side = 1, at = unique(concentrations[sel2]),
       labels = unique(x$dataTT$conc[sel2]))
  
  # fitted curve
  lines(X[sel], fNsurvtheo[sel], col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  # legend
  if (addlegend) {
    legend(legend.position, title = legend.title, pch = c(19, 1, NA),
           lty = c(0, 0, fitlty),
           lwd = c(1, 1, fitlwd),
           col = c(1, 1, fitcol),
           legend = c(legend.name.no, legend.name.yes, x$det.part),
           bty = "n")
  }
  }

survFitPlotGenericCi <- function(concentrations, response, x, X,
                                 fNsurvtheo, CI, sel3, sel2, sel, mortality, 
                                 xlab, ylab, fitcol, fitlty, fitlwd,
                                 main, addlegend, legend.position, legend.title,
                                 legend.name.no, legend.name.yes, legend.position.ci,
                                 cicol, cilty, cilwd, ...)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style with credible interval
  plot(concentrations[sel2], response[sel2],
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1.2),
       type = "n",
       ...)
  
  # axis
  axis(side = 2, at = pretty(c(0, 1.2)))
  axis(side = 1, at = unique(concentrations[sel2]),
       labels = unique(x$dataTT$conc[sel2]))

  # Plotting the theoretical curve

  # fitted curve
  lines(X[sel], fNsurvtheo[sel], col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  # points
  points(concentrations[sel2], response[sel2], pch = mortality)
  
  # segment CI
  segments(unique(concentrations[sel2]),
           CI["qinf95",][sel3], concentrations[sel2],
           CI["qsup95", ][sel3], col = cicol, lty = cilty, lwd = cilwd)
  
  # legend
  if (addlegend) { # legend yes CI yes
    legend(legend.position, title = legend.title, pch = c(19, 1, NA, NA),
           lty = c(0, 0, fitlty, cilty),
           lwd = c(1, 1, fitlwd, cilwd),
           col = c(1, 1, fitcol, cicol),
           legend = c(legend.name.no, legend.name.yes,
                      x$det.part, "Confidence interval"),
           bty = "n")
  }
}

survFitPlotGGNoCi <- function(data.one, data.two, valCols, 
                              fitlty, fitlwd, xlab, ylab, main) {
  plt_4 <- ggplot(data.one) +
    geom_point(data=data.one, aes(concentrations.sel2., response.sel2.,
                                  color = mortality)) +
    geom_line(aes(X.sel., fNsurvtheo.sel.), data.two,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(plt_4)
}

survFitPlotGGCi <- function(X, x, CI, sel, data.one, data.two, cilty, cilwd,
                            valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = unique(concentrations[sel2]),
                           qinf95 = CI["qinf95", ][sel3],
                           qsup95 = CI["qsup95", ][sel3],
                           Ci = paste("Confidence interval"))
  
  plt_3 <- ggplot(data.one) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95),
                 data.three, col = valCols$cols3, linetype = cilty,
                 dize = cilwd)
  
  # plot IC
  # final plot
  plt_4 <- ggplot(data.one) +
    geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                    color = mortality)) +
    geom_line(aes(X.sel., fNsurvtheo.sel.), data.two, linetype = fitlty,
              size = fitlwd, color = valCols$cols2) +
    geom_line(aes(X.sel., CI.qinf95.sel.), data.three, linetype = cilty,
              size = cilwd, color = valCols$cols3) +
    geom_line(aes(X.sel., CI.qsup95.sel.), data.three, linetype = cilty,
              size = cilwd, color = valCols$cols3) +
    geom_ribbon(data = data.three, aes(x = X.sel., ymin = CI.qinf95.sel.,
                                       ymax = CI.qsup95.sel.),
                alpha = 0.4) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

#' @export
#' 
#' @name survFitTT
#' 
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTT <- function(x,
                           xlab,
                           ylab,
                           main,
                           fitcol,
                           fitlty,
                           fitlwd,
                           ci = FALSE,
                           cicol,
                           cilty,
                           cilwd,
                           addlegend = TRUE,
                           log.scale = FALSE,
                           style = "generic",
                           pool.replicate = TRUE, ...) {
  # plot the fitted curve estimated by survFitTT
  # INPUTS
  # - x:  survFitTt object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - ci : credible interval, boolean
  # - cicol : color ci
  # - cilty : type line ci
  # - cilwd : width line ci
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic or ggplot
  # OUTPUT:
  # - plot of fitted regression	
  
  if (style == "ggplot") {
    # variables declarations
    concentrations.sel2. = NULL
    response.sel2. = NULL
    X.sel. = NULL
    fNsurvtheo.sel. = NULL
    qinf95.sel. = NULL
    qsup95.sel. = NULL
    CI.qinf95.sel. = NULL
    CI.qsup95.sel. = NULL
    Line = NULL
    Ci = NULL
  }
  
  # Define data
  concentrations <- x$dataTT$conc
  response <- x$dataTT$Nsurv / x$dataTT$Ninit
  
  if (pool.replicate) {
    resp.temp <- cbind(response, x$dataTT$conc)
    response <- aggregate(resp.temp, by = list(resp.temp[,2]), mean)$response
  }
  
  # Fitted curve parameters
  
  X <- logTransXaxisFit(log.scale, concentrations) # X log.scale transformation
  
  res.M <- summary(x$mcmc)
  
  # Choose median of posteriors as parameter values	
  fNsurvtheo <- survLlbinomFit(res.M, x, X)
  
  # calculate IC 95 values
  if (ci) {
    CI <- survLlbinomCi(x)
  }
  
  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  sel <- logTransConcFit(log.scale, X, x, concentrations, "sel")
  X <- logTransConcFit(log.scale, X, x, concentrations, "X")
  sel2 <- logTransConcFit(log.scale, X, x, concentrations, "sel2")
  sel3 <- logTransConcFit(log.scale, X, x, concentrations, "sel3")
  concentrations <- logTransConcFit(log.scale, X, x, concentrations, "conc")

  nomortality <- match(x$dataTT$Nsurv[sel2] == x$dataTT$Ninit[sel2],
                       c(TRUE, FALSE)) # valid if at least one replicat
  
  if (pool.replicate) {
    nomortality.temp <- cbind(nomortality, x$dataTT$conc[sel2])
    nomortality <- aggregate(nomortality.temp, by = list(nomortality.temp[,2]),
                             mean)$nomortality
    nomortality[nomortality != 1] <- 2
  }
  
  # without mortality
  mortality <- mortality[nomortality] # vector of 0 and 1
  
  # encodes mortality empty dots (1) and not mortality solid dots (19)
  if (style == "generic") {
    mortality[which(mortality == 0)] <- 19
  }
  if (style == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }
  
  # default axis parameters
  if (missing(xlab)) {
    xlab <- "Concentrations"
  }
  if (missing(ylab)) {
    ylab <- "Response"
  }
  
  # default legend parameters
  if (missing(fitcol)) {
    fitcol <- "red"
  }
  if (missing(fitlty)) {
    fitlty <- 1
  }
  if (missing(fitlwd)) {
    fitlwd <- 1
  }
  if (missing(main)) {
    main = NULL
  }
  if (style == "generic") {
    legend.position <- "bottomleft"
    legend.position.ci <- "left"
  }
  legend.title <- "Mortality"
  legend.name.no <- "No"
  legend.name.yes <- "Yes"
  
  # IC parameters
  if (missing(cicol)) {
    cicol <- "red"
  }
  if (missing(cilty)) {
    cilty <- 2
  }
  if (missing(cilwd)) {
    cilwd <- 1
  }
  
  # Plotting data
  if (style == "generic") {
    if (!ci) {
      survFitPlotGenericNoCi(concentrations, response, x, X,
                             fNsurvtheo, sel2, sel, mortality, 
                             xlab, ylab, fitcol, fitlty, fitlwd,
                             main, addlegend, legend.position, legend.title,
                             legend.name.no, legend.name.yes, ...)
    }
    if (ci) {
      survFitPlotGenericCi(concentrations, response, x, X,
                           fNsurvtheo, CI, sel3, sel2, sel, mortality, 
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend, legend.position, legend.title,
                           legend.name.no, legend.name.yes, legend.position.ci,
                           cicol, cilty, cilwd, ...)
    }

  }
  
  if (style == "ggplot") {
    if (Sys.getenv("RSTUDIO") == "") {
      dev.new() # create a new page plot
      # when not use RStudio
    }
    
    # dataframes points (one) and curve (two)
    data.one <- data.frame(conc = concentrations[sel2],
                           response = response[sel2], mortality)
    data.two <- data.frame(X = X[sel], fNsurvtheo = fNsurvtheo[sel],
                           Line = x$det.part)
    
    # colors
    valCols <- fCols(data.one, x, fitcol, cicol)
    
    # points (to create the legend)
    plt_1 <- ggplot(data.one) +
      geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                      color = mortality)) +
      scale_color_manual(values = valCols$cols1)
    
    # curve (to create the legend)
    plt_2 <- ggplot(data.one) +
      geom_line(data = data.two, aes(X.sel., fNsurvtheo.sel., color = Line),
                linetype = fitlty, size = fitlwd) +
      scale_color_manual(values = valCols$cols2)
    
    if (ci) { # IC yes
      plt_3 <- survFitPlotGGCi(X, x, CI, sel, data.one, data.two, cilty, cilwd,
                               valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3

      plt_4 <- survFitPlotGGCi(X, x, CI, sel, data.one, data.two, cilty, cilwd,
                               valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
    }
    if (!ci) { # IC no
      plt_4 <- survFitPlotGGNoCi(data.one, data.two, valCols, fitlty, fitlwd,
                                 xlab, ylab, main)
    }
    
    if (addlegend) { # legend yes
      # create legends
      mylegend_1 <- legendGgplotFit(plt_1) # points legend
      mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
      
      if (ci) mylegend_3 <- legendGgplotFit(plt_3) # CI legend
      if (log.scale) { # log.sclae yes
        plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations.sel2.),
                                            labels =  unique(x$dataTT$conc[sel2]))
      } else { # log.scale no
        plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations.sel2.))
      }
      if (!ci) { # CI no
        grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, nrow = 6),
                     ncol = 2, widths = c(7,1))
      }
      if (ci) { # CI yes
        grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, mylegend_3,
                                        nrow = 6), ncol = 2, widths = c(7,1))
      }
    } else { # legend no
      if (log.scale) { # log.scale yes
        plt_5 <- plt_4 +
          scale_x_continuous(breaks = unique(data.one$concentrations.sel2.),
                             labels = unique(x$dataTT$conc[sel2]))
      } else { # log.scale no
        plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations.sel2.))
      }
      return(plt_5)
    }
  }
}
