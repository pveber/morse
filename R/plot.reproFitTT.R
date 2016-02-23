#' Plotting method for \code{reproFitTT} objects
#' 
#' This function plots exposure-response fits for target time reproduction
#' analysis (a.k.a. \code{reproFitTT} objects).
#' 
#' The fitted curve represents the \strong{estimated reproduction rate} after
#' the target time has passed as a function of the concentration of pollutant;
#' When \code{adddata = TRUE} the black dots depict the \strong{observed reproduction
#' rate} at each tested concentration.
#' The function plots both 95 \% credible intervals for the estimated reproduction
#' rate (by default the red area around the fitted curve) and 95 \% confidence
#' intervals for the observed reproduction rate (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#' If spaghetti = TRUE, the credible intervals are represented by two dotted
#' lines limiting the credible band, and a spaghetti plot is added to this band.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (10 \% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x an object of class \code{reproFitTT}
#' @param xlab a title for the \eqn{x}-label
#' @param ylab a title for the \eqn{y}-label
#' @param main main title for the plot
#' @param fitcol color used for the fitted curve
#' @param fitlty line type for the fitted curve
#' @param fitlwd width of the fitted curve
#' @param spaghetti if \code{TRUE}, the credible interval is drawn by  multiple
#' curves
#' @param cicol color for the 95 \% credible limits of the fitted curve
#' @param cilty line type for the 95 \% credible limits of the fitted curve
#' @param cilwd width of the 95 \% credible limits of the fitted curve
#' @param adddata if \code{TRUE}, adds the observed data with confidence interval
#' to the plot
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log-scale
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @note When \code{style = "ggplot"}, the function calls function
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' @note For an example, see the paragraph on \code{\link{reproFitTT}}.
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#' @importFrom reshape2 melt
#'
#' @keywords plot
#' 
#' @export
plot.reproFitTT <- function(x,
                            xlab = "Concentration",
                            ylab = "Nb of offspring per ind.day",
                            main = NULL,
                            fitcol = "red",
                            fitlty = 1,
                            fitlwd = 1,
                            spaghetti = FALSE,
                            cicol = "pink1",
                            cilty = 1,
                            cilwd = 1,
                            adddata = FALSE,
                            addlegend = FALSE,
                            log.scale = FALSE,
                            style = "generic", ...) {
  # plot the fitted curve estimated by reproFitTT
  # INPUTS
  # - x:  reproFitTT object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - cicol : color ci
  # - cilty : type line ci
  # - cilwd : width line ci
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic ou ggplot
  # OUTPUT:
  # - plot of fitted regression
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if (log.scale) x$dataTT$conc > 0 else TRUE
  
  dataTT <- x$dataTT[sel, ]
  dataTT$resp <- dataTT$Nreprocumul / dataTT$Nindtime
  transf_data_conc <- optLogTransform(log.scale, dataTT$conc)
  
  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, dataTT$conc)
    s <- seq(min(x),max(x), length = 100)
    if(log.scale) exp(s) else s
  })()
  
  # Possibly log transformed concentration values for display
  curv_conc <- optLogTransform(log.scale, display.conc)
  
  cred.int <- reproMeanCredInt(x, display.conc)
  
  spaghetti.CI <- if (spaghetti) { reproSpaghetti(x, display.conc) } else NULL
  dataCIm <- if (spaghetti) { melt(cbind(curv_conc, spaghetti.CI),
                                   id.vars = c("curv_conc", "conc"))} else NULL
  
  curv_resp <- data.frame(conc = curv_conc, resp = cred.int[["q50"]],
                          Line = "loglogistic")
  
  if (style == "generic") {
    reproFitPlotGenericCredInt(x, dataTT$conc, transf_data_conc, dataTT$resp,
                               curv_conc, curv_resp,
                               cred.int, spaghetti.CI, dataCIm,
                               xlab, ylab, fitcol, fitlty, fitlwd,
                               main, addlegend, adddata,
                               cicol, cilty, cilwd, log.scale)
  }
  else if (style == "ggplot") {
    reproFitPlotGG(x, dataTT$conc, transf_data_conc, dataTT$resp,
                   curv_conc, curv_resp,
                   cred.int, spaghetti.CI, dataCIm,
                   xlab, ylab, fitcol, fitlty, fitlwd,
                   main, addlegend, adddata,
                   cicol, cilty, cilwd, log.scale)
  }
  else stop("Unknown style")
}

#' @importFrom stats quantile rgamma
reproMeanCredInt <- function(fit, x) {
  # create the parameters for credible interval for the log logistic model
  # INPUT:
  # - fit : object of class reproFitTT
  # - x : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit
  
  mctot <- do.call("rbind", fit$mcmc)
  k <- nrow(mctot)
  # parameters
  d2 <- mctot[, "d"]
  log10b2 <- mctot[, "log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"]
  e2 <- 10^log10e2
  
  # quantiles
  qinf95 = NULL
  q50 = NULL
  qsup95 = NULL
  
  # poisson
  if (fit$model.label == "P") {
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      # IC 95%
      qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
      q50[i] <- quantile(theomean, probs = 0.5, na.rm = TRUE)
    }
  }
  
  # gamma poisson
  else if (fit$model.label == "GP") {
    # parameters
    log10omega2 <- mctot[, "log10omega"]
    omega2 <- 10^(log10omega2)
    
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      theo <- rgamma(n = k, shape = theomean / omega2, rate = 1 / omega2)
      # IC 95%
      qinf95[i] <- quantile(theo, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theo, probs = 0.975, na.rm = TRUE)
      q50[i] <- quantile(theo, probs = 0.5, na.rm = TRUE)
    }
  }
  # values for cred.int
  ci <- data.frame(qinf95 = qinf95,
                   q50 = q50,
                   qsup95 = qsup95)
  
  return(ci)
}

reproSpaghetti <- function(fit, x) {
  mctot <- do.call("rbind", fit$mcmc)
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 10)]
  k <- nrow(mctot[sel,])
  # parameters
  d2 <- mctot[, "d"][sel]
  log10b2 <- mctot[, "log10b"][sel]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"][sel]
  e2 <- 10^log10e2
  if (fit$model.label == "GP") {
    log10omega2 <- mctot[, "log10omega"][sel]
    omega2 <- 10^(log10omega2)
  }
  
  # all theorical
  dtheo <- array(data = NA, dim = c(length(x), length(e2)))
  if (fit$model.label == "GP") dtheotemp <- dtheo
  for (i in 1:length(e2)) {
    if (fit$model.label == "P") {
      dtheo[, i] <- d2[i] / (1 + (x / e2[i])^(b2[i])) # mean curve
    }
    else if (fit$model.label == "GP") {
      dtheotemp[, i] <- d2[i] / (1 + (x / e2[i])^(b2[i])) # mean curve
      dtheo[, i] <- rgamma(n = length(x), shape = dtheotemp[, i] / omega2[i], rate = 1 / omega2[i])
    }
  }
  dtheof <- as.data.frame(cbind(x, dtheo))
  names(dtheof) <- c("conc", paste0("X", 1:length(sel)))
  
  return(dtheof)
}

reproFitPlotGenericCredInt <- function(x, data_conc, transf_data_conc, data_resp,
                                       curv_conc, curv_resp,
                                       cred.int, spaghetti.CI, dataCIm,
                                       xlab, ylab, fitcol, fitlty, fitlwd,
                                       main, addlegend, adddata,
                                       cicol, cilty, cilwd, log.scale) {
  # plot the fitted curve estimated by reproFitTT
  # with generic style with credible interval
  
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(c(data_resp, cred.int[["qsup95"]])) + 0.01),
       type = "n")
  
  # axis
  axis(side = 2, at = pretty(c(0, max(cred.int[["qsup95"]]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  # cred.int ribbon + lines
  if (!is.null(spaghetti.CI)) {
    color <- "gray"
    color_transparent <- adjustcolor(color, alpha.f = 0.05)
    by(dataCIm, dataCIm$variable, function(x) {
      lines(x$curv_conc, x$value, col = color_transparent)
    })
  } else {
    polygon(c(curv_conc, rev(curv_conc)), c(cred.int[["qinf95"]],
                                            rev(cred.int[["qsup95"]])),
            col = cicol, border = NA)
  }
  
  lines(curv_conc, cred.int[["qsup95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, cred.int[["qinf95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  
  # fitted curve
  lines(curv_conc, curv_resp[, "resp"], col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  if (adddata) {
    par(new = TRUE)
    plotDoseResponse.reproData(x = x$transformed.data,
                               ylim = c(0, max(c(data_resp,
                                                 cred.int[["qsup95"]])) + 0.01),
                               target.time = unique(x$dataTT$time),
                               style = "generic",
                               log.scale = log.scale, addlegend = FALSE,
                               remove.someLabels = FALSE,
                               axis = FALSE)
  }
  
  # legend
  if(addlegend)  {
    legend("bottomleft",
           pch = c(ifelse(adddata, 16, NA), NA, NA, NA),
           lty = c(NA, ifelse(adddata, 1, NA), cilty, fitlty),
           lwd = c(NA, ifelse(adddata, 1, NA), cilwd, fitlwd),
           col = c(ifelse(adddata, 1, NA), 1, fitcol, cicol),
           legend = c(ifelse(adddata, "Observed values", NA),
                      ifelse(adddata, "Confidence interval", NA),
                      "Credible limits", "loglogistic"),
           bty = "n")
  }
}

reproFitPlotGGCredInt <- function(curv_resp, cred.int, spaghetti.CI, dataCIm,
                                  cicol, cilty, cilwd, valCols, fitlty, fitlwd,
                                  xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = curv_resp$conc,
                           qinf95 = cred.int[["qinf95"]],
                           qsup95 = cred.int[["qsup95"]],
                           Cred.Lim = "Credible limits")

  plt_31 <- if (!is.null(spaghetti.CI)) {
    ggplot(data.three) + geom_line(data = dataCIm, aes(x = curv_conc, y = value,
                                                       group = variable),
                                   col = "gray", alpha = 0.05)
  } else {
    ggplot(data.three) + geom_ribbon(data = data.three, aes(x = conc, ymin = qinf95,
                                                            ymax = qsup95),
                                     fill = valCols$cols3, col = valCols$cols3, alpha = 0.4)
  }
  
  plt_3 <- plt_31 +
    geom_line(data = data.three, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.three, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    scale_color_discrete(name = "") +
    theme_minimal()
  
  # plot IC
  # final plot
  
  if (!is.null(spaghetti.CI)) {
    plt_40 <- ggplot(data.three) +
      geom_line(data = dataCIm, aes(x = curv_conc, y = value, group = variable),
                col = "gray", alpha = 0.05)
  } else {
    plt_40 <- ggplot(data.three) + geom_ribbon(data = data.three, aes(x = conc, ymin = qinf95,
                                                         ymax = qsup95),
                                  fill = valCols$cols3,
                                  col = valCols$cols3, alpha = 0.4)
  }

  plt_4 <- plt_40 +
    geom_line(data = data.three, aes(conc, qinf95),
              linetype = cilty, size = cilwd, color = valCols$cols3) +
    geom_line(data = data.three, aes(conc, qsup95),
              linetype = cilty, size = cilwd, color = valCols$cols3) +
    geom_line(aes(conc, resp), curv_resp,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    ylim(0, max(cred.int[["qsup95"]]) + 0.2) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

reproFitPlotGG <- function(x, data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           cred.int, spaghetti.CI, dataCIm,
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend, adddata,
                           cicol, cilty, cilwd, log.scale) {
  
  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }
  
  # dataframes points (data) and curve (curv)
  # colors
  valCols <- fCols(curv_resp, fitcol, cicol, "repro")
  
  if (adddata) {
    plt_1 <- plotDoseResponse.reproData(x = x$transformed.data,
                                        target.time = unique(x$dataTT$time),
                                        style = "ggplot",
                                        log.scale = log.scale, addlegend = TRUE)
    mylegend_1 <- legendGgplotFit(plt_1) # mean line legend
  }
  
  plt_4 <-
    reproFitPlotGGCredInt(curv_resp, cred.int, spaghetti.CI, dataCIm,
                          cicol, cilty, cilwd, valCols, fitlty, fitlwd, xlab,
                          ylab, main)$plt_4
  
  if (adddata) {
    plt_4 <- plt_4 +
      geom_segment(aes(x = jittered_conc, xend = jittered_conc,
                       y = reproRateInf, yend = reproRateSup),
                   data = plt_1$data,
                   arrow = arrow(length = unit(0.1, "cm"),
                                 angle = 90, ends = "both")) +
      geom_point(aes(x = jittered_conc, y = resp), plt_1$data)
  }
  
  if (addlegend) {
    
    # create legends
    
    # curve (to create the legend)
    plt_2 <- ggplot(curv_resp) +
      geom_line(data = curv_resp, aes(conc, resp, colour = Line),
                linetype = fitlty, size = fitlwd) +
      scale_color_manual("", values = valCols$cols2) +
      theme_minimal()
    
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
    
    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)
    
    plt_3 <- reproFitPlotGGCredInt(curv_resp, cred.int, spaghetti.CI, dataCIm,
                                   cicol, cilty, cilwd, valCols, fitlty,
                                   fitlwd, xlab, ylab, main)$plt_3
    
    mylegend_3 <- legendGgplotFit(plt_3)
    
    if (adddata) {
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, mylegend_3,
                                      nrow = 6), ncol = 2,
                   widths = c(6, 2))
    } else {
      grid.arrange(plt_5, arrangeGrob(mylegend_2, mylegend_3,
                                      nrow = 6), ncol = 2,
                   widths = c(6, 2))
    }
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = transf_data_conc,
                                        labels = data_conc)
    return(plt_5)
  }
}

