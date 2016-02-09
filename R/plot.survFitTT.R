#' Plotting method for \code{survFitTT} objects
#'
#' This function plots exposure-response fits for target time survival
#' analysis (a.k.a. \code{survFitTT} objects).
#'
#' The fitted curve represents the \strong{estimated survival rate} after
#' the target time has passed as a function of the concentration of pollutant;
#' the black dots depict the \strong{observed survival rate} at each tested
#' concentration. Note that since our model does not take inter-replicate
#' variability into consideration, replicates are systematically pooled in this
#' plot. When \code{ci = TRUE}, the function plots both credible intervals for
#' the estimated survival rate (by default the red area around the fitted
#' curve) and confidence intervals for the observed survival rate (as black
#' error bars). Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#'
#' @param x an object of class \code{survFitTT}
#' @param xlab a title for the \eqn{x}-axis
#' @param ylab a title for the \eqn{y}-axis
#' @param main main title for the plot
#' @param fitcol color of the fitted curve
#' @param fitlty line type of the fitted curve
#' @param fitlwd width of the fitted curve
#' @param cicol color of the 95 \% confidence interval limits
#' @param cilty line type for the 95 \% confidence interval limits
#' @param cilwd width of the 95 \% confidence interval limits
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' @note For an example, see the paragraph on \code{\link{reproFitTT}}.
#'
#' @keywords plot
#'
#' @import grDevices
#' @import ggplot2
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot axis legend lines par points polygon segments
#' @importFrom stats aggregate
#' @importFrom reshape2 melt
#'
#' @export
plot.survFitTT <- function(x,
                           xlab = "Concentration",
                           ylab = "Survival rate",
                           main = NULL,
                           fitcol = "red",
                           fitlty = 1,
                           fitlwd = 1,
                           spaghetti = FALSE,
                           cicol = "pink1",
                           cilty = 1,
                           cilwd = 1,
                           addlegend = FALSE,
                           log.scale = FALSE,
                           style = "generic", ...) {
  # plot the fitted curve estimated by survFitTT
  # INPUTS
  # - x:  survFitTt object
  # - xlab : label x
  # - ylab : label y
  # - main : main title
  # - fitcol : color fitted curve
  # - fitlty : type line fitted curve
  # - fitlwd : width line fitted curve
  # - cicol : color ci ribbon
  # - cilty : type line ci ribbon
  # - cilwd : width line ci ribbon
  # - addlegend : boolean
  # - log.scale : x log option
  # - style : generic or ggplot
  # OUTPUT:
  # - plot of fitted regression
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) x$dataTT$conc > 0 else TRUE
  
  dataTT <- x$dataTT[sel, ]
  dataTT$resp <- dataTT$Nsurv / dataTT$Ninit
  # data points are systematically pooled, since our model does not
  # take individual variation into account
  dataTT <- aggregate(resp ~ conc, dataTT, mean)
  transf_data_conc <- optLogTransform(log.scale, dataTT$conc)
  
  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, dataTT$conc)
    s <- seq(min(x),max(x), length = 100)
    if(log.scale) exp(s) else s
  })()
  
  # Possibly log transformed concentration values for display
  curv_conc <- optLogTransform(log.scale, display.conc)
  
  conf.int <- survLlbinomConfInt(x, log.scale)
  cred.int <- survMeanCredInt(x, display.conc)
  spaghetti.CI <- if (spaghetti) { survSpaghetti(x, display.conc) } else NULL
  dataCIm <- if (spaghetti) {melt(cbind(curv_conc, spaghetti.CI),
                                  id.vars = c("curv_conc", "conc"))} else NULL
  
  curv_resp <- data.frame(conc = curv_conc, resp = cred.int[["q50"]],
                          Line = "loglogistic")
  
  if (style == "generic") {
    survFitPlotGenericCredInt(x,
                              dataTT$conc, transf_data_conc, dataTT$resp,
                              curv_conc, curv_resp,
                              conf.int, cred.int, spaghetti.CI, dataCIm,
                              xlab, ylab, fitcol, fitlty, fitlwd,
                              main, addlegend,
                              cicol, cilty, cilwd, log.scale)
  }
  else if (style == "ggplot") {
    survFitPlotGG(x,
                  dataTT$conc, transf_data_conc, dataTT$resp,
                  curv_conc, curv_resp,
                  conf.int, cred.int, spaghetti.CI, dataCIm,
                  xlab, ylab, fitcol, fitlty, fitlwd,
                  main, addlegend,
                  cicol, cilty, cilwd / 2)
  }
  else stop("Unknown style")
}

#' @importFrom stats aggregate binom.test
survLlbinomConfInt <- function(x, log.scale) {
  # create confidente interval on observed data for the log logistic
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # - log.scale : boolean
  # OUTPUT:
  
  # - ci : confidente interval
  x <- cbind(aggregate(Nsurv ~ time + conc, x$dataTT, sum),
             Ninit = aggregate(Ninit ~ time + conc, x$dataTT, sum)$Ninit)
  
  ci <- apply(x, 1, function(x) {
    binom.test(x["Nsurv"], x["Ninit"])$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$conc
  
  if (log.scale) ci <- ci[ ,colnames(ci) != 0]
  
  return(ci)
}

#' @importFrom stats quantile
survMeanCredInt <- function(fit, x) {
  # create the parameters for credible interval for the log logistic binomial
  # model
  # INPUT:
  # - fit : object of class survFitTT
  # - x : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit
  
  mctot <- do.call("rbind", fit$mcmc)
  k <- nrow(mctot)
  # parameters
  if (fit$det.part == "loglogisticbinom_3") {
    d2 <- mctot[, "d"]
  }
  log10b2 <- mctot[, "log10b"]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"]
  e2 <- 10^log10e2
  
  # quantiles
  qinf95 = NULL
  q50 = NULL
  qsup95 = NULL
  
  for (i in 1:length(x)) {
    # llbinom 2 parameters
    if (fit$det.part == "loglogisticbinom_2") {
      theomean <- 1 / (1 + (x[i] / e2)^(b2)) # mean curve
    }
    
    # llbinom 3 parameters
    else if (fit$det.part == "loglogisticbinom_3") {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
    }
    # IC 95%
    qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
    q50[i] <- quantile(theomean, probs = 0.5, na.rm = TRUE)
    qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
  }
  
  # values for CI
  ci <- data.frame(qinf95 = qinf95,
                   q50 = q50,
                   qsup95 = qsup95)
  
  return(ci)
}

survSpaghetti <- function(fit, x) {
  mctot <- do.call("rbind", fit$mcmc)
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 10)]
  
  # parameters
  if (fit$det.part == "loglogisticbinom_3") {
    d2 <- mctot[, "d"][sel]
  }
  log10b2 <- mctot[, "log10b"][sel]
  b2 <- 10^log10b2
  log10e2 <- mctot[, "log10e"][sel]
  e2 <- 10^log10e2
  
  # all theorical
  dtheo <- array(data = NA, dim = c(length(x), length(e2)))
  for (i in 1:length(e2)) {
    # llbinom 2 parameters
    if (fit$det.part == "loglogisticbinom_2") {
      dtheo[, i] <- 1 / (1 + (x / e2[i])^(b2[i])) # mean curve
    }
    # llbinom 3 parameters
    else if (fit$det.part == "loglogisticbinom_3") {
      dtheo[, i] <- d2[i] / (1 + (x / e2[i])^(b2[i])) # mean curve
    }
  }
  dtheof <- as.data.frame(cbind(x, dtheo))
  names(dtheof) <- c("conc", paste0("X", 1:length(sel)))
  
  return(dtheof)
}

survFitPlotGenericCredInt <- function(x,
                                      data_conc, transf_data_conc, data_resp,
                                      curv_conc, curv_resp,
                                      conf.int, cred.int, spaghetti.CI, dataCIm,
                                      xlab, ylab, fitcol, fitlty, fitlwd,
                                      main, addlegend,
                                      cicol, cilty, cilwd, log.scale)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style with credible interval
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(conf.int["qsup95",]) + 0.01),
       type = "n")
  
  # axis
  axis(side = 2, at = pretty(c(0, max(conf.int["qsup95",]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  # CI ribbon + lines
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
  
  # segment CI
  
  segments(transf_data_conc, data_resp,
           transf_data_conc, conf.int["qsup95", ])
  
  Bond <- if (log.scale) {
    0.03 * (max(transf_data_conc) - min(transf_data_conc))
  } else {
    0.03 * (max(transf_data_conc) - min(transf_data_conc[which(transf_data_conc != 0)]))
  }
  
  segments(transf_data_conc - Bond,
           conf.int["qsup95", ],
           transf_data_conc + Bond,
           conf.int["qsup95", ])
  
  segments(transf_data_conc, data_resp,
           transf_data_conc, conf.int["qinf95", ])
  
  segments(transf_data_conc - Bond,
           conf.int["qinf95", ],
           transf_data_conc + Bond,
           conf.int["qinf95", ])
  
  # points
  points(transf_data_conc, data_resp, pch = 16)
  # fitted curve
  lines(curv_conc, curv_resp[, "resp"], type = "l", col = fitcol, lty = fitlty,
        lwd = fitlwd)
  
  # legend
  if (addlegend) {
    legend("bottomleft", pch = c(16, NA, NA, NA),
           lty = c(NA, 1, cilty, fitlty),
           lwd = c(NA, 1, cilwd, fitlwd),
           col = c(1, 1, cicol, fitcol),
           legend = c("Observed values", "Confidence interval",
                      "Credible limits", x$det.part),
           bty = "n")
  }
}

#' @importFrom grid arrow unit
survFitPlotGGCredInt <- function(x, data, curv_resp, conf.int, cred.int,
                                 spaghetti.CI, dataCIm, cilty, cilwd,
                                 valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = data$transf_conc,
                           qinf95 = conf.int["qinf95",],
                           qsup95 = conf.int["qsup95",],
                           Conf.Int = "Confidence interval")
  data.four <- data.frame(conc = curv_resp$conc,
                          qinf95 = cred.int[["qinf95"]],
                          qsup95 = cred.int[["qsup95"]],
                          Cred.Lim = "Credible limits")
  
  plt_3 <- ggplot(data) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95,
                     linetype = Conf.Int),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                               ends = "both"), data.three,
                 color = valCols$cols3) +
    scale_linetype(name = "") +
    theme_minimal()
  
  plt_302 <- ggplot(data) +
    geom_line(data = data.four, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.four, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = curv_resp, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(name = "") +
    theme_minimal()
  
  plt_32 <- if (!is.null(spaghetti.CI)) {
    plt_302 + geom_line(data = dataCIm, aes(x = curv_conc, y = value,
                                            group = variable),
                        alpha = 0.05)
  } else {
    plt_302 + geom_ribbon(data = data.four, aes(x = conc, ymin = qinf95,
                                                ymax = qsup95),
                          fill = valCols$cols4, col = valCols$cols4, alpha = 0.4)
  }
  
  # plot IC
  # final plot
  plt_40 <- ggplot(data) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95,
                     color = Conf.Int),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                               ends = "both"),
                 data.three, linetype = cilty,
                 size = cilwd, col = valCols$cols3) +
    geom_line(data = data.four, aes(conc, qinf95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.four, aes(conc, qsup95, color = Cred.Lim),
              linetype = cilty, size = cilwd) +
    geom_line(data = curv_resp, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd, col = valCols$cols2) +
    geom_point(data = data, aes(transf_conc, resp)) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  if (!is.null(spaghetti.CI)) {
    plt_4 <- plt_40 +
      geom_line(data = dataCIm, aes(x = curv_conc, y = value, group = variable),
                alpha = 0.05)
  } else {
    plt_4 <- plt_40 + geom_ribbon(data = data.four, aes(x = conc, ymin = qinf95,
                                                        ymax = qsup95),
                                  fill = valCols$cols4,
                                  col = valCols$cols4, alpha = 0.4)
  }
  
  return(list(plt_3 = plt_3,
              plt_32 = plt_32,
              plt_4 = plt_4))
}

survFitPlotGG <- function(x,
                          data_conc, transf_data_conc, data_resp,
                          curv_conc, curv_resp,
                          conf.int, cred.int, spaghetti.CI, dataCIm,
                          xlab, ylab, fitcol, fitlty, fitlwd,
                          main, addlegend,
                          cicol, cilty, cilwd) {
  
  
  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }
  
  # dataframes points (data) and curve (curv)
  data <- data.frame(conc = data_conc, transf_conc = transf_data_conc,
                     resp = data_resp, Points = "Observed values")
  
  # colors
  valCols <- fCols(data, fitcol, cicol, "surv")
  
  # points (to create the legend)
  plt_1 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp, fill = Points),
               col = valCols$cols1) + scale_fill_hue("") +
    theme_minimal()
  
  # curve (to create the legend)
  plt_2 <- ggplot(data) +
    geom_line(data = curv_resp, aes(conc, resp, colour = Line),
              linetype = fitlty, size = fitlwd) +
    scale_colour_manual("", values = valCols$cols2) +
    theme_minimal()
  
  plt_4 <-
    survFitPlotGGCredInt(x, data, curv_resp, conf.int, cred.int, spaghetti.CI,
                         dataCIm, cilty, cilwd,
                         valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
  
  if (addlegend) { # legend yes
    # create legends
    mylegend_1 <- legendGgplotFit(plt_1) # points legend
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
    
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    
    plt_3 <- survFitPlotGGCredInt(x, data, curv_resp, conf.int, cred.int, 
                                  spaghetti.CI, dataCIm, cilty, cilwd,
                                  valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3
    plt_32 <- survFitPlotGGCredInt(x, data, curv_resp, conf.int, cred.int, 
                                   spaghetti.CI, dataCIm, cilty, cilwd,
                                   valCols, fitlty, fitlwd, xlab, ylab, main)$plt_32
    mylegend_3 <- legendGgplotFit(plt_3)
    mylegend_32 <- legendGgplotFit(plt_32)
    grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_3, mylegend_32,
                                    mylegend_2, nrow = 6), ncol = 2,
                 widths = c(6, 2))
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    return(plt_5)
  }
}

