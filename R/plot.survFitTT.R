survEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: survFitTT object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNsurvtheo
  
  res.M <- summary(fit$mcmc)
  
  # unlog parameters
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  
  if (fit$det.part == "loglogisticbinom_3") {
    d <- res.M$quantiles["d", "50%"]
    d / (1 + (x / e)^b) # mean curve equation 3 parameters
  } else {
    1 / (1 + (x / e)^b) # mean curve equation 2 parameters
  }
}

survLlbinomCI <- function(x, log.scale, conf.level) {
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
    binom.test(x["Nsurv"], x["Ninit"], conf.level = conf.level)$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$conc
  
  if (log.scale) ci <- ci[ ,colnames(ci) != 0]
  
  return(ci)
}


survFitPlotGenericNoCI <- function(x,
                                   data_conc, transf_data_conc, data_resp,
                                   curv_conc, curv_resp,
                                   xlab, ylab, fitcol, fitlty, fitlwd,
                                   main, addlegend, ...)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style without credible interval
  
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       pch  = 16,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, 1.05),
       ...)
  
  # axis
  axis(side = 2, at = pretty(c(0, 1)))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  # legend
  if (addlegend) {
    legend("bottomleft",
           pch = c(16, NA),
           lty = c(NA, fitlty),
           lwd = c(NA, fitlwd),
           col = c(1, fitcol),
           legend = c("Observed values", x$det.part),
           bty = "n")
  }
}

#' @importFrom plotrix plotCI
survFitPlotGenericCI <- function(x,
                                 data_conc, transf_data_conc, data_resp,
                                 curv_conc, curv_resp,
                                 CI,
                                 xlab, ylab, fitcol, fitlty, fitlwd,
                                 main, addlegend,
                                 cicol, cilty, cilwd, ...)
{
  # plot the fitted curve estimated by survFitTT
  # with generic style with credible interval
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(CI["qsup95",]) + 0.2),
       type = "n",
       ...)
  
  # axis
  axis(side = 2, at = pretty(c(0, max(CI["qsup95",]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  # points
  points(transf_data_conc, data_resp, pch = 16)
  
  # segment CI
  arrows(transf_data_conc, data_resp,
         transf_data_conc, CI["qsup95", ],
         col = cicol, lty = cilty, lwd = cilwd, angle = 90, length = 0.1)
  
  arrows(transf_data_conc, data_resp,
         transf_data_conc, CI["qinf95", ],
         col = cicol, lty = cilty, lwd = cilwd, angle = 90, length = 0.1)
  
  # legend
  if (addlegend) {
    legend("bottomleft", pch = c(16, NA, NA),
           lty = c(NA, cilty, fitlty),
           lwd = c(NA, cilwd, fitlwd),
           col = c(1, cicol, fitcol),
           legend = c("Observed values", "Confidence interval", x$det.part),
           bty = "n")
  }
}


survFitPlotGeneric <- function(x,
                               data_conc, transf_data_conc, data_resp,
                               curv_conc, curv_resp,
                               CI,
                               xlab, ylab, fitcol, fitlty, fitlwd,
                               main, addlegend,
                               cicol, cilty, cilwd, ...) {
  
  
  if(!is.null(CI)) survFitPlotGenericCI(x,
                                        data_conc, transf_data_conc, data_resp,
                                        curv_conc, curv_resp,
                                        CI,
                                        xlab, ylab, fitcol, fitlty, fitlwd,
                                        main, addlegend,
                                        cicol, cilty, cilwd, ...)
  else {
    survFitPlotGenericNoCI(x,
                           data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend, ...)
  }
}


survFitPlotGGNoCI <- function(data, curv, valCols,
                              fitlty, fitlwd, xlab, ylab, main) {
  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp)) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(plt_4)
}

#' @importFrom grid arrow
survFitPlotGGCI <- function(x, data, curv, CI, cilty, cilwd,
                            valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = data$transf_conc,
                           qinf95 = CI["qinf95",],
                           qsup95 = CI["qsup95",],
                           CI = "Confidence interval")
  
  plt_3 <- ggplot(data) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95,
                     linetype = CI),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                              ends = "both"), data.three,
                 color = valCols$cols3) +
    scale_linetype_manual(values = cilty) + theme_minimal()
  
  # plot IC
  # final plot
  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp)) +
    geom_line(aes(conc, resp), curv, linetype = fitlty,
              size = fitlwd, color = valCols$cols2) +
    geom_segment(aes(x = conc, xend = conc, y = qinf95, yend = qsup95),
                 arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                               ends = "both"),
                 data.three, color = valCols$cols3, linetype = cilty,
                 size = cilwd) +
    scale_color_discrete(guide = "none") +
    ylim(0, max(CI["qsup95",]) + 0.2) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

survFitPlotGG <- function(x,
                          data_conc, transf_data_conc, data_resp,
                          curv_conc, curv_resp,
                          CI,
                          xlab, ylab, fitcol, fitlty, fitlwd,
                          main, addlegend,
                          cicol, cilty, cilwd, ...) {
  
  
  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }
  
  # dataframes points (data) and curve (curv)
  data <- data.frame(conc = data_conc, transf_conc = transf_data_conc,
                     resp = data_resp, Points = "Observed values")
  curv <- data.frame(conc = curv_conc, resp = curv_resp, Line = x$det.part)
  
  # colors
  valCols <- fCols(data, x, fitcol, cicol)
  
  # points (to create the legend)
  plt_1 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp, color = Points)) +
    scale_color_manual(values = valCols$cols1) +
    theme_minimal()
  
  # curve (to create the legend)
  plt_2 <- ggplot(data) +
    geom_line(data = curv, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd) +
    scale_color_manual(values = valCols$cols2) + theme_minimal()
  
  plt_4 <-
    if (! is.null(CI)) {
      survFitPlotGGCI(x, data, curv, CI, cilty, cilwd,
                      valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
    } else {
      survFitPlotGGNoCI(data, curv, valCols, fitlty, fitlwd,
                        xlab, ylab, main)
    }
  
  if (addlegend) { # legend yes
    # create legends
    mylegend_1 <- legendGgplotFit(plt_1) # points legend
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
    
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    
    if (is.null(CI)) {
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, nrow = 6),
                   ncol = 2, widths = c(6, 2))
    }
    else {
      plt_3 <- survFitPlotGGCI(x, data, curv, CI, cilty, cilwd,
                               valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3
      mylegend_3 <- legendGgplotFit(plt_3)
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_3, mylegend_2,
                                      nrow = 6), ncol = 2,
                   widths = c(6, 2))
    }
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    return(plt_5)
  }
}

#' Plotting method for reproFitTT objects
#' 
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
#' @keywords plot 
#'
#' @export
#'
#' @import grDevices
#' @import ggplot2
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
plot.survFitTT <- function(x,
                           xlab,
                           ylab,
                           main,
                           fitcol,
                           fitlty,
                           fitlwd,
                           ci = FALSE,
                           conf.level,
                           cicol,
                           cilty,
                           cilwd,
                           addlegend = FALSE,
                           log.scale = FALSE,
                           style = "generic",
                           ...) {
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
  # - conf.level : binom.test conf.level
  # - cicol : color ci
  # - cilty : type line ci
  # - cilwd : width line ci
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
  
  curv_resp <- survEvalFit(x, display.conc)
  
  # default axis parameters
  if (missing(xlab)) xlab <- "Concentrations"
  if (missing(ylab)) ylab <- "Survival rate"
  
  # default legend parameters
  if (missing(fitcol)) fitcol <- "red"
  if (missing(fitlty)) fitlty <- 1
  if (missing(fitlwd)) fitlwd <- 1
  
  if (missing(main)) main = NULL
  
  # CI parameters
  if (missing(cicol)) cicol <- "black"
  if (missing(cilty)) cilty <- 2
  if (missing(cilwd)) cilwd <- 1
  if (missing(conf.level)) conf.level <- 0.95
  
  CI <- if(ci) { survLlbinomCI(x, log.scale, conf.level) } else NULL
  
  if (style == "generic") {
    survFitPlotGeneric(x,
                       dataTT$conc, transf_data_conc, dataTT$resp,
                       curv_conc, curv_resp,
                       CI,
                       xlab, ylab, fitcol, fitlty, fitlwd,
                       main, addlegend,
                       cicol, cilty, cilwd, ...)
  }
  else if (style == "ggplot") {
    survFitPlotGG(x,
                  dataTT$conc, transf_data_conc, dataTT$resp,
                  curv_conc, curv_resp,
                  CI,
                  xlab, ylab, fitcol, fitlty, fitlwd,
                  main, addlegend,
                  cicol, cilty, cilwd, ...)
  }
  else stop("Unknown style")
}
