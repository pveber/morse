survEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - x: repro.fit object
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

survLlbinomCI <- function(x, X) {
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
  med = NULL
  qsup95 = NULL

  for (i in 1:length(X)) {
    if (x$det.part == "loglogisticbinom_3") {
      theomean <- d2/(1 + (X[i] / e2)^(b2))
    } else {
      theomean <- 1/(1 + (X[i] / e2)^(b2))
    }

    # IC 95%
    qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
    med[i] <- median(theomean, na.rm = TRUE)
    qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
  }
  # values for CI
  ci <- list(qinf95 = qinf95,
             med = med,
             qsup95 = qsup95)
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
  axis(side = 2, at = pretty(c(0,1)))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)

  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")

  # legend
  if (addlegend) {
    legend("bottomleft",
           lty = c(fitlty),
           lwd = c(fitlwd),
           col = c(fitcol),
           legend = c(x$det.part),
           bty = "n")
  }
}

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
       ylim = c(0, max(CI$qsup95) + 0.2),
       type = "n",
       ...)

  # axis
  axis(side = 2, at = pretty(c(0, max(CI$qsup95))))
  axis(side   = 1,
       at     = transf_data_conc,
       labels = data_conc)

  # Plotting the theoretical curve
  # CI ribbon + lines
  polygon(c(curv_conc, rev(curv_conc)), c(CI$qinf95, rev(CI$qsup95)),
          col = "grey40", border = NA)
  lines(curv_conc, CI$qsup95, type = "l", col = cicol, lty = cilty, lwd = cilwd)
  lines(curv_conc, CI$qinf95, type = "l", col = cicol, lty = cilty, lwd = cilwd)
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  # points
  points(transf_data_conc, data_resp, pch = 16)

  # legend
  if (addlegend) {
    legend("bottomleft", pch = c(NA, NA),
           lty = c(fitlty, cilty),
           lwd = c(fitlwd, cilwd),
           col = c(fitcol, cicol),
           legend = c(x$det.part, paste("Credible limits of", x$det.part,
                                        sep = " ")),
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
    geom_point(data=data, aes(transf_conc, resp)) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = valCols$cols2) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(plt_4)
}

survFitPlotGGCI <- function(x, data, curv, CI, cilty, cilwd,
                            valCols, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  data.three <- data.frame(conc = curv$conc,
                           qinf95 = CI$qinf95,
                           qsup95 = CI$qsup95,
                           CI = paste("Credible limits of", x$det.part,
                                      sep = " "))

  plt_3 <- ggplot(data) +
    geom_line(data = data.three, aes(conc, qinf95, color = CI),
              linetype = cilty, size = cilwd) +
    geom_line(data = data.three, aes(conc, qsup95, color = CI),
              linetype = cilty, size = cilwd) +
    scale_color_manual(values = valCols$cols3)

  # plot IC
  # final plot
  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(conc, resp, color = 16)) +
    geom_line(aes(conc, resp), curv, linetype = fitlty,
              size = fitlwd, color = valCols$cols2) +
    geom_line(aes(conc, qinf95), data.three, linetype = cilty,
              size = cilwd, color = valCols$cols3) +
    geom_line(aes(conc, qsup95), data.three, linetype = cilty,
              size = cilwd, color = valCols$cols3) +
    geom_ribbon(data = data.three, aes(x = conc, ymin = qinf95,
                                       ymax = qsup95),
                alpha = 0.4) +
    scale_color_discrete(guide = "none") +
    ylim(0, 1) +
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

  # dataframes points (one) and curve (two)
  data <- data.frame(conc = data_conc, transf_conc = transf_data_conc, resp = data_resp)
  curv <- data.frame(conc = curv_conc, resp = curv_resp, Line = x$det.part)

  # colors
  valCols <- fCols(data, x, fitcol, cicol)

  # points (to create the legend)
  plt_1 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp)) +
    scale_color_manual(values = valCols$cols1)

  # curve (to create the legend)
  plt_2 <- ggplot(data) +
    geom_line(data = curv, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd) +
    scale_color_manual(values = valCols$cols2)


  plt_4 <-
    if (is.null(CI))
      survFitPlotGGCI(x, data, curv, CI, cilty, cilwd,
                      valCols, fitlty, fitlwd, xlab, ylab, main)$plt_4
    else
      survFitPlotGGNoCI(data, curv, valCols, fitlty, fitlwd,
                        xlab, ylab, main)

  if (addlegend) { # legend yes
    # create legends
#    mylegend_1 <- legendGgplotFit(plt_1) # points legend
    mylegend_2 <- legendGgplotFit(plt_2) # mean line legend

    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)

    if (is.null(CI)) {
      grid.arrange(plt_5, arrangeGrob(mylegend_2, nrow = 6),
                   ncol = 2, widths = c(7,1))
    }
    else {
      plt_3 <- survFitPlotGGCI(x, data, curv, CI, cilty, cilwd,
                               valCols, fitlty, fitlwd, xlab, ylab, main)$plt_3
      mylegend_3 <- legendGgplotFit(plt_3)
      grid.arrange(plt_5, arrangeGrob(mylegend_2, mylegend_3,
                                      nrow = 6), ncol = 2, widths = c(7,1))
    }
  }
  else { # no legend
    plt_5 <- plt_4 + scale_x_continuous(breaks = data$transf_conc,
                                        labels = data$conc)
    return(plt_5)
  }
}



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
                           cicol,
                           cilty,
                           cilwd,
                           addlegend = TRUE,
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
  if (missing(ylab)) ylab <- "Response"

  # default legend parameters
  if (missing(fitcol)) fitcol <- "red"
  if (missing(fitlty)) fitlty <- 1
  if (missing(fitlwd)) fitlwd <- 1

  if (missing(main)) main = NULL

  # CI parameters
  if (missing(cicol)) cicol <- "red"
  if (missing(cilty)) cilty <- 2
  if (missing(cilwd)) cilwd <- 1
  CI <- if(ci) {
    survLlbinomCI(x, display.conc)
  }
  else NULL

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
