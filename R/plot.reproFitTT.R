reproEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: reproFitTT object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNcumulpidtheo
  
  res.M <- summary(fit$mcmc)
  
  # unlog parameters
  d <- res.M$quantiles["d", "50%"]
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  fNcumulpidtheo <- d / (1 + ( x / e)^b) # mean curve equation
  
  return(fNcumulpidtheo)
}

reproLlmCI <- function(fit, x) {
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
  qsup95 = NULL
  
  # poisson
  if (fit$model.label == "P") {
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      # IC 95%
      qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
    }
  }
  
  # gamma poisson
  if (fit$model.label == "GP") {
    # parameters
    log10omega2 <- mctot[, "log10omega"]
    omega2 <- 10^(log10omega2)
    
    for (i in 1:length(x)) {
      theomean <- d2 / (1 + (x[i] / e2)^(b2)) # mean curve
      theo <- rgamma(n = k, shape = theomean / omega2, rate = 1 / omega2)
      # IC 95%
      qinf95[i] <- quantile(theo, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theo, probs = 0.975, na.rm = TRUE)
    }
  }
  # values for CI
  ci <- list(qinf95 = qinf95,
             qsup95 = qsup95)
  
  return(ci)
}

reproFitPlotGenericNoCI <- function(data_conc, transf_data_conc, data_resp,
                                    curv_conc, curv_resp, mortality,
                                    xlab, ylab, fitcol, fitlty, fitlwd,
                                    main, addlegend, ...) {
  # plot the fitted curve estimated by reproFitTT
  # with generic style without credible interval
  
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       pch = mortality,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(data_resp) + 0.2),
       ...)
  # axis
  axis(side = 2, at = pretty(c(0, max(data_resp))))
  axis(side = 1, at = transf_data_conc,
       labels = data_conc)
  
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  # legend
  if (addlegend) {
    legend("bottomleft", pch = c(19, 1, NA),
           lty = c(0, 0, fitlty),
           lwd = c(1, 1, fitlwd),
           col = c(1, 1, fitcol),
           legend = c("No mortality", "Mortality", "loglogistic"),
           bty = "n")
  }
}

reproFitPlotGenericCI <- function(data_conc, transf_data_conc, data_resp,
                                  curv_conc, curv_resp,
                                  CI, mortality,
                                  xlab, ylab, fitcol, fitlty, fitlwd,
                                  main, addlegend,
                                  cicol, cilty, cilwd, ...) {
  # plot the fitted curve estimated by reproFitTT
  # with generic style with credible interval
  
  plot(transf_data_conc, data_resp,
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n",
       ylim = c(0, max(CI[["qsup95"]]) + 0.2),
       type = "n",
       ...)
  
  # axis
  axis(side = 2, at = pretty(c(0, max(CI[["qsup95"]]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  # CI ribbon + lines
  polygon(c(curv_conc, rev(curv_conc)), c(CI[["qinf95"]], rev(CI[["qsup95"]])),
          col = "pink1", border = NA)
  lines(curv_conc, CI[["qsup95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, CI[["qinf95"]], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  
  # fitted curve
  lines(curv_conc, curv_resp, col = fitcol,
        lty = fitlty, lwd = fitlwd, type = "l")
  
  # points
  points(transf_data_conc, data_resp, pch = mortality)
  
  # legend
  if(addlegend)
    legend("bottomleft", pch = c(19, 1, NA, NA),
           lty = c(0, 0, fitlty, cilty),
           lwd = c(1, 1, fitlwd, cilwd),
           col = c(1, 1, fitcol, cicol),
           legend = c("No mortality", "Mortality",
                      "loglogistic", "Credible limits of loglogistic"),
           bty = "n")
}

reproFitPlotGeneric <- function(data_conc, transf_data_conc, data_resp,
                                curv_conc, curv_resp,
                                CI, mortality,
                                xlab, ylab, fitcol, fitlty, fitlwd,
                                main, addlegend,
                                cicol, cilty, cilwd, ...) {
  
  if(!is.null(CI)) reproFitPlotGenericCI(data_conc, transf_data_conc,
                                         data_resp,
                                         curv_conc, curv_resp,
                                         CI, mortality,
                                         xlab, ylab, fitcol, fitlty, fitlwd,
                                         main, addlegend,
                                         cicol, cilty, cilwd, ...)
  else {
    reproFitPlotGenericNoCI(data_conc, transf_data_conc, data_resp,
                            curv_conc, curv_resp, mortality,
                            xlab, ylab, fitcol, fitlty, fitlwd,
                            main, addlegend, ...)
  }
}

reproFitPlotGGNoCI <- function(data, curv, cols2,
                               fitlty, fitlwd, xlab, ylab, main) {
  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp, shape = Mortality),
               size = 3) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = cols2) +
    scale_color_discrete(guide = "none") +
    scale_shape(guide = "none") +
    ylim(0, max(data$resp) + 1) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()
  
  return(plt_4)
}

reproFitPlotGGCI <- function(data, curv, CI, cicol, cilty, cilwd,
                             cols2, fitlty, fitlwd, xlab, ylab, main) {
  # IC
  cri <- data.frame(conc = curv$conc,
                           qinf95 = CI[["qinf95"]],
                           qsup95 = CI[["qsup95"]],
                           CI = "Credible limits of loglogistic")
  
  # colors 
  cols3 <- cicol
  names(cols3) <- "Credible limits of loglogistic"
  
  plt_3 <- ggplot(data) +
    geom_line(data = cri, aes(conc, qinf95, color = CI),
              linetype = cilty, size = cilwd) +
    geom_line(data = cri, aes(conc, qsup95, color = CI),
              linetype = cilty, size = cilwd) +
    geom_ribbon(data = cri, aes(x = conc, ymin = qinf95,
                                ymax = qsup95), fill = "pink", alpha = 0.4) +
    scale_color_manual(values = cols3) + theme_minimal()

  plt_4 <- ggplot(data) +
    geom_point(data = data, aes(transf_conc, resp,
                                shape = Mortality), size = 3) +
    geom_line(aes(conc, resp), curv,
              linetype = fitlty, size = fitlwd, color = cols2) +
    geom_line(data = cri, aes(conc, qinf95),
              linetype = cilty, size = cilwd, color = cols3) +
    geom_line(data = cri, aes(conc, qsup95),
              linetype = cilty, size = cilwd, color = cols3) +
    geom_ribbon(data = cri, aes(x = conc, ymin = qinf95,
                                ymax = qsup95), fill = "pink", alpha = 0.4) +
    scale_color_discrete(guide = "none") +
    scale_shape(guide = "none") +
    ylim(0, max(CI[["qsup95"]]) + 0.2) +
    labs(x = xlab, y = ylab) +
    ggtitle(main) + theme_minimal()

  return(list(plt_3 = plt_3,
              plt_4 = plt_4))
}

reproFitPlotGG <- function(data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           CI, mortality,
                           xlab, ylab, fitcol, fitlty, fitlwd,
                           main, addlegend,
                           cicol, cilty, cilwd, ...) {
  
  if (Sys.getenv("RSTUDIO") == "") {
    dev.new() # create a new page plot
    # when not use RStudio
  }
  
  # dataframes points (data) and curve (curv)
  data <- data.frame(conc = data_conc, transf_conc = transf_data_conc,
                     resp = data_resp, Mortality = mortality)
  curv <- data.frame(conc = curv_conc, resp = curv_resp, Line = "loglogistic")
  
  # colors
  # fitted curve
  cols2 <- fitcol
  names(cols2) <- "loglogistic"
  
  # points (to create the legend)
  plt_1 <- ggplot(data) +
    geom_point(data = data, aes(conc, resp,
                                shape = Mortality), size = 3) +
    theme_minimal()
  
  # curve (to create the legend)
  plt_2 <- ggplot(data) +
    geom_line(data = curv, aes(conc, resp, color = Line),
              linetype = fitlty, size = fitlwd) +
    scale_color_manual(values = cols2) +
    theme_minimal()
  
  plt_4 <-
    if (! is.null(CI)) {
      reproFitPlotGGCI(data, curv, CI, cicol, cilty, cilwd,
                       cols2, fitlty, fitlwd, xlab, ylab, main)$plt_4
    } else {
      reproFitPlotGGNoCI(data, curv, cols2, fitlty, fitlwd,
                         xlab, ylab, main)
    }
  
  if (addlegend) {
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
      plt_3 <- reproFitPlotGGCI(data, curv, CI, cicol, cilty, cilwd,
                                cols2, fitlty, fitlwd, xlab, ylab, main)$plt_3
      mylegend_3 <- legendGgplotFit(plt_3)
      grid.arrange(plt_5, arrangeGrob(mylegend_1, mylegend_2, mylegend_3,
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
#' @param x An object of class \code{reproFitTT}.
#' @param xlab A label for the \eqn{X}-label, by default \code{Concentrations}.
#' @param ylab A label for the \eqn{Y}-label, by default \code{Response}.
#' @param main A main title for the plot.
#' @param fitcol A single color to plot the fitted curve, by default
#' \code{red}.
#' @param fitlty A single line type to plot the fitted curve, by default
#' \code{1}.
#' @param fitlwd A single numeric which controls the width of the fitted curve,
#' by default \code{1}.
#' @param ci If \code{TRUE}, the 95 \% credible limits are draw for the model.
#' @param cicol A single color to plot the 95 \% credible limits, by default
#' \code{blue}.
#' @param cilty A single line type to plot 95 \% credible limits, by default
#' \code{1}.
#' @param cilwd A single numeric which controls the width of the 95 \% credible
#' limits, by default \code{2}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param log.scale Log option for the \eqn{X}-axis.
#' @param style Graphical package method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
#' @keywords plot 
#' 
plot.reproFitTT <- function(x,
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
                            addlegend = FALSE,
                            log.scale = FALSE,
                            style = "generic",
                            ...) {
  # plot the fitted curve estimated by reproFitTT
  # INPUTS
  # - x:  reproFitTT object
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
  # - style : generic ou ggplot
  # OUTPUT:
  # - plot of fitted regression
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) x$dataTT$conc > 0 else TRUE
  
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
  
  curv_resp <- reproEvalFit(x, display.conc)
  
  # default axis parameters
  if (missing(xlab)) xlab <- "Concentrations"
  if(missing(ylab)) ylab <- "Nb. of offspring / NID"
  
  # default legend parameters	
  if (missing(fitcol)) fitcol <- "red"
  if (missing(fitlty)) fitlty <- 1
  if (missing(fitlwd)) fitlwd <- 1
  
  if (missing(main)) main = NULL
  
  # IC parameters
  if (missing(cicol)) cicol <- "red"
  if (missing(cilty)) cilty <- 2
  if(missing(cilwd)) cilwd <- 1
  
  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  # valid if at least one replicat
  mortality <- mortality[match(dataTT$Nsurv == dataTT$Ninit,
                               c(TRUE, FALSE))] # vector of 0 and 1
  
  # encodes mortality empty dots (1) and not mortality solid dots (19)
  if (style == "generic")  mortality[which(mortality == 0)] <- 19
  else if (style == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }
  
  CI <- if (ci) { CI <- reproLlmCI(x, display.conc) } else NULL
  
  if (style == "generic") {
    reproFitPlotGeneric(dataTT$conc, transf_data_conc, dataTT$resp,
                        curv_conc, curv_resp,
                        CI, mortality,
                        xlab, ylab, fitcol, fitlty, fitlwd,
                        main, addlegend,
                        cicol, cilty, cilwd, ...)
  }
  else if (style == "ggplot") {
    reproFitPlotGG(dataTT$conc, transf_data_conc, dataTT$resp,
                   curv_conc, curv_resp,
                   CI, mortality,
                   xlab, ylab, fitcol, fitlty, fitlwd,
                   main, addlegend,
                   cicol, cilty, cilwd, ...)
  }
  else stop("Unknown style")
}