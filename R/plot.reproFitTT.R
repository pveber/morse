reproEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: reproFitTT object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNcumulpidtheo
  
  res.M <- summary(x$mcmc)
  
  # unlog parameters
  d <- res.M$quantiles["d", "50%"]
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  fNcumulpidtheo <- d / (1 + ( x / e)^b) # mean curve equation
  
  return(fNcumulpidtheo)
}

reproLlmCI <- function(x, X) {
  # create the parameters for credible interval for the log logistic model
  # INPUT:
  # - x : object of class reproFitTT
  # - X : vector of concentrations values (x axis)
  # OUTPUT:
  # - ci : credible limit
  
  mctot <- do.call("rbind", x$mcmc)
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
  if (x$model.label == "P") {
    for (i in 1:length(X)) {
      theomean <- d2 / (1 + (X[i] / e2)^(b2)) # mean curve
      # IC 95%
      qinf95[i] <- quantile(theomean, probs = 0.025, na.rm = TRUE)
      qsup95[i] <- quantile(theomean, probs = 0.975, na.rm = TRUE)
    }
  }
  
  # gamma poisson
  if (x$model.label == "GP") {
    # parameters
    log10omega2 <- mctot[, "log10omega"]
    omega2 <- 10^(log10omega2)
    
    for (i in 1:length(X)) {
      theomean <- d2 / (1 + (X[i] / e2)^(b2)) # mean curve
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

reproFitPlotGenericNoCI <- function(x,
                        data_conc, transf_data_conc, data_resp,
                        curv_conc, curv_resp,
                        xlab, ylab, mortality, fitcol, fitlty, fitlwd,
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
       ylim = c(0, 1.05),
       ...)
  # axis
  axis(side = 2, at = pretty(c(0, 1)))
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
           legend = c("No mortality", "Mortality", x$det.part),
           bty = "n")
  }
}

reproFitPlotGenericCI <- function(x, 
                                  data_conc, transf_data_conc, data_resp,
                                  curv_conc, curv_resp,
                                  CI,
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
       ylim = c(0, max(CI["qsup95",]) + 0.2),
       type = "n",
       ...)
  
  # axis
  axis(side = 2, at = pretty(c(0, max(CI["qsup95",]))))
  axis(side = 1,
       at = transf_data_conc,
       labels = data_conc)
  
  # Plotting the theoretical curve
  # CI ribbon + lines
  polygon(c(curv_conc, rev(curv_conc)), c(CI["qinf95",], rev(CI["qsup95",])),
          col = "grey40")
  lines(curv_conc, CI["qsup95",], type = "l", col = cicol, lty = cilty,
        lwd = cilwd)
  lines(curv_conc, CI["qinf95",], type = "l", col = cicol, lty = cilty,
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
                    x$det.part, paste("Credible limits of", x$det.part,
                                      sep = " ")),
         bty = "n")
}


reproFitPlotGeneric <- function(x,
                    data_conc, transf_data_conc, data_resp,
                    curv_conc, curv_resp,
                    CI,
                    xlab, ylab, fitcol, fitlty, fitlwd,
                    main, addlegend,
                    cicol, cilty, cilwd, ...) {
  
  if(!is.null(CI)) reproFitPlotGenericCI(x,
                                        data_conc, transf_data_conc, data_resp,
                                        curv_conc, curv_resp,
                                        CI,
                                        xlab, ylab, fitcol, fitlty, fitlwd,
                                        main, addlegend,
                                        cicol, cilty, cilwd, ...)
  else {
    reproFitPlotGenericNoCI(x,
                           data_conc, transf_data_conc, data_resp,
                           curv_conc, curv_resp,
                           xlab, ylab, mortality, fitcol, fitlty, fitlwd,
                           main, addlegend, ...)
  }
}

#' @export
#' 
#' @name reproFitTt
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.reproFitTt <- function(x,
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
  if(missing(ylab)) ylab <- "Response"
  
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
  if (type == "generic")  mortality[which(mortality == 0)] <- 19
  else if (type == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }
  
  CI <- if (ci) { CI <- reproLlmCI(x, X) } else NULL

  if (type == "generic") {
    reproFitPlotGeneric(x,
                        dataTT$conc, transf_data_conc, dataTT$resp,
                        curv_conc, curv_resp,
                        CI,
                        xlab, ylab, fitcol, fitlty, fitlwd,
                        main, addlegend,
                        cicol, cilty, cilwd, ...)
  }
  else if (style == "ggplot") {
    reproFitPlotGG(x,
                   dataTT$conc, transf_data_conc, dataTT$resp,
                   curv_conc, curv_resp,
                   CI,
                   xlab, ylab, fitcol, fitlty, fitlwd,
                   main, addlegend,
                   cicol, cilty, cilwd, ...)
  }
  else stop("Unknown style")
}

  
  # posterior predictive check
  if (ppc) {
    if (Sys.getenv("RSTUDIO") == "") {
      dev.new() # create a new page plot
      # when not use RStudio
    }
    plot(Nreprocumul, Nreprocumul,
         bty = "n",
         type = "n",
         main = "Posterior predictive check",
         ylab = "Predicted Nr of offsprings for each replicate",
         xlab = "Observed Nr of offsprings for each replicate")
    abline(a = 0, b = 1, lty = 2) # add x = y
    
    # create segments for post predictive check plot
    fNcumulpidtheo.ppc <- reproLlmFit(res.M, concentrations)
    for (i in 1:length(concentrations)) {
      nbi <- rpois(n = nrow(res.Mtot),
                   lambda = fNcumulpidtheo.ppc[i] * Nindtime[i])
      qinf95nbi <- quantile(nbi, probs = 0.025)
      qsup95nbi <- quantile(nbi, probs = 0.975)
      segments(Nreprocumul[i], qinf95nbi, Nreprocumul[i], qsup95nbi, col = "red")
    }
  }
  
  if (type == "ggplot") {
    if (Sys.getenv("RSTUDIO") == "") {
      dev.new() # create a new page plot
      # when not use RStudio
    }
    
    # dataframes points (one) and curve (two)
    data.one <- data.frame(concentrations[sel2], response[sel2], mortality)
    data.two <- data.frame(X[sel], fNcumulpidtheo[sel], Line = x$det.part)
    
    # colors
    # points vector
    n <- length(unique(data.one$mortality))
    cols <- hcl(h = seq(15, 375 - 360 / n, length = n) %% 360, c = 100, l = 65)
    cols1 <- cols[1:n]
    names(cols1) <- sort(unique(data.one$mortality))
    # fitted curve
    cols2 <- fitcol
    names(cols2) <- c(x$det.part)
    
    # points (to create the legend)
    plt_1 <- ggplot(data.one) +
      geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                      color = mortality)) +
      scale_color_manual(values = cols1)
    
    # curve (to create the legend)
    plt_2 <- ggplot(data.one) +
      geom_line(data = data.two, aes(X.sel., fNcumulpidtheo.sel.,
                                     color = Line), linetype = fitlty,
                size = fitlwd) +
      scale_color_manual(values = cols2)
    if (ci) { # IC yes
      # IC
      data.three <- data.frame(X[sel], CI$qinf95[sel], CI$qsup95[sel],
                               Ci = paste("Credible limits of", x$det.part,
                                          sep = " "))
      # colors 
      cols3 <- cicol
      names(cols3) <- c(paste("Credible limits of", x$det.part, sep = " "))
      
      plt_3 <- ggplot(data.one) +
        geom_line(data = data.three, aes(X.sel., CI.qinf95.sel., color = Ci),
                  linetype = cilty, size = cilwd) +
        geom_line(data = data.three, aes(X.sel., CI.qsup95.sel., color = Ci),
                  linetype = cilty,size = cilwd) +
        scale_color_manual(values = cols3)
      
      # plot IC
      # final plot
      plt_4 <- ggplot(data.one) +
        geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                        color = mortality)) +
        geom_line(aes(X.sel., fNcumulpidtheo.sel.), data.two,
                  linetype = fitlty, size = fitlwd, color = cols2) +
        geom_line(aes(X.sel., CI.qinf95.sel.), data.three, linetype = cilty,
                  size = cilwd, color = cols3) +
        geom_line(aes(X.sel., CI.qsup95.sel.), data.three, linetype = cilty,
                  size = cilwd, color = cols3) +
        geom_ribbon(data = data.three, aes(x = X.sel., ymin = CI.qinf95.sel.,
                                           ymax = CI.qsup95.sel.),
                    alpha = 0.4) +
        scale_color_discrete(guide = "none") +
        labs(x = xlab, y = ylab)
      
      if (!is.null(main)) { # main title
        plt_4 <- plt_4 + ggtitle(main)
      }
    }
    if (!ci) { # IC no
      plt_4 <- ggplot(data.one) +
        geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                        color = mortality)) +
        geom_line(aes(X.sel., fNcumulpidtheo.sel.), data.two,
                  linetype = fitlty, size = fitlwd, color = cols2) +
        scale_color_discrete(guide = "none") +
        labs(x = xlab, y = ylab)
      
      if (!is.null(main)) { # personal title
        plt_4 <- plt_4 + ggtitle(main)
      }
    }
    
    if (addlegend) { # legend yes
      # create legends
      mylegend_1 <- legendGgplotFit(plt_1) # points legend
      mylegend_2 <- legendGgplotFit(plt_2) # mean line legend
      
      if (ci) mylegend_3 <- legendGgplotFit(plt_3) # CI legend
      
      if (log.scale) { # log.sacle yes
        plt_5 <- plt_4 +
          scale_x_continuous(breaks = unique(data.one$concentrations.sel2.),
                             labels =  unique(x$dataTT$conc[sel2]))
      } else { # log.scale no
        plt_5 <- plt_4 +
          scale_x_continuous(breaks = unique(data.one$concentrations.sel2.))
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
