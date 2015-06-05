reproEvalFit <- function(fit, x) {
  # eval the fitted function on x
  # INPUT :
  # - fit: repro.fit object
  # - x: vector of concentrations
  # OUTPUT :
  # - fNsurvtheo
  
  # unlog parameters
  d <- res.M$quantiles["d", "50%"]
  b <- 10^res.M$quantiles["log10b", "50%"]
  e <- 10^res.M$quantiles["log10e", "50%"]
  fNcumulpidtheo <- d / (1 + ( X / e)^b) # mean curve equation
  
  return(fNcumulpidtheo)
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
                            type = "generic",
                            ppc = FALSE, ...) {
  # plot the fitted curve estimated by repro.tt.fit
  # INPUTS
  # - x:  repro.tt.fit object
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
  # - type : generic ou ggplot
  # - ppc : plot posterior predictive check
  # OUTPUT:
  # - plot of fitted regression
  
  if (type == "ggplot") {
    # variables declarations
    concentrations.sel2. = NULL
    response.sel2. = NULL
    X.sel. = NULL
    fNcumulpidtheo.sel. = NULL
    qinf95.sel. = NULL
    qsup95.sel. = NULL
    CI.qinf95.sel. = NULL
    CI.qsup95.sel. = NULL
    Line = NULL
    Ci = NULL
  }
  
  # Define data
  concentrations <- x$dataTt$conc
  response <- x$dataTt$Nreprocumul / x$dataTt$Nindtime
  if (ppc) {
    Nreprocumul <- x$dataTt$Nreprocumul
    Nindtime <- x$dataTt$Nindtime
    res.Mtot <- do.call("rbind", x$mcmc)
  }
  
  # Fitted curve parameters
  X <- logTransXaxisFit(log.scale, concentrations) # X log.scale transformation
  
  res.M <- summary(x$mcmc)
  
  # Choose median of posteriors as parameter values
  fNcumulpidtheo <- reproLlmFit(res.M, X)
  
  # calculate IC 95 values
  if (ci) {
    CI <- reproLlmCi(x, X)
  }
  
  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  temp.conc.lt <- logTransConcFit(log.scale, X, x, concentrations)
  sel <- temp.conc.lt$sel
  X <- temp.conc.lt$X
  sel2 <- temp.conc.lt$sel2
  concentrations <- temp.conc.lt$concentrations
  
  rm(temp.conc.lt)
  
  nomortality <- match(x$dataTt$Nsurv[sel2] == x$dataTt$Ninit[sel2],
                       c(TRUE, FALSE)) # valid if at least one replicat
  
  # without mortality
  mortality <- mortality[nomortality] # vector of 0 and 1
  
  # encodes mortality empty dots (1) and not mortality solid dots (19)
  if (type == "generic") {
    mortality[which(mortality == 0)] <- 19
  }
  if (type == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }
  
  # default axis parameters
  if (missing(xlab)) {
    xlab <- "Concentrations"
  }
  if(missing(ylab)) {
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
  if (type == "generic") {
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
  if(missing(cilwd)) {
    cilwd <- 1
  }
  
  # Plotting data
  if (type == "generic") {
    if (!ci) { # CI no
      plot(concentrations[sel2], response[sel2],
           xlab = xlab,
           ylab = ylab,
           main = main,
           pch = mortality,
           xaxt = "n",
           yaxt = "n",
           ylim = c(0, max(response)),
           ...)
      # axis
      axis(side = 2, at = pretty(c(0, max(response))))
      axis(side = 1, at = unique(concentrations[sel2]),
           labels = unique(x$dataTt$conc[sel2]))
      
      # fitted curve
      lines(X[sel], fNcumulpidtheo[sel], col = fitcol,
            lty = fitlty, lwd = fitlwd, type = "l")
    }
    if (ci) {# CI yes
      # plotting data
      plot(concentrations[sel2], response[sel2],
           xlab = xlab,
           ylab = ylab,
           main = main,
           xaxt = "n",
           yaxt = "n",
           ylim = c(0, max(CI$qsup95)),
           type = "n",
           ...)
      # axis
      axis(side = 2, at = pretty(c(0, max(CI$qsup95))))
      axis(side = 1, at = unique(concentrations[sel2]),
           labels = unique(x$dataTt$conc[sel2]))
      
      # Plotting the theoretical curve
      # CI ribbon + lines
      polygon(c(X[sel], rev(X[sel])), c(CI$qinf95[sel], rev(CI$qsup95[sel])),
              col = "grey40")
      lines(X[sel], CI$qsup95[sel], type = "l", col = cicol, lty = cilty,
            lwd = cilwd)
      lines(X[sel], CI$qinf95[sel], type = "l", col = cicol, lty = cilty,
            lwd = cilwd)
      
      # fitted curve
      lines(X[sel], fNcumulpidtheo[sel], col = fitcol,
            lty = fitlty, lwd = fitlwd, type = "l")
      
      # points
      points(concentrations[sel2], response[sel2], pch = mortality)
    }
    
    # legend
    if (addlegend && !ci) { # legend yes CI no
      legend(legend.position, title = legend.title, pch = c(19, 1, NA),
             lty = c(0, 0, fitlty),
             lwd = c(1, 1, fitlwd),
             col = c(1, 1, fitcol),
             legend = c(legend.name.no, legend.name.yes, x$det.part),
             bty = "n")
    }
    if (addlegend && ci) { # legend yes CI yes
      legend(legend.position, title = legend.title, pch = c(19, 1, NA, NA),
             lty = c(0, 0, fitlty, cilty),
             lwd = c(1, 1, fitlwd, cilwd),
             col = c(1, 1, fitcol, cicol),
             legend = c(legend.name.no, legend.name.yes,
                        x$det.part, paste("Credible limits of", x$det.part,
                                          sep = " ")),
             bty = "n")
    }
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
                             labels =  unique(x$dataTt$conc[sel2]))
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
                             labels = unique(x$dataTt$conc[sel2]))
      } else { # log.scale no
        plt_5 <- plt_4 + scale_x_continuous(breaks = unique(data.one$concentrations.sel2.))
      }
      return(plt_5)
    }
  }
}
