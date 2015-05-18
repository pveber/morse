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

survFitPlotGenericNoCi <- function(concentrations, response, x, X,
                                   fNsurvtheo, sel2, sel, mortality, addlegend,
                                   ...)
{
opt_args <- list(...)
xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Concentration"
ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Survival Rate"
main <- if("main" %in% names(opt_args)) opt_args[["main"]] else "NULL"
fitcol <- if("fitcol" %in% names(opt_args)) opt_args[["fitcol"]] else "red"
fitlty <- if("fitlty" %in% names(opt_args)) opt_args[["fitlty"]] else 1
fitlwd <- if("fitlwd" %in% names(opt_args)) opt_args[["fitlwd"]] else 1

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
                     fNsurvtheo, CI, sel2, sel, mortality, addlegend,
                     ...)
{
opt_args <- list(...)
xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Concentration"
ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Survival Rate"
main <- if("main" %in% names(opt_args)) opt_args[["main"]] else "NULL"
fitcol <- if("fitcol" %in% names(opt_args)) opt_args[["fitcol"]] else "red"
fitlty <- if("fitlty" %in% names(opt_args)) opt_args[["fitlty"]] else 1
fitlwd <- if("fitlwd" %in% names(opt_args)) opt_args[["fitlwd"]] else 1
cicol <- if("cicol" %in% names(opt_args)) opt_args[["cicol"]] else "red"
cilty <- if("fitlwd" %in% names(opt_args)) opt_args[["cilty"]] else 2
cilwd <- if("cilwd" %in% names(opt_args)) opt_args[["cilwd"]] else 1

plot(concentrations[sel2], response[sel2],
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
axis(side = 1, at = unique(concentrations[sel2]),
     labels = unique(x$dataTT$conc[sel2]))

# Plotting the theoretical curve
# CI ribbon + lines
polygon(c(X[sel], rev(X[sel])), c(CI$qinf95[sel], rev(CI$qsup95[sel])),
        col = "grey40", border = NA)
lines(X[sel], CI$qsup95[sel], type = "l", col = cicol, lty = cilty,
      lwd = cilwd)
lines(X[sel], CI$qinf95[sel], type = "l", col = cicol, lty = cilty,
      lwd = cilwd)

# fitted curve

lines(X[sel], fNsurvtheo[sel], col = fitcol,
      lty = fitlty, lwd = fitlwd, type = "l")

# points
points(concentrations[sel2], response[sel2], pch = mortality)

# legend
if (addlegend) {
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
  # plot the fitted curve estimated by survFitTt
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
    CI <- survLlbinomCi(x, X)
  }
  
  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  temp.conc.lt <- logTransConcFit(log.scale, X, x, concentrations)
  sel <- temp.conc.lt$sel
  X <- temp.conc.lt$X
  sel2 <- temp.conc.lt$sel2
  concentrations <- temp.conc.lt$concentrations
  
  rm(temp.conc.lt)
  
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
#   if (missing(xlab)) {
#     xlab <- "Concentrations"
#   }
#   if (missing(ylab)) {
#     ylab <- "Response"
#   }
#   
#   # default legend parameters
#   if (missing(fitcol)) {
#     fitcol <- "red"
#   }
#   if (missing(fitlty)) {
#     fitlty <- 1
#   }
#   if (missing(fitlwd)) {
#     fitlwd <- 1
#   }
#   if (missing(main)) {
#     main = NULL
#   }
  if (style == "generic") {
    legend.position <- "bottomleft"
    legend.position.ci <- "left"
  }
  legend.title <- "Mortality"
  legend.name.no <- "No"
  legend.name.yes <- "Yes"
  
#   # IC parameters
#   if (missing(cicol)) {
#     cicol <- "red"
#   }
#   if (missing(cilty)) {
#     cilty <- 2
#   }
#   if (missing(cilwd)) {
#     cilwd <- 1
#   }
  
  # Plotting data
  if (style == "generic") {
    if (!ci) {
      survFitPlotGenericNoCi(concentrations, response, x, X,
                                         fNsurvtheo, sel2, sel, mortality, ...)
      }
    if (ci) {
      survFitPlotGenericCi(concentrations, response, x, X,
                             fNsurvtheo, CI, sel2, sel, mortality, ...)
    }
  }
  
  if (style == "ggplot") {
    if (Sys.getenv("RSTUDIO") == "") {
      dev.new() # create a new page plot
      # when not use RStudio
    }
    
    # dataframes points (one) and curve (two)
    data.one <- data.frame(concentrations[sel2], response[sel2], mortality)
    data.two <- data.frame(X[sel], fNsurvtheo[sel], Line = x$det.part)
    
    # colors
    # points vector
    n <- length(unique(data.one$mortality))
    cols <- hcl(h=seq(15, 375 - 360 / n, length = n) %% 360, c = 100, l = 65)
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
      geom_line(data = data.two, aes(X.sel., fNsurvtheo.sel., color = Line),
                linetype = fitlty, size = fitlwd) +
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
                  linetype = cilty, size = cilwd) +
        scale_color_manual(values = cols3)
      
      # plot IC
      # final plot
      plt_4 <- ggplot(data.one) +
        geom_point(data = data.one, aes(concentrations.sel2., response.sel2.,
                                        color = mortality)) +
        geom_line(aes(X.sel., fNsurvtheo.sel.), data.two, linetype = fitlty,
                  size = fitlwd, color = cols2) +
        geom_line(aes(X.sel., CI.qinf95.sel.), data.three, linetype = cilty,
                  size = cilwd, color = cols3) +
        geom_line(aes(X.sel., CI.qsup95.sel.), data.three, linetype = cilty,
                  size = cilwd, color = cols3) +
        geom_ribbon(data = data.three, aes(x = X.sel., ymin = CI.qinf95.sel.,
                                           ymax = CI.qsup95.sel.),
                    alpha = 0.4) +
        scale_color_discrete(guide = "none") +
        ylim(0, 1) +
        labs(x = xlab, y = ylab)
      
      if (!is.null(main)) { # main title
        plt_4 <- plt_4 + ggtitle(main)
      }
    }
    if (!ci) { # IC no
      plt_4 <- ggplot(data.one) +
        geom_point(data=data.one, aes(concentrations.sel2., response.sel2.,
                                      color = mortality)) +
        geom_line(aes(X.sel., fNsurvtheo.sel.), data.two,
                  linetype = fitlty, size = fitlwd, color = cols2) +
        scale_color_discrete(guide = "none") +
        ylim(0, 1) +
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
