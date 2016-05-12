#' Plotting method for \code{survFitTKTD} objects
#' 
#' This is the generic \code{plot} S3 method for the \code{survFitTKTD}.
#' It plots time-exposure-response fits for each concentration of
#' survival analysis.
#' 
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration;
#' When \code{adddata = TRUE} the black dots depict the \strong{observed survival
#' rate} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95 \% credible intervals for the estimated survival
#' rate (by default the red area around the fitted curve) and 95 \% confidence
#' intervals for the observed survival rate (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted ines limiting the credible band, and a spaghetti plot is added to this
#' band.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (2 \% of the MCMC chains are randomly
#' taken for this sample).
#' 
#' @param x An object of class \code{survFitTKTD}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param concentration A numeric value corresponding to some concentration in
#' \code{data}. If \code{concentration = NULL}, draws a plot for each concentration.
#' @param spaghetti if \code{TRUE}, the credible interval is represented by 
#' multiple curves
#' @param one.plot if \code{TRUE}, draws all the estimeted curves in one plot.
#' @param adddata if \code{TRUE}, adds the observed data with confidence interval
#' to the plot
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @keywords plot 
#' 
#' @examples
#' 
#' # (1) Load the survival data
#' data(propiconazole)
#' 
#' # (2) Create an object of class "survData"
#' dat <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#'
#' # (4) Plot the fitted curves in one plot 
#' plot(out)
#'
#' # (5) Plot one fitted curve by concentration with credible limits as
#' # spaghetti, data and confidence intervals
#' # and with a ggplot style
#' plot(out, spaghetti = TRUE , adddata = TRUE, one.plot = FALSE,
#'      style = "ggplot")
#'
#' # (6) Plt fitted curve for one specific concentration
#' plot(out, concentration = 36, style = "ggplot")
#' }
#' 
#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTKTD <- function(x,
                             xlab = "Time",
                             ylab = "Survival rate",
                             main = NULL,
                             concentration = NULL,
                             spaghetti = FALSE,
                             one.plot = FALSE,
                             adddata = FALSE,
                             addlegend = FALSE,
                             style = "generic", ...) {
  
  if (one.plot && !is.null(concentration))
    one.plot <- FALSE
  
  if ((addlegend && is.null(concentration)) ||
      (addlegend && !one.plot))
    warning("The legend is available only if [one.plot] is TRUE or if [concentration] is not NULL !", call. = FALSE)
  
  if (!is.null(concentration) && !any(x$transformed.data$conc == concentration))
    stop("The [concentration] argument is not one of the possible concentration !")
  
  if (one.plot)
    warning("The credible limits and confidence intervals are not drawn when 'one.plot' = TRUE.", call. = FALSE)
  
  conf.int <- survTKTDConfInt(x)
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     time = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init,
                     Points = "Observed values",
                     color = as.numeric(as.factor(x$transformed.data$conc)),
                     conf.int)
  
  
  # remove time 0 in dobs
  dobs <- filter(dobs, time != 0)
  
  data.credInt <- survFitPlotCITKTD(x)
  
  if (style == "generic") {
    survFitPlotTKTDGeneric(data.credInt, dobs, xlab, ylab, main, concentration,
                           one.plot, spaghetti,
                           adddata, addlegend)
  }
  else if (style == "ggplot") {
    survFitPlotTKTDGG(data.credInt, dobs, xlab, ylab, main, concentration,
                      one.plot, spaghetti, adddata, addlegend)
  }
  else stop("Unknown style")
}

Surv <- function(Cw, time, ks, kd, NEC, m0)
{
  S <- exp(-m0*time)
  x <- ifelse(Cw > NEC, 1 - NEC/Cw, NA)
  tNEC <- -(1/kd)*log(x)
  y <- ifelse(time > tNEC,
              exp( ks/kd*Cw*(exp(-kd*tNEC) -exp(-kd*time))
                   - ks*(Cw-NEC)*(time - tNEC)),
              NA)
  return(ifelse(!is.na(x) & !is.na(y), S * y, S))
}

#' @importFrom stats aggregate binom.test
survTKTDConfInt <- function(x) {
  # create confidente interval on observed data
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # OUTPUT:
  # - ci : confidente interval
  
  ci <- apply(x$transformed.data, 1, function(x) {
    binom.test(x["N_alive"], x["N_init"])$conf.int
  })
  ci <- as.data.frame(t(ci))
  colnames(ci) <- c("qinf95", "qsup95")
  ci$Conf.Int <- "Confidence interval"
  
  return(ci)
}

survFitPlotCITKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  ks <- 10^mctot[, "log10ks"]
  kd <- 10^mctot[, "log10kd"]
  m0 <- 10^mctot[, "log10m0"]
  nec <- 10^mctot[, "log10NEC"]
  
  # all theorical
  
  k <- 1:length(concobs)
  j <- 1:npoints
  
  dtheo <- lapply(k, function(y) { # conc
    sapply(j, function(z) { # time
      Surv(Cw = concobs[y], time = tfin[z],
           ks = ks, kd = kd,
           NEC = nec,
           m0 = m0)
    })
  })
  
  # transpose dtheo
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  # quantile
  qinf95 <- apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE)
  qsup95 <- apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  q50 <- apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE)
  
  
  dtheo <- as.data.frame(cbind(rep(concobs, rep(npoints, length(concobs))),
                               rep(tfin, length(concobs)),
                               dtheo))

  names(dtheo) <- c("conc", "time", paste0("X", 1:length(nec)))
  
  # divide number of mcmc by 50
  sel <- sample(ncol(dtheo[3:ncol(dtheo)]))[1:ceiling(ncol(dtheo) / 50)]
  dtheo <- cbind(dtheo[, 1:2], dtheo[, sel])
  
  # add credible limits
  dtheo$qinf95 <- qinf95
  dtheo$qsup95 <- qsup95
  dtheo$q50 <- q50
  
  # names for legend plots
  dtheo$Cred.Lim <- "Credible limits"
  dtheo$Mean.C <- "Mean curve"
  # vector color
  dtheo$color <- as.numeric(as.factor(dtheo$conc))
  
  return(dtheo)
}

survFitPlotTKTDGeneric <- function(data, dobs, xlab, ylab, main, concentration,
                                   one.plot, spaghetti, adddata,
                                   addlegend) {
  
  if (one.plot) {
    survFitPlotTKTDGenericOnePlot(data, dobs, xlab, ylab, main, adddata,
                                  addlegend)
  } else if (!one.plot && is.null(concentration)) {
    par(mfrow = plotMatrixGeometry(length(unique(dobs$conc))))
    
    survFitPlotTKTDGenericNoOnePlot(data, dobs, xlab, ylab, spaghetti,
                                    adddata, concentration)
    
    par(mfrow = c(1, 1))
  } else {
    survFitPlotTKTDGenericNoOnePlot(data, dobs, xlab, ylab, spaghetti,
                                    adddata, concentration, addlegend)
  }
}

survFitPlotTKTDGenericOnePlot <- function(data, dobs, xlab, ylab, main, adddata,
                                          addlegend) {

  plot(c(0, dobs$time),
       c(dobs$psurv, 1),
       xlab = xlab,
       ylab = ylab,
       type = "n",
       main = main)
  
  by(data, list(data$conc),
     function(x) {
       lines(x$time, x$q50, # lines
             col = x$color)
     })
  
  # points
  if (adddata) {
    points(dobs$time,
           dobs$psurv,
           pch = 20,
           col = dobs$color)
  }
  if (addlegend) {
    legend("bottomleft",
           legend = c(ifelse(adddata, "Observed values", NA),
                      unique(dobs$conc)),
           pch = c(ifelse(adddata, 20, NA),
                   rep(NA, length(unique(dobs$conc)))),
           lty = c(NA, rep(1, length(unique(dobs$conc)))),
           bty = "n",
           cex = 1,
           ncol = 2,
           col = c(ifelse(adddata, "black", NA),
                   unique(dobs$color)),
           title = "Concentrations")
  }
}

#' importFrom reshape2 melt
survFitPlotTKTDGenericNoOnePlot <- function(data, dobs, xlab, ylab, spaghetti,
                                            adddata, concentration,
                                            addlegend = FALSE) {
  
  delta <- 0.01 * (max(dobs$time) - min(dobs$time))
  
  #remove unused column before melt
  dataTmp <- data
  dataTmp[, c("Cred.Lim", "Mean.C")] <- list(NULL)
  dtheoQm <- melt(dataTmp,
                  id.vars = c("conc", "time", "color"))
  rm(dataTmp)
  
  if (is.null(concentration)) {
    dobs <- split(dobs, dobs$conc)
    dtheoQ <- split(data, data$conc)
    dtheoQm <- split(dtheoQm, dtheoQm$conc)
    
  } else {
    dobs <- list(filter(dobs, conc == concentration))
    dtheoQ <- list(filter(data, conc == concentration))
    dtheoQm <- list(filter(dtheoQm, conc == concentration))
  }
  
  mapply(function(x, y, z) {
    plot(x[, "time"],
         x[, "q50"],
         xlab = xlab,
         ylab = ylab,
         type = "n",
         ylim = c(0, 1),
         main = paste0("Concentration = ", unique(x[, "conc"])))
    
    if (spaghetti) {
      color <- "gray"
      color_transparent <- adjustcolor(color, alpha.f = 0.05)
      by(z, z$variable, function(x) {
        lines(x[, "time"], x[, "value"], col = color_transparent)
      })
    } else {
      polygon(c(x[, "time"], rev(x[, "time"])), c(x[, "qinf95"],
                                                  rev(x[, "qsup95"])),
              col = "pink", border = NA)
    }
    
    lines(x[, "time"], x[, "q50"], # lines
          col = "red")
    lines(x[, "time"], x[, "qinf95"],
          col = "pink")
    lines(x[, "time"], x[, "qsup95"], 
          col = "pink")
    
    if (adddata) {
      points(y[, "time"],
             y[, "psurv"],
             pch = 20,
             col = "black") # points
      segments(y[, "time"], y[, "qinf95"],
               y[, "time"], y[, "qsup95"],
               col = "black")
    }
    
    if (addlegend) {
      legend("bottomleft", pch = c(ifelse(adddata, 20, NA), NA, NA, NA),
             lty = c(NA, ifelse(adddata, 1, NA), 1, 1),
             col = c(ifelse(adddata, "black", NA),
                     ifelse(adddata, "black", NA),
                     "red", "pink"),
             legend = c(ifelse(adddata, "Observed values", NA),
                        ifelse(adddata, "Confidence interval", NA),
                        "Mean curve", "Credible limits"),
             bty = "n")
    }
  }, x = dtheoQ, y = dobs, z = dtheoQm)
}

survFitPlotTKTDGG <- function(data, dobs, xlab, ylab, main, concentration,
                              one.plot, spaghetti, adddata, addlegend) {
  
  if (one.plot) {
    survFitPlotTKTDGGOnePlot(data, dobs, xlab, ylab, main, adddata, addlegend)
  } else if (!one.plot && is.null(concentration)) {
    survFitPlotTKTDGGNoOnePlot(data, dobs, xlab, ylab, main, spaghetti,
                               adddata, concentration)
  } else {
    survFitPlotTKTDGGNoOnePlot(data, dobs, xlab, ylab, main, spaghetti,
                               adddata, concentration, addlegend)
  }
}

survFitPlotTKTDGGOnePlot <- function(data, dobs, xlab, ylab, main, adddata,
                                     addlegend) {
  gf <- ggplot(dobs) +
    geom_line(aes(x = time, y = q50, colour = factor(conc)),
              data = data) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
  
  if (adddata) {
    gf <- gf + geom_point(aes(x = time, y = psurv, colour = factor(conc)),
                          data = dobs)
  }
  
  if (addlegend) {
    gf + scale_color_discrete("Concentrations")
  } else {
    gf + scale_color_discrete(guide = "none")
  }
}

#' importFrom dplyr filter
#' importFrom reshape2 melt
survFitPlotTKTDGGNoOnePlot <- function(data, dobs, xlab, ylab, main, spaghetti,
                                       adddata, concentration,
                                       addlegend = FALSE) {
  
  # colors
  valCols <- fCols(data, "red", "pink")
  dataTmp <- data
  dataTmp[, c("Cred.Lim", "Mean.C")] <- list(NULL)
  dtheoQm <- melt(dataTmp,
                  id.vars = c("conc", "time", "color"))
  rm(dataTmp)
  
  if (!is.null(concentration)) {
    dobs <- filter(dobs, conc == concentration)
    data <- filter(data, conc == concentration)
    dtheoQm <- filter(dtheoQm, conc == concentration)
    curv_resp <- data.frame(time = data$time,
                            resp = data$q50,
                            Line = "Mean curve")
  }
  
  if (spaghetti) {
    gf <- ggplot(dobs) +
      geom_line(data = dtheoQm, aes(x = time, y = value, group = variable),
                alpha = 0.05)
  } else {
    gf <- ggplot(dobs) +
      geom_ribbon(data = data, aes(x = time, ymin = qinf95,
                                   ymax = qsup95),
                  fill = valCols$cols4, col = valCols$cols4, alpha = 0.4)
  }
  
  if (!is.null(concentration)) {
    gf <- gf + geom_line(data = curv_resp, aes(x = time, y = resp,
                                               color = Line)) +
      scale_colour_hue(name = "")
  } else {
    gf <- gf + geom_line(data = data, aes(x = time, y = q50),
                         color = valCols$cols5)
  }
  gf <- gf + geom_line(data = data, aes(x = time, y = qinf95,
                                        color = Cred.Lim)) +
    geom_line(data = data, aes(x = time, y = qsup95,
                               color = Cred.Lim)) +
    facet_wrap(~conc) +
    scale_linetype(name = "") +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
  
  if (adddata) {
    # dataframes points (data) and curve (curv)
    gf <- gf +
      geom_point(aes(x = time, y = psurv, fill = Points),
                 data = dobs, col = valCols$cols1) +
      geom_segment(aes(x = time, xend = time, y = qinf95, yend = qsup95,
                       linetype = Conf.Int),
                   dobs, col = valCols$cols3,
                   size = 0.5) +
      scale_fill_hue("")
  } else {
    gf <- gf
  }
  
  if (addlegend) {
    gf
  } else {
    gf + theme(legend.position = "none")
  }
  
}
