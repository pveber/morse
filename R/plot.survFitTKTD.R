#' Plotting method for survFitTKTD objects
#' 
#' @param x An object of class \code{survFitTKTD}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Concentrations}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @keywords plot 
#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTKTD <- function(x,
                             xlab = "Time",
                             ylab = "Survival rate",
                             main = NULL,
                             fitlty = 1,
                             fitlwd = 1,
                             ci = FALSE,
                             cicol = "pink1",
                             cilty = 1,
                             cilwd = 1,
                             one.plot = TRUE,
                             addlegend = FALSE,
                             style = "generic", ...) {
  
  data <- survFitPlotDataTKTD(x)
  
  dataCI <- if (ci && !one.plot) { survFitPlotCITKTD(x) } else NULL
  dataCIm <- if (ci && !one.plot) { melt(dataCI,
                                         id.vars = c("conc", "time")) } else NULL
  
  if (style == "generic") {
    survFitPlotTKTDGeneric(data, xlab, ylab, main, one.plot, ci, dataCIm)
  }
  
  if (style == "ggplot") {
    survFitPlotTKTDGG(data, xlab, ylab, main, one.plot, ci, dataCI, dataCIm)
  }
  else stop("Unknown style")
}

Surv <- function (Cw, time, ks, ke, NEC, m0)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- exp(-m0*time) # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    tNEC <- -(1/ke)*log(1 - NEC/Cw)
    if (time > tNEC) {
      # ajoute de la mortalite due au toxique
      S <- S * exp( ks/ke*Cw*(exp(-ke*tNEC) -exp(-ke*time))
                          - ks*(Cw-NEC)*(time - tNEC) )
    }
  }
  return(S)
}

survFitPlotDataTKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  dtheo <- data.frame(conc = numeric(), t = numeric(), psurv = numeric())
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # parameters
  ks <- x$estim.par["ks", "median"]
  ke <- x$estim.par["ke", "median"]
  nec <- x$estim.par["nec", "median"]
  m0 <- x$estim.par["m0", "median"]
  
  for (i in 1:length(concobs)) {
    for (j in 1:npoints) {
      psurv <- Surv(Cw = concobs[i], time = tfin[j],
                    ks = ks, ke = ke,
                    NEC = nec,
                    m0 = m0)
      dtheo <- rbind(dtheo, data.frame(conc = concobs[i],
                                       t = tfin[j],
                                       psurv = psurv))
    }
  }
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     t = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init)
  
  return(list(dtheo = dtheo,
              dobs = dobs))
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
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 10)]
  ks <- 10^mctot[, "log10ks"][sel]
  ke <- 10^mctot[, "log10ke"][sel]
  m0 <- 10^mctot[, "log10m0"][sel]
  nec <- 10^mctot[, "log10NEC"][sel]
  
  # all theorical
  dtheo = list()
  for (k in 1:length(concobs)) {
    dtheo[[k]] <- array(data = NA, dim = c(npoints, length(nec)))
    for (i in 1:length(nec)) {
      for (j in 1:npoints) {
        dtheo[[k]][j, i] <- Surv(Cw = concobs[k], time = tfin[j],
                                ks = ks[i], ke = ke[i],
                                NEC = nec[i],
                                m0 = m0[i])
      }
    }
  }
  
  dtheof <- do.call("rbind", dtheo)
  dtheof <- as.data.frame(cbind(rep(concobs, rep(npoints, length(concobs))),
                                rep(tfin, length(concobs)),
                                dtheof))
  names(dtheof) <- c("conc", "time", paste0("X", 1:length(sel)))
  
  # quantile
  qinf95 = NULL
  qsup95 = NULL
  
  for (i in 1:dim(dtheof)[1]) {
    qinf95[i] <- quantile(dtheof[i, 3:length(dtheof)],
                          probs = 0.025, na.rm = TRUE)
    qsup95[i] <- quantile(dtheof[i, 3:length(dtheof)],
                          probs = 0.975, na.rm = TRUE)
  }
  
  dtheof <- cbind(qinf95, qsup95, dtheof)
  
  return(dtheof)
}

survFitPlotTKTDGeneric <- function(data, xlab, ylab, main, one.plot, ci, dataCIm) {
  # vector color
  data[["dobs"]]$color <- as.numeric(as.factor(data[["dobs"]][["conc"]]))
  data[["dtheo"]]$color <- as.numeric(as.factor(data[["dtheo"]][["conc"]]))
  
  if (one.plot) {
    survFitPlotTKTDGenericOnePlot(data, xlab, ylab, main)
  } else {
    par(mfrow = plotMatrixGeometry(length(unique(data[["dobs"]][["conc"]]))))
    
    survFitPlotTKTDGenericNoOnePlot(data, xlab, ylab, main, ci, dataCIm)
    
    par(mfrow = c(1, 1))
  }
}

survFitPlotTKTDGenericOnePlot <- function(data, xlab, ylab, main) {
  plot(data[["dobs"]][["t"]],
       data[["dobs"]][["psurv"]],
       xlab = xlab,
       ylab = ylab,
       pch = 16,
       col = data[["dobs"]]$color,
       main = main)
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       lines(x$t, x$psurv, # lines
             col = x$color)
     })
}

survFitPlotTKTDGenericNoOnePlot <- function(data, xlab, ylab, main, ci, dataCIm) {
  if (ci) {
    survFitPlotTKTDGenericNoOnePlotCi(data, xlab, ylab, main, dataCIm)
  } else {
    survFitPlotTKTDGenericNoOnePlotNoCi(data, xlab, ylab, main)
  }
}

survFitPlotTKTDGenericNoOnePlotCi <- function(data, xlab, ylab, main, dataCIm) {
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       plot(x[, "t"],
            x[, "psurv"],
            xlab = xlab,
            ylab = ylab,
            type = "n",
            ylim = c(0, 1),
            col = x[, "color"],
            main = main)
       lines(x[, "t"], x[, "psurv"], # lines
             col = x[, "color"])
     })
}

survFitPlotTKTDGenericNoOnePlotNoCi <- function(data, xlab, ylab, main) {
  # one line by replicate
  by(data[["dtheo"]], list(data[["dtheo"]]$conc),
     function(x) {
       plot(x[, "t"],
            x[, "psurv"],
            xlab = xlab,
            ylab = ylab,
            type = "n",
            ylim = c(0, 1),
            col = x[, "color"],
            main = main)
       lines(x[, "t"], x[, "psurv"], # lines
             col = x[, "color"])
     })
}

survFitPlotTKTDGG <- function(data, xlab, ylab, main, one.plot, ci, dataCI,
                              dataCIm) {

  if (one.plot) {
    if (ci) warning("Credible intervals are only evalables in grid plot !")
    survFitPlotTKTDGGOnePlot(data, xlab, ylab, main)
  } else {
    survFitPlotTKTDGGNoOnePlot(data, xlab, ylab, main, ci, dataCI, dataCIm)
  }
}

survFitPlotTKTDGGOnePlot <- function(data, xlab, ylab, main) {
  ggplot(data$dobs, aes(x = t, y = psurv, colour = factor(conc))) +
    geom_point() + geom_line(data = data$dtheo) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
}

survFitPlotTKTDGGNoOnePlot <- function(data, xlab, ylab, main, ci, dataCI,
                                       dataCIm) {
  if (ci) {
    survFitPlotTKTDGGNoOnePlotCi(data, xlab, ylab, main, dataCI, dataCIm)
  } else {
    survFitPlotTKTDGGNoOnePlotNoCi(data, xlab, ylab, main)
  }
}

survFitPlotTKTDGGNoOnePlotCi <- function(data, xlab, ylab, main, dataCI,
                                         dataCIm) {
  ggplot(data$dobs,
         aes(x = t, y = psurv, colour = factor(conc))) +
    geom_line(data = dataCIm, aes(x = time, y = value, group = variable),
              alpha = 0.05) +
    geom_line(data = data$dtheo, color = "red") +
    geom_line(data = dataCI, aes(x = time, y = qinf95), linetype = 'dashed', color = "black") +
    geom_line(data = dataCI, aes(x = time, y = qsup95), linetype = 'dashed', color = "black") +
    facet_wrap(~conc) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal() +
    scale_color_discrete(guide = "none")
}

survFitPlotTKTDGGNoOnePlotNoCi <- function(data, xlab, ylab, main) {
  ggplot(data$dobs,
         aes(x = t, y = psurv, colour = factor(conc))) +
    geom_point() +
    geom_line(data = data$dtheo, colour = "red") +
    facet_wrap(~conc) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal() +
    scale_color_discrete(guide = "none")
}

