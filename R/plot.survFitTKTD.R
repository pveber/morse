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
  k <- nrow(mctot)
  ks <- 10^mctot[, "lks"]
  ke <- 10^mctot[, "lke"]
  m0 <- 10^mctot[, "lm0"]
  nec <- 10^mctot[, "lNEC"]
  
  dtheo <- array(data = NA, dim = c(npoints, length(nec), length(concobs)))
  for (k in 1:length(concobs)) {
    for (i in 1:length(nec)) {
      for (j in 1:npoints) {
        dtheo[j, i, k] <- Surv(Cw = concobs[k], time = tfin[j],
                                ks = ks[i], ke = ke[i],
                                NEC = nec[i],
                                m0 = m0[i])
#         
#         
#         dtheo[[k]] <- rbind(dtheo[[k]], data.frame(conc = concobs[i],
#                                                      t = tfin[j],
#                                                      psurv = psurv,
#                                                      q = k))
      }
    }
  }
  
  dtheof <- rbind(dtheo[,, 1], dtheo[,, 2], dtheo[,, 3],
                  dtheo[,, 4], dtheo[,, 5])
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     t = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init)
  
  return(list(dtheo = dtheof,
              dobs = dobs))
}

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
  
  if (CI) dataCI <- survFitPlotCITKTD(x)
  
  
  if (style == "generic") {
    # vector color
    data[["dobs"]]$color <- as.numeric(as.factor(data[["dobs"]][["conc"]]))
    data[["dtheo"]]$color <- as.numeric(as.factor(data[["dtheo"]][["conc"]]))
    
    if (one.plot) {
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
    } else {
      par(mfrow = plotMatrixGeometry(length(unique(data[["dobs"]][["conc"]]))))
      
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
      par(mfrow = c(1, 1))
    }
  }
  
  if (style == "ggplot") {
    if (one.plot) {
      plt1 <- ggplot(data$dobs, aes(x = t, y = psurv, colour = factor(conc))) +
        labs(x = xlab, y = ylab) + ggtitle(main) +
        ylim(c(0, 1)) +
        geom_point() + geom_line(data = data$dtheo) + theme_minimal()
    } else {
      plt1 <- ggplot(data$dobs, aes(x = t, y = psurv, colour = factor(conc))) +
        facet_wrap(~conc) +
        labs(x = xlab, y = ylab) + ggtitle(main) +
        ylim(c(0, 1)) +
        
        geom_point() + geom_line(data = data$dtheo) + theme_minimal()
    }
    
    plt1
  }
}
