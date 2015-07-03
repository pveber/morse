EvalsurvPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)
  
  if (x$det.part == "loglogisticbinom_3") {
    d <- tot.mcmc[1:5000, "d"]
  }
  b <- 10^tot.mcmc[1:5000, "log10b"]
  e <- 10^tot.mcmc[1:5000, "log10e"]
  
  n <- x$jags.data$n
  xconc <- x$jags.data$xconc
  Ninit <- x$jags.data$Ninit
  NsurvObs <- x$jags.data$Nsurv
  NsurvPred <- matrix(NA, nrow = 5000, ncol = n)
  
  if (x$det.part == "loglogisticbinom_2") {
    for (i in 1:n) {
      p <- 1 / (1 + (xconc[i]/e)^b)
      NsurvPred[, i] <- rbinom(5000, Ninit[i], p)
    }
  }
  if (x$det.part == "loglogisticbinom_3") {
    for (i in 1:n) {
      p <- d / (1 + (xconc[i]/e)^b)
      NsurvPred[, i] <- rbinom(5000, Ninit[i], p)
    }
  }
  QNsurvPred <- t(apply(NsurvPred, 2, quantile,
                        probs = c(2.5, 50, 97.5) / 100))
  tab <- data.frame(QNsurvPred,
                    Ninit, NsurvObs,
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs | QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Ninit", "Obs", "col")
  
  return(tab)
}

EvalreproPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)
  
  if (x$model.label == "GP") {
    omega <- 10^tot.mcmc[1:5000, "log10omega"]
  }
  b <- 10^tot.mcmc[1:5000, "log10b"]
  d <- tot.mcmc[1:5000, "d"]
  e <- 10^tot.mcmc[1:5000, "log10e"]
  
  n <- x$jags.data$n
  xconc <- x$jags.data$xconc
  Nindtime <- x$jags.data$Nindtime
  NcumulObs <- x$jags.data$Ncumul
  NcumulPred <- matrix(NA, nrow = 5000, ncol = n)
  
  if (x$model.label == "GP") {
    for (i in 1:n) {
      theomean <- d / (1 + (xconc[i]/e)^b)
      nbtheo <- theomean * Nindtime[i]
      NcumulPred[, i] <- rgamma(n = 5000, shape = nbtheo / omega, rate = 1 / omega)
    }
    
  }
  if (x$model.label == "P") {
    for (i in 1:n) {
      ytheo <- d / (1 + (xconc[i]/e)^b)
      nbtheo <- ytheo * Nindtime[i]
      NcumulPred[, i] <- rpois(5000, nbtheo)
    }
  }
  QNreproPred <- t(apply(NcumulPred, 2, quantile,
                         probs = c(2.5, 50, 97.5) / 100))
  tab <- data.frame(QNreproPred,
                    Nindtime, NcumulObs,
                    col = ifelse(QNreproPred[,"2.5%"] > NcumulObs | QNreproPred[,"97.5%"] < NcumulObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nindtime", "Obs", "col")
  
  return(tab)
}

PpcGeneric <- function(tab, xlab, ylab) {
  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
       xaxt = "n",
       yaxt = "n",
       xlab = xlab,
       ylab = ylab)
  
  points(tab[, "Obs"], tab[, "P50"],
         pch = 16)
  
  abline(0, 1)
  
  arrows(tab[, "Obs"], tab[, "P50"],
         tab[, "Obs"], tab[, "P2.5"],
         angle = 90, length = 0.05,
         col = as.character(tab[, "col"]))
  arrows(tab[, "Obs"], tab[, "P50"],
         tab[, "Obs"], tab[, "P97.5"],
         angle = 90, length = 0.05,
         col = as.character(tab[, "col"]))
  
  # axis
  axis(side = 1, at = c(0, unique(tab[, "Obs"])))
  axis(side = 2, at = c(0, unique(tab[, "P50"])))
}

#' @import ggplot2
#' @import grid
PpcGG <- function(tab, xlab, ylab) {
  
  ggplot(tab, aes(x = Obs, y = P50)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    xlim(c(0, max(tab[, "P97.5"]))) +
    ylim(c(0, max(tab[, "P97.5"]))) +
    scale_x_continuous(breaks = c(0, unique(tab[, "Obs"]))) +
    scale_y_continuous(breaks = c(0, unique(tab[, "P50"]))) +
    geom_segment(aes(x = Obs, xend = Obs,
                     y = P2.5, yend = P97.5),
                 arrow = arrow(length = unit(0.25, "cm"), angle = 90,
                               ends = "both"),
                 tab, color = tab$col) +
    labs(x = xlab, y = ylab) +
    theme_minimal()
}

#' Posterior predictive check plot
#' 
#' The \code{ppc} functions plot the observed versus predicted values for the
#' \code{survFitTT} and \code{reporFitTT} objects.
#' 
#' @param x An object of class \code{reproFitTT} or \code{survFitTT}
#' @param style Graphical package method: \code{generic} or \code{ggplot}.
#' 
#' @examples
#' 
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create an object of class "reproData"
#' dat <- reproData(cadmium1)
#' 
#' \dontrun{
#' # (3) Run the reproFitTT function with the log-logistic gamma-poisson model
#' out <- reproFitTT(dat, stoc.part = "gammapoisson", 
#' ecx = c(5, 10, 15, 20, 30, 50, 80), quiet = TRUE)
#' 
#' # (4) Plot observed versus predicted values
#' ppc(out)
#' }
#' 
#' @export
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#' 
ppc <- function(x, style = "generic") {
  
  if (!is(x, "survFitTT") & !is(x, "reproFitTT"))
    stop("x is not of class 'survFitTT' or 'reproFitTT' !")
  
  if (is(x, "survFitTT")) {
    tab <- EvalsurvPpc(x)
    xlab <- "Observed Nbr. of survivor"
    ylab <- "Predicted Nbr. of survivor"
    
    if (style == "generic")
      PpcGeneric(tab, xlab, ylab)
    else if (style == "ggplot")
      PpcGG(tab, xlab, ylab)
    else stop("Unknown style")
  }
  
  else if (is(x, "reproFitTT")) {
    tab <- EvalreproPpc(x)
    xlab <- "Observed Nbr. of cumulated offspring"
    ylab <- "Predicted Nbr. of cumulated offspring"
    
    if (style == "generic")
      PpcGeneric(tab, xlab, ylab)
    else if (style == "ggplot")
      PpcGG(tab, xlab, ylab)
    else stop("Unknown style")
  }
}
