survPpcGeneric <- function(tab, xlab, ylab) {
  plot(c(0, max(tab[, "Ninit"])),
       c(0, max(tab[, "Ninit"])),
       type = "n",
       xaxt = "n",
       yaxt = "n",
       xlab = xlab,
       ylab = ylab)
  
  points(tab[, "NsurvObs"], tab[, "P50"],
         pch = 16)
  
  abline(0, 1)
  
  for (i in 1:length(tab[, "NsurvObs"])) {
    arrows(tab[i, "NsurvObs"], tab[i, "P50"],
           tab[i, "NsurvObs"], tab[i, "P2.5"],
           angle = 90, length = 0.1,
           col = if(tab[i, "P2.5"] > tab[i, "NsurvObs"] | tab[i, "P97.5"] < tab[i, "NsurvObs"]) { "red" } else {"green"})
    arrows(NsurvObs[i], tab[i, "P50"],
           NsurvObs[i], tab[i, "P97.5"],
           angle = 90, length = 0.1,
           col = if(tab[i, "P2.5"] > tab[i, "NsurvObs"] | tab[i, "P97.5"] < tab[i, "NsurvObs"]) { "red" } else {"green"})
  }
  
  # axis
  axis(side = 1, at = pretty(c(0, max(Ninit))))
  axis(side = 2, at = pretty(c(0, max(Ninit))))
}

#' @import ggplot2
#' @import grid
survPpcGG <- function(tab, xlab, ylab) {

  ggplot(tab, aes(x = NsurvObs, y = P50)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    xlim(c(0, max(tab$Ninit))) +
    ylim(c(0, max(tab$Ninit))) +
    geom_segment(aes(x = NsurvObs, xend = NsurvObs,
                     y = P50, yend = P2.5),
                 arrow = arrow(length = unit(0.5, "cm"), angle = 90),
                 tab, color = tab$col) +
    geom_segment(aes(x = NsurvObs, xend = NsurvObs,
                     y = P50, yend = P97.5),
                 arrow = arrow(length = unit(0.5, "cm"), angle = 90),
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
  if (is(x, "survFitTT")) {
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
    colnames(tab) <- c("P2.5", "P50", "P97.5", "Ninit", "NsurvObs", "col")
    
    
    xlab <- "Observed Nbr. of survivor"
    ylab <- "Predicted Nbr. of survivor"
    
    if (style == "generic") {
      survPpcGeneric(tab, xlab, ylab)
    }
    else if (style == "ggplot") {
      survPpcGG(tab, xlab, ylab)
    }
    else stop("Unknown style")
  }
  else if (is(x, "reproFitTT")) {
    tot.mcmc <- do.call("rbind", x$mcmc)
    
    if (x$dmodel.label == "GP") {
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
        rate <- d / (1 + (xconc[i]/e)^b) / omega
        p <- 1 / (Nindtime[i] * omega + 1)
        NcumulPred[, i] ~ rnbinom(5000, Ninit[i], p, rate)
      }
    }
    
  }
  else stop("x is not of class 'survFitTT' or 'reproFitTT' !")
}
