#' @importFrom stats rbinom quantile
EvalsurvPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)
  
  if (x$det.part == "loglogisticbinom_3") {
    d <- sample(tot.mcmc[, "d"], 5000)
  }
  
  b <- 10^sample(tot.mcmc[, "log10b"], 5000)
  e <- 10^sample(tot.mcmc[, "log10e"], 5000)
  
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

#' @importFrom stats rgamma rpois quantile
EvalreproPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)
  
  if (x$model.label == "GP") {
    omega <- 10^sample(tot.mcmc[, "log10omega"], 5000)
  }
  b <- 10^sample(tot.mcmc[, "log10b"], 5000)
  d <- sample(tot.mcmc[, "d"], 5000)
  e <- 10^sample(tot.mcmc[, "log10e"], 5000)
  
  n <- x$jags.data$n
  xconc <- x$jags.data$xconc
  Nindtime <- x$jags.data$Nindtime
  NcumulObs <- x$jags.data$Ncumul
  NcumulPred <- matrix(NA, nrow = 5000, ncol = n)
  
  if (x$model.label == "GP") {
    for (i in 1:n) {
      popmean <- d / (1 + (xconc[i]/e)^b)
      indmean <- rgamma(n = 5000, shape = popmean / omega, rate = 1 / omega)
      NcumulPred[, i] <- rpois(5000, indmean * Nindtime[i])
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

#' @importFrom graphics abline arrows
PpcGeneric <- function(tab, xlab, ylab) {
  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
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
}

#' @import ggplot2
#' @importFrom  grid arrow unit
PpcGG <- function(tab, xlab, ylab) {
  
  ggplot(tab, aes(x = Obs, y = P50)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1) +
    geom_segment(aes(x = Obs, xend = Obs,
                     y = P2.5, yend = P97.5),
                 arrow = arrow(length = unit(0.25, "cm"), angle = 90,
                               ends = "both"),
                 tab, color = tab$col) +
    xlim(0, max(tab[, c("P97.5", "Obs")]) + 1) +
    ylim(0, max(tab[, c("P97.5", "Obs")]) + 1) +
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
    tab[, "Obs"] <- jitter(tab[, "Obs"])
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
    tab[, "Obs"] <- jitter(tab[, "Obs"])
    xlab <- "Observed Nbr. of cumulated offspring"
    ylab <- "Predicted Nbr. of cumulated offspring"
    
    if (style == "generic")
      PpcGeneric(tab, xlab, ylab)
    else if (style == "ggplot")
      PpcGG(tab, xlab, ylab)
    else stop("Unknown style")
  }
}
