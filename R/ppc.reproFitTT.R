#' Posterior predictive check plot for reproFitTT objects
#'
#' The \code{ppc} function plots the predicted values with 95 \% credible intervals
#' versus the observed values for \code{reproFitTT} objects.
#' 
#' The coordinates of black points are the observed values of the cumulated number
#' of reproduction outputs for a given concentration (x-scale) and the corresponding 
#' predicted values (y-scale). 95 \% prediction intervals are added to each predicted
#' value, colored in green if this interval contains the observed value and in red
#' in the other case.
#' As replicates are not pooled in this plot, overlapped points
#' are shifted on the x-axis to help the visualization of replicates.
#' The bisecting line (y = x), is represented by steps, and is added to the plot
#' in order to see if each prediction interval contains each observed value.
#'
#' @param x An object of class \code{reproFitTT}
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid the label overlap
#' @param style Graphical package method: \code{generic} or \code{ggplot}
#' @param \dots Further arguments to be passed to generic methods
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
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
#' @export
ppc.reproFitTT <- function(x, remove.someLabels = FALSE,
                           style = "generic", ...) {
  if (!is(x, "reproFitTT"))
    stop("x is not of class 'reproFitTT'!")
  
  xlab <- "Observed Cumul. Nbr. of offspring"
  ylab <- "Predicted Cumul. Nbr. of offspring"

  ppc_gen(EvalreproPpc(x), style, xlab, ylab, remove.someLabels)
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
