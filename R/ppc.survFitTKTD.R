#' Posterior predictive check plot for \code{survFitTKTD} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitTKTD} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{survFitTKTD} objects.
#' 
#' The black points show the observed number of survivors (pooled
#' replicates, on \eqn{X}-axis) against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise. Samples with equal observed
#' value are shifted on the \eqn{X}-axis. For that reason, the
#' bisecting line (y = x), is represented by steps when observed
#' values are low. That way we ensure green intervals do intersect the
#' bisecting line.
#'
#' @param x An object of class \code{survFitTKTD}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param main main title for the plot
#' @param \dots Further arguments to be passed to generic methods
#'
#' @examples
#'
#' # (1) Load the data
#' data(propiconazole)
#'
#' # (2) Create an object of class "survData"
#' dat <- survData(propiconazole)
#'
#' \dontrun{
#' # (3) Run the survFitTKTD function with the TKTD model ('SD' only)
#' out <- survFitTKTD(dat)
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
ppc.survFitTKTD <- function(x, style = "ggplot", main = NULL,...) {
  if (!is(x, "survFitTKTD"))
    stop("x is not of class 'survFitTKTD'!")
  
  xlab <- "Observed nb of survivors"
  ylab <- "Predicted nb of survivors"
  
  ppc_gen(EvalsurvTKTDPpc(x), style, xlab, ylab, main)
}

#' @importFrom stats rbinom quantile
EvalsurvTKTDPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)

  kd <- 10^(tot.mcmc[, "log10kd"])
  ks <- 10^(tot.mcmc[, "log10ks"])
  nec <- 10^(tot.mcmc[, "log10NEC"])
  m0 <- 10^(tot.mcmc[, "log10m0"])
  
  niter <- nrow(tot.mcmc)
  n <- x$jags.data$ndat
  xconc <- x$jags.data$x
  t <- x$jags.data$t
  tprec <- x$jags.data$tprec
  NsurvObs <- x$jags.data$y
  Nprec <- x$jags.data$Nprec
  bigtime <- x$jags.data$bigtime
  NsurvPred <- matrix(NA, nrow = niter, ncol = n)
  psurv = NULL
  for (i in 1:n) {
    for (j in 1:length(kd)) {
      xcor <- ifelse(xconc[i] > 0, xconc[i], 10)
      R <- ifelse(xconc[i] > nec[j], nec[j]/xcor, 0.1)
      tNEC <- ifelse(xconc[i] > nec[j], -1 / kd[j] * log(1 - R), bigtime)
      tref <- max(tprec[i], tNEC)
      psurv[j] <- exp(-m0[j] * (t[i] - tprec[i]) +
                        if (t[i] > tNEC) {
                          -ks[j] * ((xconc[i] - nec[j]) * (t[i] - tref) +
                                      xconc[i]/kd[j] * (exp(-kd[j] * t[i]) - exp(-kd[j] * tref)))
                        } else {
                          0
                        })
    }
    NsurvPred[, i] <- rbinom(niter, Nprec[i], psurv)
  }

  QNsurvPred <- t(apply(NsurvPred, 2, quantile,
                        probs = c(2.5, 50, 97.5) / 100, na.rm = TRUE))
  tab <- data.frame(QNsurvPred,
                    Nprec, NsurvObs,
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs |
                                   QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nprec", "Obs", "col")
  
  return(tab)
}

