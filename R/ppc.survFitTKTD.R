#' Posterior predictive check plot for \code{survFitTKTD} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitTKTD} class. It
#' plots the predicted values with 95 \% credible intervals
#' versus the observed values for \code{survFitTKTD} objects.
#' 
#' The coordinates of black points are the observed values of the number of survivor
#' (poolled replicates) for a given concentration (x-scale) and the corresponding 
#' predicted values (y-scale). 95 \% prediction intervals are added to each predicted
#' value, colored in green if this interval contains the observed value and in red
#' in the other case.
#' As replicates are shifted on the x-axis, the bisecting line (y = x), is
#' represented by steps, and is added to the plot in order to see if each
#' prediction interval contains each observed value. 
#'
#' @param x An object of class \code{survFitTKTD}
#' \code{'ggplot'} style to avoid the label overlap
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
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
#' # (3) Run the survFitTKTD function with the TKTD model
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
ppc.survFitTKTD <- function(x, meth = "for", style = "generic", ...) {
  if (!is(x, "survFitTKTD"))
    stop("x is not of class 'survFitTKTD'!")
  
  xlab <- "Observed Nbr. of survivor"
  ylab <- "Predicted Nbr. of survivor"
  
  ppc_gen(EvalsurvTKTDPpc(x, meth), style, xlab, ylab)
}

#' @importFrom stats rbinom quantile
EvalsurvTKTDPpc <- function(x, meth) {
  tot.mcmc <- do.call("rbind", x$mcmc)

  kd <- 10^sample(tot.mcmc[, "log10kd"], 5000)
  ks <- 10^sample(tot.mcmc[, "log10ks"], 5000)
  nec <- 10^sample(tot.mcmc[, "log10NEC"], 5000)
  m0 <- 10^sample(tot.mcmc[, "log10m0"], 5000)
  
  n <- x$jags.data$ndat
  xconc <- x$jags.data$x
  t <- x$jags.data$t
  tprec <- x$jags.data$tprec
  Nprec <- x$jags.data$Nprec
  NsurvObs <- x$jags.data$y
  Nprec <- x$jags.data$Nprec
  bigtime <- x$jags.data$bigtime
  NsurvPred <- matrix(NA, nrow = 5000, ncol = n)
  
  if (meth == "for") {
    for (i in 1:n) {
      for (j in 1:length(kd)) {
        xcor <- ifelse(xconc[i] > 0, xconc[i], 10)
        R <- ifelse(xconc[i] > nec[j], nec[j]/xcor, 0.1)
        tNEC <- ifelse(xconc[i] > nec[j], -1 / kd[j] * log(1 - R), bigtime)
        tref <- max(tprec[i], tNEC)
        psurv <- exp(-m0 * (t[i] - tprec[i]) +
                       if (t[i] > tNEC) {
                         -ks * ((xconc[i] - nec[j]) * (t[i] - tref) +
                                  xconc[i]/kd[j] * (exp(-kd[j] * t[i]) - exp(-kd[j] * tref)))
                       } else {
                         0
                       })
      }
      NsurvPred[, i] <- rbinom(5000, Nprec[i], psurv)
    }
  } else {
    i <- 1:n
    j <- 1:length(kd)
    NsurvPred[, 1] <- sapply(j, function(x) {
        xcor <- ifelse(xconc[1] > 0, xconc[1], 10)
        R <- ifelse(xconc[1] > nec[x], nec[x]/xcor, 0.1)
        tNEC <- ifelse(xconc[1] > nec[x], -1 / kd[x] * log(1 - R), bigtime)
        tref <- max(tprec[1], tNEC)
        psurv <- exp(-m0 * (t[1] - tprec[1]) +
                       if (t[1] > tNEC) {
                         -ks * ((xconc[1] - nec[x]) * (t[1] - tref) +
                                  xconc[1]/kd[x] * (exp(-kd[x] * t[1]) - exp(-kd[x] * tref)))
                       } else {
                         0
                       })
        rbinom(5000, Nprec[1], psurv)
    }, y = i, x = j)
    return(NsurvPred)
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

