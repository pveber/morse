#' Posterior predictive check plot for survFitTT objects
#'
#' The \code{ppc} functions plot the observed versus predicted values for the
#' \code{survFitTT} and \code{reporFitTT} objects.
#'
#' @param x An object of class \code{reproFitTT} or \code{survFitTT}
#' @param style Graphical package method: \code{generic} or \code{ggplot}.
#' @param \dots Further arguments to be passed to generic methods.
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
ppc.survFitTT <- function(x, style = "generic", ...) {
  if (!is(x, "survFitTT"))
    stop("x is not of class 'survFitTT'!")

  ppc_gen(EvalsurvPpc(x), style)
}

ppc_gen <- function(tab, style) {
  xlab <- "Observed Nbr. of survivor"
  ylab <- "Predicted Nbr. of survivor"

  if (style == "generic") PpcGeneric(tab, xlab, ylab)
  else if (style == "ggplot") PpcGG(tab, xlab, ylab)
  else stop("Unknown style")
}

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
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs |
                                   QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Ninit", "Obs", "col")

  return(tab)
}

stepCalc <- function(obs_val) {
  # calculation of steps coordinate
  sObs <- sort(c(0, obs_val))
  stepX <- c(0, sapply(2:length(sObs), function(i) {
    sObs[i-1] + (sObs[i] - sObs[i-1]) / 2}), max(sObs))
  return(list(sObs = sObs, stepX = stepX))
}

jitterObsGenerator <- function(stepX, tab, obs_val) {
  # uniform jittering of observed values
  spaceX <- min(sapply(2:length(stepX),
                       function(i) { stepX[i] - stepX[i-1] }))/2
  lengthX <- table(tab[, "Obs"])
  return(unlist(mapply(function(x, y) {
    seq(x - spaceX, x + spaceX, length.out = y)
  }, x = sort(obs_val), y = lengthX)))
}

#' @importFrom graphics abline segments
PpcGeneric <- function(tab, xlab, ylab) {
  obs_val <- unique(tab[, "Obs"])
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  jittered_obs <- jitterObsGenerator(stepX, tab, obs_val)

  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
       xlab = xlab,
       ylab = ylab)
  
  sapply(2:length(stepX), function(i) {
    segments(stepX[i-1], sObs[i-1], stepX[i], sObs[i-1])
    segments(stepX[i], sObs[i-1], stepX[i], sObs[i])
  })
  
  tab0 <- tab[order(tab$Obs),]
  delta <- 0.01 * (max(obs_val) - min(obs_val))
  segments(jittered_obs, tab0[, "P2.5"],
           jittered_obs, tab0[, "P97.5"],
           col = as.character(tab0[, "col"]))
  segments(jittered_obs - delta, tab0[, "P2.5"],
           jittered_obs + delta, tab0[, "P2.5"],
           col = as.character(tab0[, "col"]))
  segments(jittered_obs - delta, tab0[, "P97.5"],
           jittered_obs + delta, tab0[, "P97.5"],
           col = as.character(tab0[, "col"]))

  points(jittered_obs, tab0[, "P50"],
         pch = 16)
}

#' @import ggplot2
#' @importFrom  grid arrow unit
PpcGG <- function(tab, xlab, ylab) {
  obs_val <- unique(tab[, "Obs"])
  tab$jittered_obs <- jitter(tab[, "Obs"])
  
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  
  df <- data.frame(sObs = c(sObs, max(sObs)),
                   stepX)
  
  ggplot(df, aes(x = stepX, y = sObs)) +
    geom_step() +
    geom_segment(aes(x = jittered_obs, xend = jittered_obs,
                     y = P2.5, yend = P97.5), tab,
                 arrow = arrow(length = unit(0.1, "cm"), angle = 90,
                               ends = "both"),
                 color = tab$col) +
    geom_point(aes(x = jittered_obs, y = P50), tab) +
    xlim(0, max(tab[, c("P97.5", "Obs", "jittered_obs")]) + 1) +
    ylim(0, max(tab[, c("P97.5", "Obs", "jittered_obs")]) + 1) +
     labs(x = xlab, y = ylab) +
     theme_minimal()
}

