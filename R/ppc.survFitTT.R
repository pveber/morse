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

  if (style == "generic") {
    PpcGeneric(tab, xlab, ylab)
  } else if (style == "ggplot") {
    PpcGG(tab, xlab, ylab)
  } else stop("Unknown style")
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

#' @importFrom graphics abline segments
PpcGeneric <- function(tab, xlab, ylab) {
  obs_val <- unique(tab[, "Obs"])
  jittered_obs <- jitter(tab[, "Obs"])

  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
       xlab = xlab,
       ylab = ylab)

  segments(obs_val - 0.5, obs_val,
           obs_val + 0.5, obs_val)

  delta <- 0.01 * (max(obs_val) - min(obs_val))
  segments(jittered_obs,tab[, "P2.5"],
           jittered_obs,tab[, "P97.5"],
           col = as.character(tab[, "col"]))
  segments(jittered_obs - delta, tab[, "P2.5"],
           jittered_obs + delta, tab[, "P2.5"],
           col = as.character(tab[, "col"]))
  segments(jittered_obs - delta, tab[, "P97.5"],
           jittered_obs + delta, tab[, "P97.5"],
           col = as.character(tab[, "col"]))

  points(jittered_obs, tab[, "P50"],
         pch = 16)
}

#' @import ggplot2
#' @importFrom  grid arrow unit
PpcGG <- function(tab, xlab, ylab) {
  tab$jittered_obs <- jitter(tab[, "Obs"])
  delta <- 0.01 * (max(tab[, "Obs"]) - min(tab[, "Obs"]))

  ggplot(tab) +
    geom_segment(aes(x = jittered_obs, xend = jittered_obs,
                     y = P2.5, yend = P97.5), tab,
                 arrow = arrow(length = unit(0.25, "cm"), angle = 90,
                               ends = "both"),
                 color = tab$col) +
    geom_point(aes(x = jittered_obs, y = P50), tab) +
    xlim(0, max(tab[, c("P97.5", "Obs", "jittered_obs")]) + 1) +
    ylim(0, max(tab[, c("P97.5", "Obs", "jittered_obs")]) + 1) +
     labs(x = xlab, y = ylab) +
     theme_minimal()
}

