#' Posterior predictive check plot for \code{survFitTT} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitTT} class. It
#' plots the predicted values with 95 \% credible intervals versus the observed
#' values for \code{survFitTT} objects.
#' 
#' The coordinates of black points are the observed values of the number of survivors
#' (pooled replicates) for a given concentration (\eqn{X}-axis) and the corresponding 
#' predicted values (\eqn{Y}-axis). 95\% prediction intervals are added to each predicted
#' value, colored in green if this interval contains the observed value and in red
#' otherwise.
#' The bisecting line (y = x) is added to the plot in order to see if each
#' prediction interval contains each observed value. As replicates are shifted
#' on the x-axis, this line is represented by steps.
#'
#' @rdname PPC
#'
#' @param x An object of class \code{survFitTT}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param main main title for the plot
#' @param \dots Further arguments to be passed to generic methods
#'
#' @examples
#'
#' # (1) Load the data
#' data(cadmium1)
#'
#' # (2) Create an object of class "survData"
#' dat <- survData(cadmium1)
#'
#' \dontrun{
#' # (3) Run the survFitTT function with the log-logistic binomial model
#' out <- survFitTT(dat, lcx = c(5, 10, 15, 20, 30, 50, 80),
#' quiet = TRUE)
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
ppc.survFitTT <- function(x, style = "ggplot", main = NULL, ...) {
  if (!is(x, "survFitTT"))
    stop("x is not of class 'survFitTT'!")
  
  xlab <- "Observed nb of survivors"
  ylab <- "Predicted nb of survivors"
  
  ppc_gen(EvalsurvPpc(x), style, xlab, ylab, main)
}

ppc_gen <- function(tab, style, xlab, ylab, main) {
  
  if (style == "generic") PpcGeneric(tab, xlab, ylab, main)
  else if (style == "ggplot") PpcGG(tab, xlab, ylab, main)
  else stop("Unknown style")
}

#' @importFrom stats rbinom quantile
EvalsurvPpc <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)
  
  if (x$det.part == "loglogisticbinom_3") {
    d <- tot.mcmc[, "d"]
  }
  
  b <- 10^tot.mcmc[, "log10b"]
  e <- 10^tot.mcmc[, "log10e"]
  
  niter <- nrow(tot.mcmc)
  n <- x$jags.data$n
  xconc <- x$jags.data$xconc
  Ninit <- x$jags.data$Ninit
  NsurvObs <- x$jags.data$Nsurv
  NsurvPred <- matrix(NA, nrow = niter, ncol = n)
  
  if (x$det.part == "loglogisticbinom_2") {
    for (i in 1:n) {
      p <- 1 / (1 + (xconc[i]/e)^b)
      NsurvPred[, i] <- rbinom(niter, Ninit[i], p)
    }
  }
  if (x$det.part == "loglogisticbinom_3") {
    for (i in 1:n) {
      p <- d / (1 + (xconc[i]/e)^b)
      NsurvPred[, i] <- rbinom(niter, Ninit[i], p)
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
PpcGeneric <- function(tab, xlab, ylab, main) {
  obs_val <- unique(tab[, "Obs"])
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  jittered_obs <- jitterObsGenerator(stepX, tab, obs_val, ppc = TRUE)$jitterObs
  spaceX <- jitterObsGenerator(stepX, tab, obs_val, ppc = TRUE)$spaceX
  
  plot(c(0, max(tab[, "P97.5"])),
       c(0, max(tab[, "P97.5"])),
       type = "n",
       xlab = xlab,
       ylab = ylab,
       main = main,
       xaxt = "n",
       yaxt = "n")
  
  # axis
  axis(side = 2, at = if (max(tab[, "Obs"]) == 1) {
    c(0, 1)
  } else {
    pretty(c(0, max(tab[, "P97.5"])))
  })
  axis(side = 1, at = if (max(tab[, "Obs"]) == 1) {
    c(0, 1)
  } else {
    pretty(c(0, max(tab[, "P97.5"])))
  })
  
  if (max(sObs) < 20) {
    sapply(1:length(sObs), function(i) {
      segments(sObs[i] - (spaceX * 1.25), sObs[i],
               sObs[i] + (spaceX * 1.25), sObs[i])
    })
  } else {
    abline(0, 1)
  }
  
  tab0 <- tab[order(tab$Obs),]
  segments(jittered_obs, tab0[, "P2.5"],
           jittered_obs, tab0[, "P97.5"],
           col = as.character(tab0[, "col"]))
  
  points(jittered_obs, tab0[, "P50"],
         pch = 20)
}

#' @import ggplot2
#' @importFrom  grid arrow unit
PpcGG <- function(tab, xlab, ylab, main) {
  obs_val <- unique(tab[, "Obs"])
  sObs <- stepCalc(obs_val)$sObs
  stepX <- stepCalc(obs_val)$stepX
  jittered_obs <- jitterObsGenerator(stepX, tab, obs_val, ppc = TRUE)$jitterObs
  spaceX <- jitterObsGenerator(stepX, tab, obs_val, ppc = TRUE)$spaceX
  
  tab0 <- cbind(tab[order(tab$Obs),], jittered_obs)
  
  df <- data.frame(sObs, spaceX)

  if (max(sObs) < 20) {
    gf1 <- ggplot(df) +
      geom_segment(aes(x = sObs - (spaceX * 1.25),
                       xend = sObs + (spaceX * 1.25),
                       y = sObs, yend = sObs))
  } else {
    gf1 <- ggplot(tab0) +
      geom_abline(intercept = 0, slope = 1)
  }
  
  gf2 <- gf1 +
    geom_segment(aes(x = jittered_obs, xend = jittered_obs,
                     y = P2.5, yend = P97.5), data = tab0,
                 color = tab0$col) +
    geom_point(aes(x = jittered_obs, y = P50), tab0) +
    expand_limits(y = 0) +
    expand_limits(x = 0) +
    labs(x = xlab, y = ylab, title = main) +
    theme_minimal()
  
  return(gf2)
}
