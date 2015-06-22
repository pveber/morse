survPpcGeneric <- function(data, CI) {
  plot(data$obs, data$obs,
       type = "n",
       bty = "n",
       xaxt = "n",
       yaxt = "n",
       xlim = c(0, max(pretty(c(0, max(CI["qsup95",]))))),
       ylim = c(0, max(pretty(c(0, max(CI["qsup95",]))))),
       ylab = "Predicted survival rate",
       xlab = "Observed survival rate")
  abline(a = 0, b = 1)
  
  # axis
  axis(side = 2, at = pretty(c(0, max(CI["qsup95",]))))
  axis(side = 1, at = pretty(c(0, max(CI["qsup95",]))))
  
  points(data$obs, data$pred, pch = 16)
  for(i in 1:length(data$conc)) {
    segments(data$obs[i], CI["qinf95",][i], data$obs[i], CI["qsup95",][i],
             col = if (CI["qinf95",][i] > data$obs[i] | CI["qsup95",][i] < data$obs[i]) { "red" } else { "green" })
  }
}

#' @import ggplot2
survPpcGG <- function(data, CI) {
  data <- data.frame(data, qinf95 = CI["qinf95",],
                     qsup95 = CI["qsup95",],
                     color = ifelse(CI["qinf95",] > data$obs | CI["qsup95",] < data$obs,
                                    "red", "green"))

  ggplot(data, aes(x = obs, y = pred)) +
    geom_point() +
    geom_segment(aes(x = obs, xend = obs, y = qinf95, yend = qsup95),
                     data, color = data$color) +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Observed survival rate") +
    ylab("Predicted survival rate") +
    xlim(c(0, max(pretty(c(0, max(CI["qsup95",])))))) +
    ylim(c(0, max(pretty(c(0, max(CI["qsup95",])))))) +
    theme_minimal()
}

survPpc <- function(x, data, style) {
  data$obs <- data$Nsurv / data$Ninit
  # data points are systematically pooled, since our model does not
  # take individual variation into account
  data <- aggregate(obs ~ conc, data, mean)
  data$pred <- survEvalFit(x, data$conc)
  CI <- survLlbinomCI(x, log.scale = FALSE)
  
  if (style == "generic") {
    survPpcGeneric(data, CI)
  }
  else if (style == "ggplot") {
    survPpcGG(data, CI)
  }
  else stop("Unknown style")
}

reproPpcGeneric <- function(data, CI) {
  plot(data$obs, data$obs,
       type = "n",
       bty = "n",
       xaxt = "n",
       yaxt = "n",
       xlim = c(0, max(pretty(c(0, max(CI$qsup95))))),
       ylim = c(0, max(pretty(c(0, max(CI$qsup95))))),
       ylab = "Predicted reproduction rate",
       xlab = "Observed reproduction rate")
  abline(a = 0, b = 1)
  
  # axis
  axis(side = 2, at = pretty(c(0, max(CI$qsup95))))
  axis(side = 1, at = pretty(c(0, max(CI$qsup95))))
  
  points(data$obs, data$pred, pch = 16)
  for (i in 1:length(data$conc)) {
    segments(data$obs[i], CI$qinf95[i], data$obs[i], CI$qsup95[i],
             col = if (CI$qinf95[i] > data$obs[i] | CI$qsup95[i] < data$obs[i]) { "red" } else { "green" })
  }
}

#' @import ggplot2
reproPpcGG <- function(data, CI) {
  data <- data.frame(data, qinf95 = CI$qinf95,
                     qsup95 = CI$qsup95,
                     color = ifelse(CI$qinf95 > data$obs | CI$qsup95 < data$obs, "red", "green"))
  ggplot(data, aes(x = obs, y = pred)) +
    geom_point() +
    geom_segment(aes(x = obs, xend = obs, y = qinf95, yend = qsup95),
                 data, color = data$color) +
    geom_abline(intercept = 0, slope = 1) +
    xlab("Observed reproduction rate") +
    ylab("Predicted reproduction rate") +
    xlim(c(0, max(pretty(c(0, max(CI$qsup95)))))) +
    ylim(c(0, max(pretty(c(0, max(CI$qsup95)))))) +
    theme_minimal()
}

reproPpc <- function(x, data, style) {
  data$obs <- data$Nreprocumul / data$Nindtime
  data$pred <- reproEvalFit(x, data$conc)
  CI <- reproLlmCI(x, data$conc)
  
  if (style == "generic") {
    reproPpcGeneric(data, CI)
  }
  else if (style == "ggplot") {
    reproPpcGG(data, CI)
  }
  else stop("Unknown style")
}

#' Predicted vs obseved plot function
#' 
#' The \code{ppc} function plots the predicted versus observed response values.
#' 
#' @param x An object of class \code{survFitTT} or \code{reproFitTT}
#' @param style generic or ggplot
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
#' out <- reproFitTT(dat2, stoc.part = "gammapoisson", 
#' ecx = c(5, 10, 15, 20, 30, 50, 80), quiet = TRUE)
#' 
#' # (4) Plot the observed vs predicted values with generic style
#' ppc(out)
#' 
#' # (4) Plot the observed vs predicted values with ggplot2 style
#' ppc(out, style = "ggplot")
#' }
#' 
#' @export
#' @import ggplot2
#' 
ppc <- function(x, style = "generic") {
  # test class object
  if (! is(x, "reproFitTT") && ! is(x, "survFitTT"))
    stop("The object passed in argument [out] is not of class 'reproFitTT' or 'survFitTT' !\n")
  
  dataTT <- x$dataTT
  
  if (is(x, "reproFitTT")) {
    reproPpc(x, dataTT, style)
  }
  else if (is(x, "survFitTT")) {
    survPpc(x, dataTT, style)
  }
}