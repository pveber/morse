#' Plotting method for \code{reproData} objects
#'
#' This is the generic \code{plotDoseResponse} S3 method for the \code{reproData}
#' class. It plots the number of offspring per individual-days as a function of
#' concentration (for a given target time).
#' 
#' The function plots the observed values of the reproduction rate (number of
#' reproduction outputs per individual-day) for a given time as a function of
#' concentration. The 95 \% Poisson confidence interval is added to each reproduction
#' rate. It is calculated using function \code{\link[epitools]{pois.exact}}
#' from package \code{epitools}.
#' As replicates are not pooled in this plot, overlapped points are shifted on
#' the x-axis to help the visualization of replicates.
#'
#' @param x an object of class \code{reproData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a label for the \eqn{y}-axis
#' @param main main title for the plot
#' @param ylim Y-axis limits
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid the label overlap
#' @param axis if \code{TRUE} displays ticks and label axis
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param rect if \code{TRUE}, confidence interval on data are represented by 
#' rectangles, else by segments
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @note When \code{style = "ggplot"}, the function calls function
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' 
#' @seealso \code{\link[epitools]{pois.exact}}
#'
#' @keywords plot
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(zinc)
#' zinc <- reproData(zinc)
#'
#' # (2) Plot dose-response
#' plotDoseResponse(zinc)
#'
#' # (3) Plot dose-response with a ggplot style
#' plotDoseResponse(zinc, style = "ggplot")
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom graphics plot axis lines points
#' title
#' @importFrom methods is
#' @importFrom stats aggregate
#' @importFrom epitools pois.exact
#'
#' @export
plotDoseResponse.reproData <- function(x,
                                       xlab = "Concentration",
                                       ylab = "Nb of offspring / Nb individual-days",
                                       main = NULL,
                                       ylim = NULL,
                                       target.time = NULL,
                                       style = "generic",
                                       log.scale = FALSE,
                                       remove.someLabels = FALSE,
                                       axis = TRUE,
                                       addlegend = TRUE,
                                       rect = FALSE,
                                       ...) {
  if (is.null(target.time)) target.time <- max(x$time)
  
  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")
  
  x$resp <- x$Nreprocumul / x$Nindtime
  
  # select the target.time
  xf <- filter(x, x$time == target.time)
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) xf$conc > 0 else TRUE
  x <- xf[sel, ]
  transf_data_conc <- optLogTransform(log.scale, x$conc)
  
  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, x$conc)
    if(log.scale) exp(x) else x
  })()
  
  ICpois <- pois.exact(x$Nreprocumul, x$Nindtime)
  x$reproRateInf <- ICpois$lower
  x$reproRateSup <- ICpois$upper
  conc_val <- unique(transf_data_conc)
  x$Obs <- x$conc
  stepX <- stepCalc(conc_val)$stepX
  jittered_conc <- jitterObsGenerator(stepX, x, conc_val, log.scale)$jitterObs
  
  if (style == "generic")
    reproDoseResponseCIGeneric(x, conc_val, jittered_conc, transf_data_conc,
                               display.conc, ylim, axis, main, addlegend, rect)
  else if (style == "ggplot")
    reproDoseResponseCIGG(x, conc_val, jittered_conc, transf_data_conc,
                          display.conc, main, addlegend, remove.someLabels,
                          rect)
  else stop("Unknown style")
}

reproDoseResponseCIGeneric <- function(x, conc_val, jittered_conc,
                                       transf_data_conc, display.conc, ylim,
                                       axis, main, addlegend, rect) {
  
  if (is.null(ylim)) ylim <- c(0, max(x$reproRateSup))
  plot(jittered_conc, x$resp,
       ylim = ylim,
       type = "n",
       xaxt = "n",
       yaxt = "n",
       main = main,
       xlab = if (axis) { "Concentration" } else "",
       ylab = if (axis) { "Reproduction rate"} else "")
  
  # axis
  if (axis) {
    axis(side = 2, at = pretty(c(0, max(x$resp))))
    axis(side = 1, at = transf_data_conc,
         labels = display.conc)
  }
  
  x0 <- x[order(x$conc),]
  delta <- 0.01 * (max(conc_val) - min(conc_val))
  if (rect) {
    color <- "gray"
    color_transparent <- adjustcolor(color, alpha.f = 0.3)
    
    points(sort(transf_data_conc), x0$resp, pch = 20)
    
    rect(sort(transf_data_conc) - delta, x0[, "reproRateInf"],
         sort(transf_data_conc) + delta, x0[, "reproRateSup"],
         col = color_transparent, border = NA)
  } else {
    segments(jittered_conc, x0[, "reproRateInf"],
             jittered_conc, x0[, "reproRateSup"])
    segments(jittered_conc - delta, x0[, "reproRateInf"],
             jittered_conc + delta, x0[, "reproRateInf"])
    segments(jittered_conc - delta, x0[, "reproRateSup"],
             jittered_conc + delta, x0[, "reproRateSup"])
    
    points(jittered_conc, x0$resp, pch = 20)
  }
  
  if (addlegend) {
    legend("bottomleft", pch = c(20, NA),
           lty = c(NA, 1),
           lwd = c(NA, ifelse(rect, 5, 1)),
           col = c(1, ifelse(rect, "gray80", 1)),
           legend = c("Observed values", "Confidence interval"),
           bty = "n")
  }
}

reproDoseResponseCIGG <- function(x, conc_val, jittered_conc, transf_data_conc,
                                  display.conc, main, addlegend,
                                  remove.someLabels, rect) {
  
  x0 <- cbind(x[order(x$conc),], jittered_conc = as.vector(jittered_conc))
  
  df <- data.frame(x0,
                   transf_data_conc,
                   display.conc,
                   Points = "Observed values")
  
  dfCI <- data.frame(x0,
                     transf_data_conc,
                     display.conc,
                     Conf.Int = "Confidence interval")
  
  # colors
  valCols <- fCols(df, fitcol = NA, cicol = NA)
  
  if (rect) {
    delta <- 0.01 * (max(conc_val) - min(conc_val))
    
    gf_1 <- ggplot(dfCI) + geom_rect(aes(xmin = sort(transf_data_conc) - delta,
                                         xmax = sort(transf_data_conc) + delta,
                                         ymin = reproRateInf,
                                         ymax = reproRateSup,
                                         fill = Conf.Int),
                                     alpha = 0.2) +
      geom_point(aes(x = sort(transf_data_conc), y = resp, shape = Points), df,
                 col = valCols$cols1)
  } else {
    gf_1 <- ggplot(dfCI) + geom_segment(aes(x = jittered_conc,
                                            xend = jittered_conc,
                                            y = reproRateInf,
                                            yend = reproRateSup,
                                            linetype = Conf.Int),
                                        arrow = arrow(length = unit(0.1, "cm"),
                                                      angle = 90, ends = "both"),
                                        col = valCols$cols3) +
      geom_point(aes(x = jittered_conc, y = resp, shape = Points), df,
                 col = valCols$cols1)
  }
  
  gf <- gf_1 +
    scale_shape(name = "") +
    scale_linetype(name = "") +
    scale_fill_hue("") +
    expand_limits(y = 0) +
    ggtitle(main) +
    labs(x = "Concentration", y = "Reproduction rate") +
    scale_x_continuous(breaks = unique(transf_data_conc),
                       labels = if (remove.someLabels) {
                         exclude_labels(unique(display.conc))
                       } else {
                         unique(display.conc)
                       }
    ) +
    theme_minimal()
  
  if (addlegend) {
    gf
  } else {
    gf + theme(legend.position = "none") # remove legend
  }
}
