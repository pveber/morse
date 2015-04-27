#' Plot of survival data
#' 
#' This function plots the survival rate in function of concentration per timepoints.
#' 
#' @param x an object of class \code{survData}.
#' @param target.time a numeric value corresponding to some observed time in
#' \code{x}.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param log.scale If \code{TRUE}, a log-scale is used on the \eqn{X}-axis.
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' pooled together for a same concentration by mean.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter
survDataPlotTargetTime <- function(x,
                                   target.time = NULL,
                                   log.scale = FALSE,
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE,
                                   pool.replicate) {
  
  # response variable
  x$response <- x$Nsurv / x$Ninit
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Concentration"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(target.time)) {
    target.time <- max(x$time)
  }
  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")
  
  # select only no null concentration datapoint for log reprsentation 
  sel <- if (log.scale) {
    x$conc > 0
  } else {
    rep(TRUE, length(x$conc))
  }
  if (log.scale) {
    x$conc2[sel] <- log(x$conc[sel])
  } else {
    x$conc2[sel] <- x$conc[sel]
  }
  x <- na.omit(x)
  
  # select the target.time
  x <- filter(x, x$time == target.time)
  
  # vector color
  x$color <- if (pool.replicate) {
    1
  } else {
    as.numeric(as.factor(x$replicate))
  }
  
  if (style == "generic") {
    plot(x$conc2, seq(0, 1, length.out = length(x$conc2)),
         type = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab)
    
    axis(side = 1, at = unique(x$conc2),
         labels = unique(x$conc))
    
    # points
    if (pool.replicate) {
      # points
        points(x$conc2, x$response,
               pch = 16)
    } else {
      by(x, list(x$replicate),
         function(x) {
           points(x$conc2, x$response,
                  pch = 16,
                  col = x$color)
         })
      if (addlegend) {
        legend("bottomleft", legend = unique(x$replicate) ,
               title = "Replicate",
               col = unique(x$color),
               pch = 16, ncol = 3)
      }
    }
  }
  if (style == "ggplot") {
    if (pool.replicate) {
      df <- ggplot(x, aes(x = conc2, y = response))
    } else {
      df <- ggplot(x, aes(x = conc2, y = response,
                                      color = factor(replicate),
                                      group = replicate))
    }
    fd <- df + geom_point() + theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = unique(x$conc2),
                         labels = unique(x$conc)) +
      scale_color_hue("Replicate")
    
    # legend option
    if (addlegend) {
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
