#' Plot of survival data
#' 
#' This function plots the survival rate in function of concentration per timepoints.
#' 
#' @param data an object of class \code{survData}.
#' @param target.time a numeric value corresponding to some observed time in
#' \code{data}.
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
survDataPlotTargetTime <- function(data,
                                   target.time = NULL,
                                   log.scale = FALSE,
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE,
                                   pool.replicate) {
  
  # response variable
  data$response <- data$Nsurv / data$Ninit
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Concentration"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(target.time)) {
    target.time <- max(data$time)
  }
  if (!target.time %in% data$time)
    stop("[target.time] is not one of the possible time !")
  
  # select only no null concentration datapoint for log reprsentation 
  sel <- if (log.scale) {
    data$conc > 0
  } else {
    rep(TRUE, length(data$conc))
  }
  if (log.scale) {
    data$conc2[sel] <- log(data$conc[sel])
  } else {
    data$conc2[sel] <- data$conc[sel]
  }
  data <- na.omit(data)
  
  # select the target.time
  data <- filter(data, data$time == target.time)
  
  # vector color
  data$color <- if (pool.replicate) {
    1
  } else {
    as.numeric(as.factor(data$replicate))
  }
  
  if (style == "generic") {
    plot(data$conc2, seq(0, 1, length.out = length(data$conc2)),
         type = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab)
    
    axis(side = 1, at = unique(data$conc2),
         labels = unique(data$conc))
    
    # points
    if (pool.replicate) {
      # points
        points(data$conc2, data$response,
               pch = 16)
    } else {
      by(data, list(data$replicate),
         function(x) {
           points(x$conc2, x$response,
                  pch = 16,
                  col = x$color)
         })
      if (addlegend) {
        legend("bottomleft", legend = unique(data$replicate) ,
               title = "Replicate",
               col = unique(data$color),
               pch = 16, ncol = 3)
      }
    }
  }
  if (style == "ggplot") {
    if (pool.replicate) {
      df <- ggplot(data, aes(x = conc2, y = response))
    } else {
      df <- ggplot(data, aes(x = conc2, y = response,
                                      color = factor(replicate),
                                      group = replicate))
    }
    fd <- df + geom_point() + theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = unique(data$conc2),
                         labels = unique(data$conc)) +
      scale_color_hue("Replicate")
    
    # legend option
    if (addlegend) {
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
