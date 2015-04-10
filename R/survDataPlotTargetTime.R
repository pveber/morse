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
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE,
                                   pool.replicate = FALSE,
                                   log.scale = FALSE) {
  
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
  
  #pool replicate
  if (pool.replicate) {
    responsetable <- aggregate(data$response, by = list(data$conc,
                                                        data$conc2), mean)
    colnames(responsetable) <- c("conc", "conc2", "response")
  } else {
    responsetable <- data
  }
  
  # vector color
  responsetable$color <- if (pool.replicate) {
    1
  } else {
    as.numeric(as.factor(responsetable$replicate))
  }
  
  if (style == "generic") {
    plot(responsetable$conc2, seq(0, 1, length.out = length(responsetable$conc2)),
         type = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab)
    
    axis(side = 1, at = unique(responsetable$conc2),
         labels = unique(responsetable$conc))
    
    # points
    if (pool.replicate) {
      # points
        points(responsetable$conc2, responsetable$response,
               pch = 16)
    } else {
      by(responsetable, list(responsetable$replicate),
         function(x) {
           points(x$conc2, x$response,
                  pch = 16,
                  col = x$color)
         })
      if (addlegend) {
        legend("bottomleft", legend = unique(responsetable$replicate) ,
               title = "Replicate",
               col = unique(responsetable$color),
               pch = 16, ncol = 3)
      }
    }
  }
  if (style == "ggplot") {
    if (pool.replicate) {
      df <- ggplot(responsetable, aes(x = conc2, y = response))
    } else {
      df <- ggplot(responsetable, aes(x = conc2, y = response,
                                      color = factor(replicate),
                                      group = replicate))
    }
    fd <- df + geom_point() + theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = unique(responsetable$conc2),
                         labels = unique(responsetable$conc)) +
      scale_color_hue("Replicate")
    
    # legend option
    if (addlegend) {
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
