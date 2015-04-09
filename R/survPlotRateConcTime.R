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
                                   xlab,
                                   ylab,
                                   target.time = NULL,
                                   log.scale = FALSE,
                                   style = "generic",
                                   addlegend = TRUE,
                                   pool.replicate = TRUE) {
  
  # response variable
  data$response <- data$Nsurv / data$Ninit
  
  # default argument
  if (missing(xlab)) {
    xlab <- "Concentration"
  }
  if (missing(ylab)) {
    ylab <- "Survival rate"
  }
  if (missing(target.time)) {
    target.time <- max(data$time)
  }
  
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
         xlab = "Concentrations",
         ylab = "Survival rate")
    
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
      df <- ggplot(responsetable, aes(x = conc2, y = response,
                                      color = time))
    } else {
      df <- ggplot(responsetable, aes(x = conc2, y = resp,
                                      color = time, group = interaction(time,
                                                                        replicate)))
    }
    fd <- df + geom_point() + theme_minimal() +
      labs(x = "Concentrations",
           y = "Survival rate") +
      scale_color_gradient("Time", limits = range(responsetable$time),
                           high = "#132B43", low = "#56B1F7") +
      scale_x_continuous(breaks = unique(responsetable$conc2),
                         labels = unique(responsetable$conc))
    
    # legend option
    if (addlegend){
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
