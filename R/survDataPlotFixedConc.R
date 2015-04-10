#' Plot of survival data
#' 
#' This function plots the survival rate in function of time per concentration.
#' 
#' @param data an object of class \code{survData}.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' pooled together for a same concentration by mean.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
survDataPlotFixedConc <- function(data,
                                  concentration = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  style = "generic",
                                  addlegend = TRUE,
                                  pool.replicate = FALSE) {
  
  
  # response variable
  data$response <- data$Nsurv / data$Ninit
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Time"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(concentration)) {
    concentration <- max(data$conc)
  }
  
  # check concentration value
  if (!concentration %in% data$conc)
    stop("The choosen value for [concentration] is not possible !")
  
  # select the concentration
  data <- filter(data, data$conc == concentration)
  
  #pool replicate
  if (pool.replicate) {
    responsetable <- aggregate(data$response, by = list(data$time),
                               mean)
    colnames(responsetable) <- c("time", "response")
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
    plot(responsetable$time, seq(0, 1, length.out = length(responsetable$time)),
         type = "n",
         xlab = xlab,
         ylab = ylab)
    
    if (pool.replicate) {
      lines(responsetable$time, responsetable$response) # lines
      points(responsetable$time, responsetable$response, # points 
             pch = 16)
    } else {
      # one line by replicate
      by(responsetable, list(responsetable$replicate),
         function(x) {
           lines(x$time, x$response, # lines
                 col = x$color)
           points(x$time, x$response, # points
                  pch = 16,
                  col = x$color)
         })
      
      if (addlegend) { # only if pool.replicate == FALSE
        legend("bottomleft", legend = unique(responsetable$replicate) ,
               col = unique(responsetable$color),
               pch = 16,
               lty = 1)
      }
    }
  }
  if (style == "ggplot") {
    if (pool.replicate) {
      df <- ggplot(responsetable, aes(x = time, y = response))
    } else {
      df <- ggplot(responsetable, aes(x = time, y = response,
                                      color = factor(replicate),
                                      group = replicate))
    }
    fd <- df + geom_line() + geom_point() + theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_color_hue("Replicate")
    
    if (addlegend) {# only if pool.replicate == FALSE
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
