#' Plot of survival data
#' 
#' This function plots the survival rate in function of time per concentration.
#' 
#' @param data an object of class \code{survData}.
#' @param type Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' pooled together for a same concentration by mean.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter
#' @importFrom plyr ldply
survPlotRateTimeConc <- function(data, concentration,
                                 xlab,
                                 ylab,
                                 style = "generic",
                                 addlegend = TRUE,
                                 pool.replicate = FALSE) {
  
  
  # response variable
  data$response <- data$Nsurv / data$Ninit
  
  # default argument
  if (missing(xlab)) {
    xlab <- "Time"
  }
  if (missing(ylab)) {
    ylab <- "Survival rate"
  }
  if (missing(concentration)) {
    concentration <- max(data$conc)
  }
  
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
  
  ifelse(test = pool.replicate, yes = 1,
         no =  'bli')
  
  if (style == "generic") {
    plot(responsetable$time, seq(0, 1, length.out = length(responsetable$time)),
         type = "n",
         xlab = xlab,
         ylab = ylab)
    
    if (pool.replicate) {
      # lines
      lines(responsetable$time, responsetable$response)
      # points
      points(responsetable$time, responsetable$response,
             pch = 16)
    } else {
      by(responsetable, list(responsetable$replicate),
         function(x) {
           # lines
           lines(x$time, x$response,
                 col = x$color)
           # points
           points(x$time, x$response,
                  pch = 16,
                  col = x$color)
         })
      if (addlegend) {
        legend("bottomleft", legend = unique(responsetable$replicate) ,
               col = unique(responsetable$color),
               pch = 16,
               lty = 1)
      }
    }
  }
  if (style == "ggplot") {
    if (pool.replicate) {
      df <- ggplot(responsetable, aes(x = time, y = response,
                                      color = factor(conc)))
    } else {
      df <- ggplot(responsetable, aes(x = time, y = resp,
                                      color = factor(conc),
                                      group = interaction(conc, replicate)))
    }
    fd <- df + geom_line() + geom_point() + theme_minimal() +
      labs(x = "Time",
           y = "Survival rate") +
      scale_color_hue("Concentrations")
    
    # legend option
    if (addlegend){
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
}
