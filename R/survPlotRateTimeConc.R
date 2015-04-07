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
#' @importFrom dplyr %>%
#' @importFrom plyr ldply
survPlotRateTimeConc <- function(data, type = "generic",
                                 addlegend = TRUE,
                                 pool.replicate = FALSE) {
  
  data$resp <- data$Nsurv / data$Ninit
  
  #pool replicate
  if (pool.replicate) {
    responsetable <- aggregate(data$resp, by = list(data$time, data$conc), mean)
    colnames(responsetable) <- c("time", "conc", "response")
  } else {
    responsetable <- data
  }
  
  responsetable$color <- as.numeric(as.factor(responsetable$conc))
  
  if (type == "generic") {
    plot(responsetable$time, seq(0, 1, length.out = length(responsetable$time)),
         type = "n",
         xlab = "Time",
         ylab = "Survival rate")
    
    # lines
    if (pool.replicate) {
      by(responsetable, responsetable$conc, function(x) {
        lines(x$time, x$response,
              col = x$color)
        points(x$time, x$response,
               pch = 16,
               col = x$color)
      })
    } else {
      by(responsetable, list(responsetable$replicate, responsetable$conc),
         function(x) {
           lines(x$time, x$resp,
                 col = x$color)
           points(x$time, x$resp,
                  pch = 16,
                  col = x$color)
         })
    }
    
    if (addlegend) {
      legend("bottomleft", legend = unique(responsetable$conc) ,
             col = unique(responsetable$color),
             pch = 16,
             lty = 1)
    }
  }
  if (type == "ggplot") {
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
