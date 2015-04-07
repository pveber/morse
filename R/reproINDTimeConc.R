#' Plot of reproduction data
#' 
#' This function plots the reproduction rate in function of time per concentration.
#' 
#' @param data an object of class \code{reproData}.
#' @param FirstTimeClutch the first day when individuals are mature.
#' @param cumul If \code{TRUE}, the response is in cumulative number of offspring
#' by individuals-day.
#' @param type Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' pooled together for a same concentration by mean.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr %>%
#' @importFrom plyr ldply
reproINDTimeConc <- function(data, FirstTimeClutch = NULL,
                                cumul = TRUE, type = "generic",
                                pool.replicate = FALSE) {
  
  if (cumul) { # response cumulate of Nrepro and NID
    data$response <- data$Nreprocumul / data$Nindtime
    data[is.na(data)] <- 0
    ylab <- "Cumulated Nb of offspring / NID"
  } else {
    data$response <- data$Nrepro / data$NindtimenoncumulTi
    ylab <- "Nb of offsring / Delta NID"
    data[is.na(data)] <- 0
  }
  # pool replicate mean of replicate
  if (pool.replicate) {
    responsetable <- aggregate(data$response, by = list(data$time, data$timeOrigine,
                                                        data$conc),
                               function(x) { mean(x, na.rm = T) })
    colnames(responsetable) <- c("time", "timeOrigine", "conc",
                                 "response")
  } else {
    responsetable <- data
  }
  
  
  responsetable$color <- as.numeric(as.factor(responsetable$conc))
  
  if (is.null(FirstTimeClutch)) {
    FirstTimeClutch <- 0
  }
  
  if (type == "generic") {
    plot(responsetable$timeOrigine, seq(0, max(sapply(responsetable$response, max)),
                                        length.out = length(responsetable$timeOrigine)),
         type = "n",
         xlab = "Time",
         ylab = ylab)
    
    # lines
    if (pool.replicate) {
      by(responsetable, responsetable$conc, function(x) {
        lines(x$timeOrigine, x$response,
              col = x$color)
        points(x$timeOrigine, x$response,
               pch = 16,
               col = x$color)
      })
    } else {
      by(responsetable, list(responsetable$replicate, responsetable$conc), function(x) {
        lines(x$timeOrigine, x$response,
              col = x$color)
        points(x$timeOrigine, x$response,
               pch = 16,
               col = x$color)
      })
    }
    legend("topright", legend = unique(responsetable$conc) ,
           col = unique(responsetable$color),
           pch = 16,
           lty = 1)
  }
  if (type == "ggplot") {
    require("ggplot2")
    if (pool.replicate) {
      df <- ggplot(responsetable, aes(x = timeOrigine, y = response,
                                      color = factor(conc)))
    } else {
      df <- ggplot(responsetable, aes(x = timeOrigine, y = response,
                                      color = factor(conc),
                                      group = interaction(conc, replicate)))
    }
    
    df + geom_line() + geom_point() + theme_minimal() +
      labs(x = "Time",
           y = ylab) +
      scale_color_hue("Concentrations")
  }
}