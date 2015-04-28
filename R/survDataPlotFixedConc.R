#' Plot of survival data
#' 
#' This function plots the survival rate in function of time per concentration.
#' 
#' @param x an object of class \code{survData}.
#' @param concentration a numeric value corresponding to some observed concentration
#' in \code{x}.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
survDataPlotFixedConc <- function(x,
                                  concentration = NULL,
                                  xlab = NULL,
                                  ylab = NULL,
                                  style = "generic",
                                  addlegend = TRUE) {
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Time"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(concentration)) {
    concentration <- max(x$conc)
  }
  
  # check concentration value
  if (!concentration %in% x$conc)
    stop("The choosen value for [concentration] is not possible !")
  
  # select the concentration
  x <- filter(x, x$conc == concentration)
  
  # vector color
  x$color <- as.numeric(as.factor(x$replicate))
  
  if (style == "generic") {
    plot(x$time, seq(0, 1, length.out = length(x$time)),
         type = "n",
         xlab = xlab,
         ylab = ylab)
    
    # one line by replicate
    by(x, list(x$replicate),
       function(x) {
         lines(x$time, x$response, # lines
               col = x$color)
         points(x$time, x$response, # points
                pch = 16,
                col = x$color)
         })
      
      if (addlegend) {
        legend("bottomleft", legend = unique(x$replicate) ,
               col = unique(x$color),
               pch = 16,
               lty = 1)
      }
    }
  if (style == "ggplot") {
    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = time, y = response))
    } else {
      df <- ggplot(x, aes(x = time, y = response,
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
