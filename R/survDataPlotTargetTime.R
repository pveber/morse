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
#' 
#' @export
#' @importFrom dplyr %>% filter
survDataPlotTargetTime <- function(x,
                                   target.time = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE) {
  
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
  
  # select the target.time
  x <- filter(x, x$time == target.time)
  
  # vector color
  x$color <- as.numeric(as.factor(x$replicate))
  
  if (style == "generic") {
    plot(x$conc, seq(0, max(x$response), length.out = length(x$conc)),
         type = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab)
    
    axis(side = 1, at = unique(x$conc),
         labels = unique(x$conc))
    
    # points
    if (length(unique(x$replicate)) == 1) {
      # points
        points(x$conc, x$response,
               pch = 16)
    } else {
      by(x, list(x$replicate),
         function(x) {
           points(x$conc, x$response,
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
    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = conc, y = response))
    } else {
      df <- ggplot(x, aes(x = conc, y = response,
                                      color = factor(replicate),
                                      group = replicate))
    }
    fd <- df + geom_point() + theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = unique(x$conc),
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
