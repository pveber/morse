#' Plot of survival data
#' 
#' This function plots the survival rate in function of replicates for a fixed
#' concentration and a fixed timepoint.
#' 
#' @param x an object of class \code{survData}.
#' @param target.time a numeric value corresponding to some observed time in
#' \code{x}.
#' @param concentration a numeric value corresponding to some observed concentration
#' in \code{x}.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter
survDataPlotReplicates <- function(x,
                                   target.time = NULL,
                                   concentration = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE) {
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Replicates"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(concentration)) {
    concentration <- max(x$conc)
  }
  if (is.null(target.time)) {
    target.time <- max(x$time)
  }
  
  # check [target.time] and [concentration]
  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")
  
  if (!concentration %in% x$conc)
    stop("The choosen value for [concentration] is not possible !")
  
  # select the concnetration and the target.time
  x <- filter(x, conc == concentration & time == target.time)
  
  if (style == "generic") {
    plot(factor(x$replicate), x$response,
         type = "n",
         xlab = xlab,
         ylab = ylab)
  }
  
  if (style == "ggplot") {
    df <- ggplot(x, aes(x = replicate, y = response))
    df + geom_point() + labs(x = xlab, y = ylab)
  }
}
