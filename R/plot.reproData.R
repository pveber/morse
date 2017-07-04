#' Plotting method for \code{reproData} objects
#'
#' This is the generic \code{plot} S3 method for the \code{reproData} class.
#' It plots the cumulated number of offspring as a function of time.
#'
#' @param x an object of class \code{reproData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a title for the \eqn{y}-axis
#' @param main main title for the plot
#' @param concentration a numeric value corresponding to some concentration in
#' \code{data}. If \code{concentration = NULL}, draws a plot for each concentration
#' @param style graphical backend, can be \code{'ggplot'} or \code{'generic'}
#' @param pool.replicate if \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid the label overlap
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @note When \code{style = "generic"}, the function calls the generic function
#' \code{\link[graphics]{plot}}
#' @note When \code{style = "ggplot"}, the function return an object of class
#'  \code{gg} and \code{ggplot}, see function \code{\link[ggplot2]{ggplot}} 
#' 
#' @keywords plot
#'
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' cadmium1 <- reproData(cadmium1)
#'
#' # (2) Plot the reproduction data
#' plot(cadmium1)
#'
#' # (3) Plot the reproduction data for a fixed concentration
#' plot(cadmium1, concentration = 4.36, style = "generic")
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom methods is
#' @importFrom stats aggregate
#' 
#' @export
plot.reproData <- function(x,
                           xlab,
                           ylab = "Cumulated Number of offspring",
                           main = NULL,
                           concentration = NULL,
                           style = "ggplot",
                           pool.replicate = FALSE,
                           addlegend = FALSE,
                           remove.someLabels = FALSE, ...) {
  if (!is(x, "reproData"))
    stop("plot.reproData: object of class reproData expected")
  
  if (style == "generic" && remove.someLabels)
    warning("'remove.someLabels' argument is valid only in 'ggplot' style.",
            call. = FALSE)
  
  if (is.null(concentration) && addlegend)
    warning("'addlegend' argument is valid only when 'concentration' is not null.",
            call. = FALSE)
  
  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(cbind(Nreprocumul, Nsurv, Ninit) ~ time + conc, x, sum),
               replicate = 1)
  }
  
  if (is.null(concentration)) {
    reproDataPlotFull(x, xlab, ylab, style, remove.someLabels)
  }
  else {
    reproDataPlotFixedConc(x, xlab, ylab, main, concentration, style, addlegend,
                           remove.someLabels)
  }
}

reproDataPlotFull <- function(data, xlab, ylab, style = "generic",
                              remove.someLabels) {
  dataPlotFull(data, xlab, ylab, "Nreprocumul", style,
               remove.someLabels)
}


reproDataPlotFixedConc <- function(x,
                                   xlab,
                                   ylab,
                                   main,
                                   concentration,
                                   style = "generic",
                                   addlegend = FALSE,
                                   remove.someLabels = FALSE) {
  dataPlotFixedConc(x, xlab, ylab, main, "Nreprocumul",
                    concentration, style, addlegend, remove.someLabels)
}

