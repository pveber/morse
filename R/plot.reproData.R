#' Plotting method for reproData objects
#'
#' Plots the cumulated number of offspring as a
#' function of either time and concentration, time only (for a fixed
#' concentration), concentration only (for a given target time). If both
#' concentration and target time are fixed, the function additionally plots
#' the experimental values for the minimum available concentration.
#'
#' @param x an object of class \code{reproData}
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param concentration a numeric value corresponding to some concentration in \code{data}
#' @param style graphical backend, can be \code{'generic'}, \code{'lattice'} or
#' \code{'ggplot'}
#' @param addlegend if \code{TRUE}, a default legend is added to the plot
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param \dots further arguments to be passed to generic methods (xlab, ylab, ...).
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot2}} and returns an object of class \code{ggplot}.
#' When \code{style = "lattice"}, the function returns an object of class
#' \code{trellis}.
#'
#'#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(zinc)
#' zinc <- reproData(zinc)
#'
#' # (2) Plot the reproduction data
#' plot(zinc, style = "generic", addlegend = TRUE)
#'
#' # (3) Plot the survival data with a lattice type
#' plot(zinc, style = "lattice", addlegend = TRUE)
#'
#' # (4) Plot the survival data with a ggplot type
#' plot(zinc, style = "ggplot", addlegend = FALSE)
#'
#' # (5) To build a specific legend with a ggplot type
#' fu <- plot(zinc, style = "ggplot", addlegend = FALSE)
#' fu + theme(legend.position = "left") + scale_colour_hue("Replicate")
#'
#' # (6) Plot the reproduction rate in function of replicates for one concentration at
#' # one target.time with a generic type
#' plot(zinc, style = "generic", target.time = 21, concentration = 0.66)
#'
#' # (7) Plot the reproduction rate in function of replicates for one concentration at
#' # one target.time with a ggplot type
#' plot(zinc, style = "ggplot", target.time = 21, concentration = 0.66)

#' @export
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
plot.reproData <- function(x,
                           target.time = NULL,
                           concentration = NULL,
                           style = "generic",
                           pool.replicate = FALSE,
                           ...) {
  if(! is(x, "reproData"))
    stop("plot.reproData: object of class reproData expected")
  
  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(Nreprocumul ~ time + conc, x, sum),
               replicate = 1)
  }
  
  if (! is.null(target.time) && is.null(concentration))
    reproDataPlotTargetTime(x, target.time, style, addlegend, ...)
  else if (is.null(target.time) && ! is.null(concentration))
    reproDataPlotFixedConc(x, concentration, style, addlegend, ...)
}
