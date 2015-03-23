#' Plot of survival data
#' 
#' The \code{survFullPlot} function plots the number of survivors as a
#' function of time for each concentration and each replicate. These function
#' is also used by \code{\link{reproDataCheck}}.
#' 
#' 
#' @param data Raw dataframe with 4 columns: "\code{replicate}", "\code{conc}",
#' "\code{time}", "\code{Nsurv}". See \code{\link{reproData}} for details.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Number of
#' survivors}.
#' @param type Graphical method: \code{generic}, \code{lattice} or
#' \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' pooled together for a same concentration by mean.
#' @note When \code{type = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot2}} and returns an object of class \code{ggplot}.
#' When \code{type = "lattice"}, the function returns an object of class
#' \code{trellis}.
#' 
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link[ggplot2]{ggplot}}, \code{\link[lattice]{xyplot}},
#' \code{\link{reproDataCheck}}
#' 
#' @keywords plot
#' 
#' @examples
#' 
#' # (1) Load the data
#' data(zinc)
#' 
#' # reproduction data
#' 
#' # (2) Plot the survival data
#' survFullPlot(zinc, type = "generic", addlegend = TRUE)
#' 
#' # (3) Plot the survival data with a lattice type
#' survFullPlot(zinc, type = "lattice", addlegend = TRUE)
#' 
#' # (4) Plot the survival data with a ggplot type
#' survFullPlot(zinc, type = "ggplot", addlegend = FALSE)
#' 
#' # (5) To build a specific legend with a ggplot type
#'  fu <- survFullPlot(zinc, type = "ggplot", addlegend = FALSE)
#'  fu + theme(legend.position = "left") + scale_colour_hue("Replicate")
#' 
#' @export
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#' @importFrom lattice lattice.options xyplot
#' 
survFullPlot <- function(data,
                         xlab,
                         ylab,
                         type = "generic",
                         addlegend = TRUE,
                         pool.replicate = FALSE) {
  # check data
  if (!is.data.frame(data))
    stop("data.frame expected !")
  
  # check data
  ref.names <- c("replicate","conc","time","Nsurv")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]
  
  if (length(missing.names) != 0)
    stop(paste("\nThe column" , missing.names,
               "is missing or wrong name in the data.frame 'data' !\n",
               sep = ""))
  rm(ref.names)
  rm(missing.names)
  
  #reorder dataset
  data <- data[order(data$replicate, data$conc, data$time), ]
  
  # default argument
  if (missing(xlab)) {
    xlab <- "Time"
  }
  if (missing(ylab)) {
    ylab <- "Number of survivors"
  }
  
  if (pool.replicate) { # pool.replicate no: 1 plot by conc
    # agregate by sum of replicate
    data <- cbind(aggregate(Nsurv ~ time + conc, data, sum),
                  replicate = 1)
  }
  
  # generic option
  if (type == "generic") {
    survFullPlotGeneric(data, xlab, ylab, addlegend)
  } else {
    # lattice option
    if (type == "lattice") {
      survFullPlotL(data, xlab, ylab, addlegend)
    } else {
      # ggplot2 option
      if (type == "ggplot") {
        survFullPlotGG(data, xlab, ylab, addlegend)
      }
    }
  }
}
