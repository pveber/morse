reproDataPlotFull <- function(data, style = "generic", addlegend = TRUE, ...) {
  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Time"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Cumulated Number of offsprings"
  dataPlotFull(data, "Nreprocumul", xlab, ylab, style, addlegend, ...)
}


#' @import ggplot2
#' @importFrom dplyr filter
reproDataPlotTargetTime <- function(x,
                                    target.time,
                                    style,
                                    addlegend, ...) {
  # plot of cumulated number of offspring as a funtion of concentration
  # for a fixed time

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Concentration"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Cumulated Number of offsprings"

  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")

  # select the target.time
  x <- filter(x, x$time == target.time)

  # Define visual parameters
  mortality <- c(0, 1) # code 0/1 mortality
  nomortality <- match(x$Nsurv == x$Ninit, c(TRUE, FALSE))

  # without mortality
  mortality <- mortality[nomortality] # vector of 0 and 1

  # encodes mortality empty dots (1) and not mortality solid dots (19)
  if (style == "generic") {
    mortality[which(mortality == 0)] <- 19
  }
  if (style == "ggplot") {
    mortality[which(mortality == 0)] <- "No"
    mortality[which(mortality == 1)] <- "Yes"
  }

  # default legend argument
  legend.position <- "right"
  legend.title <- "Mortality"
  legend.name.no <- "No"
  legend.name.yes <- "Yes"

  # generic
  if (style == "generic") {
    plot(x$conc,
         x$Nreprocumul,
         xlab = xlab,
         ylab = ylab,
         pch = mortality,
         yaxt = "n",
         xaxt = "n",
         ...)
    # axis
    axis(side = 2, at = pretty(c(0, max(x$Nreprocumul))))
    axis(side = 1, at = unique(x$conc), labels = unique(x$conc))

    # legend
    if (addlegend) {
      legend(legend.position,title = legend.title, pch = c(19, 1), bty = "n",
             legend = c(legend.name.no, legend.name.yes))
    }
  }

  #ggplot2
  if (style == "ggplot") {
    df <- data.frame(x, mortality = mortality)

    # plot
    gp <- ggplot(df, aes(conc, Nreprocumul, colour = factor(mortality))) +
      geom_point(size = 2.5) +
      labs(x = xlab, y = ylab) +
      theme_minimal() +
      scale_colour_hue(legend.title, breaks = c("No","Yes"),
                       labels = c(legend.name.no, legend.name.yes)) +
      scale_x_continuous(breaks = unique(df$conc))
    if (addlegend) {
      return(gp)
    } else {
      gp <- gp + theme(legend.position = "none")
      return(gp)
    }
  }
}

reproDataPlotFixedConc <- function(x,
                                   concentration,
                                   style = "generic",
                                   addlegend = TRUE,
                                   ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Time"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Cumulated Number of offsprings"
  dataPlotFixedConc(x, "Nreprocumul", concentration, xlab, ylab, style, addlegend, ...)
}

reproDataPlotReplicates <- function(x,
                                   target.time,
                                   concentration,
                                   style,
                                   addlegend,
                                   ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Replicate"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Cumulated Number of offsprings"
  dataPlotReplicates(x, "Nreprocumul", target.time, concentration, xlab, ylab, style, addlegend)
}

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
#' @param style graphical backend, can be \code{'generic'} or
#' \code{'ggplot'}
#' @param addlegend if \code{TRUE}, a default legend is added to the plot
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param \dots further arguments to be passed to generic methods (xlab, ylab, ...).
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot2}} and returns an object of class \code{ggplot}.
#' @keywords plot
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(cadmium1)
#' cadmium1 <- reproData(cadmium1)
#'
#' # (2) Plot the reproduction data for a fixed time
#' plot(cadmium1, target.time = 21, style = "generic", addlegend = TRUE)
#'
#' # (3) Plot the reproduction data for a fixed time with a ggplot type
#' plot(cadmium1, target.time = 21, style = "ggplot", addlegend = FALSE)
#'
#' # (4) Plot the reproduction data for a fixed concentration
#' plot(cadmium1, concentration = 4.36, style = "generic", addlegend = TRUE)
#'
#' # (5) Plot the reproduction data for a fixed concentration with a ggplot type
#' plot(cadmium1, concentration = 4.36, style = "ggplot", addlegend = FALSE)
#'
#' # (6) Plot the reproduction data in function of replicates for one concentration at
#' # one target.time with a generic type
#' plot(cadmium1, style = "generic", target.time = 21, concentration = 0.86)
#'
#' # (7) Plot the reproduction data in function of replicates for one concentration at
#' # one target.time with a ggplot type
#' plot(cadmium1, style = "ggplot", target.time = 21, concentration = 0.86)
#'
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
                           addlegend = TRUE,
                           pool.replicate = FALSE,
                           ...) {
  if(! is(x, "reproData"))
    stop("plot.reproData: object of class reproData expected")

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(Nreprocumul ~ time + conc, x, sum),
               replicate = 1)
  }

  if (is.null(target.time) && is.null(concentration))
    reproDataPlotFull(x, style, addlegend, ...)
  else if (! is.null(target.time) && is.null(concentration))
    reproDataPlotTargetTime(x, target.time, style, addlegend, ...)
  else if (is.null(target.time) && ! is.null(concentration))
    reproDataPlotFixedConc(x, concentration, style, addlegend, ...)
  else
    reproDataPlotReplicates(x, target.time, concentration, style,
                            addlegend, ...)
}
