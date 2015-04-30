survDataPlotFullGeneric <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate
  # for generic graphics

  .convert <- function(x) {
    # conversion of a replicate name in a number coding for color
    # INPUT
    # - x: name of replicate
    # OUTPUT
    #  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
    # - a position of replicate in vector to display color
    replicate <- unique(data$replicate)
    mat <- matrix(ncol = 2, nrow = length(replicate))
    mat[, 1] <- as.character(replicate)
    mat[, 2] <- as.character(seq(1, length(replicate)))
    return(as.integer(mat[mat[, 1] == as.character(x)][2]))
  }

  .NbPlot <- function(conc) {
    # definition of the number of subplots
    # INPUT
    # - conc: vector of tested concentrations
    # OUTPUT
    # - vector defining the number of columns and rows in par(mfrow)

    nbconc <- length(conc)
    PlotPar <- c(c(2, 2), c(2, 3), c(2, 4), c(3, 3), c(2, 5), c(3, 4), c(3, 5),
                 c(4, 4))
    NbPlotTheo <- matrix(ncol = 2, nrow = 8)
    NbPlotTheo[, 1] <- c(1, 3, 5, 7, 9, 11, 13, 15)
    NbPlotTheo[, 2] <- c(4, 6, 8, 9, 10, 12, 15, 16)
    if (nbconc < 15) {
      PlotTrue <- NbPlotTheo[NbPlotTheo[, 2] - nbconc > 0, ][1, 1]
    } else {
      PlotTrue = 15
    }
    return(c(PlotPar[PlotTrue], PlotPar[PlotTrue + 1]))
  }


  # creation of a vector of colors
  colors <- rainbow(length(unique(data$replicate)))
  pchs <- as.numeric(unique(data$replicate))
  # split of the graphical window in subplots
  par(mfrow = .NbPlot(unique(data$conc)))

  by(data, data$conc, function(x) {
    # bakground
    plot(x$time, rep(0, length(x$time)),
         xlab = xlab,
         ylab = ylab,
         ylim = c(0,max(x$Nsurv)),
         type = "n",
         col = 'white',
         yaxt = 'n')

    # axis
    axis(side = 2, at = pretty(c(0,max(x$Nsurv))))

    # lines and points
    by(x, x$replicate, function(y) {
      lines(y$time, y$Nsurv,
            type = "l",
            col = colors[.convert(unique(y$replicate))])
      points(y$time, y$Nsurv,
             pch = pchs[.convert(unique(y$replicate))],
             col = colors[.convert(unique(y$replicate))])
    })

    # title
    title(paste("Conc: ", unique(x$conc), sep = ""))
  })

  if (addlegend) {
    # creation of an empty plot to display legend
    plot(0, 0,
         xlab = "",
         ylab = "",
         xlim = c(0,5),
         ylim = c(0,5),
         type = "n",
         xaxt = "n",
         yaxt = "n")

    # Display legend
    title.legend <- "Replicate"
    mat <- matrix(nrow = length(unique(data$replicate)), ncol = 2)
    mat[, 1] <- rep(title.legend, length(unique(data$replicate)))
    mat[, 2] <- unique(as.character(data$replicate))
    name <- apply(mat, 1, function(x) {paste(x[1], x[2], sep = ": ")})

    legend("top", name,
           lty = rep(1, length(unique(data$replicate))),
           pch = pchs,
           col = colors,
           bty = "n",
           cex = 1)
  }
  par(mfrow = c(1, 1))
}


#' @importFrom lattice lattice.options xyplot
survDataPlotFullLattice <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate for lattice graphics

  # change order of reading concentrations
  lattice.options(default.args = list(as.table = TRUE))

  if (addlegend) {
    xyplot(Nsurv ~ time | factor(conc),
           data = data,
           group = replicate,
           index.cond = list(c((round(length(unique(data$conc)) / 2) + 1):length(unique(data$conc)),
                               1:(round(length(unique(data$conc)) / 2)))),
           type = "b",
           pch = 16,
           xlab = xlab,
           ylab = ylab,
           auto.key = list(space = "right",
                           title = "Replicate",
                           cex.title = 1,
                           lines = TRUE,
                           points = FALSE))
  } else {
    xyplot(Nsurv ~ time | factor(conc),
           data = data,
           group = replicate,
           index.cond = list(c((round(length(unique(data$conc)) / 2) + 1):length(unique(data$conc)),
                               1:(round(length(unique(data$conc))/2)))),
           type = "b",
           pch = 16,
           xlab = xlab,
           ylab = ylab)
  }
}

survDataPlotFullGG <- function(data, xlab, ylab, addlegend) {
  # plot of survival data: one subplot for each concentration, and one color for
  # each replicate for ggplot graphics

  time = NULL
  Nsurv = NULL
  title.legend <- "Replicate"

  # create ggplot object Nsurv / time / replicate / conc
  fg <- ggplot(data, aes(time, Nsurv, colour = factor(replicate))) +
    geom_point() +
    geom_line() +
    labs(x = xlab, y = ylab) +
    facet_wrap(~conc, nrow = 2) +
    scale_x_continuous(breaks = unique(data$time)) +
    ylim(0, max(data$Nsurv))

  # legend option
  if (addlegend){
    fd <- fg + scale_colour_hue(title.legend) # the default legend
  } else {
    fd <- fg + theme(legend.position = "none") # remove legend
  }
  return(fd)
}


survDataPlotFull <- function(data,
                             xlab = NULL,
                             ylab = NULL,
                             style = "generic",
                             addlegend = TRUE) {
  if (style == "generic")
    survDataPlotFullGeneric(data, xlab, ylab, addlegend)
  else if (style == "lattice")
    survDataPlotFullLattice(data, xlab, ylab, addlegend)
  else if (style == "ggplot")
    survDataPlotFullGG(data, xlab, ylab, addlegend)
  else stop("Unknown plot style")
}

#' Plotting method for survData objects
#'
#' Plots the number of survivors as a
#' function of either time and concentration, time only (for a fixed
#' concentration), concentration only (for a given target time). If both
#' concentration and target time are fixed, the function additionally plots
#' the experimental values for the minimum available concentration.
#'
#' @param x an object of class \code{survData}
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param concentration a numeric value corresponding to some concentration in \code{data}
#' @param xlab a label for the \eqn{X}-axis, by default \code{Time}
#' @param ylab a label for the \eqn{Y}-axis, by default \code{Number of
#' survivors}.
#' @param style graphical backend, can be \code{'generic'}, \code{'lattice'} or
#' \code{'ggplot'}
#' @param addlegend if \code{TRUE}, a default legend is added to the plot
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param \dots further arguments to be passed to generic methods.
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
#' zinc <- survData(zinc)
#'
#' # (2) Plot the survival data
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

#' @export
#'
#' @import grDevices
# FIXME: delete imports if really not needed
# @importFrom gridExtra grid.arrange arrangeGrob
# @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#'
plot.survData <- function(x,
                          target.time = NULL,
                          concentration = NULL,
                          xlab = NULL,
                          ylab = NULL,
                          style = "generic",
                          addlegend = TRUE,
                          pool.replicate = FALSE,
                          ...) {

  if(! is(x,"survData"))
    stop("plot.survData: object of class survData expected")

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(Nsurv ~ time + conc, x, sum),
                  replicate = 1)
  }

  if (is.null(target.time) && is.null(concentration))
    survDataPlotFull(x, xlab, ylab, style, addlegend)
  else if (! is.null(target.time) && is.null(concentration))
    survDataPlotTargetTime(x, target.time, xlab, ylab, style, addlegend)
  else if (is.null(target.time) && ! is.null(concentration))
    survDataPlotFixedConc(x, concentration, xlab, ylab, style, addlegend)
  else
    survDataPlotReplicates(x, target.time, concentration, xlab, ylab, style, addlegend)
}
