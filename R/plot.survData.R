#' Plotting method for \code{survData} objects
#'
#' Plots the number of survivors as a function of time (for a fixed
#' concentration).
#'
#' @param x an object of class \code{survData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a label for the \eqn{y}-axis
#' @param main main title for the plot
#' @param concentration a numeric value corresponding to some concentration in \code{data}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param pool.replicate if \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid label overlap
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot2}} and returns an object of class \code{ggplot}.
#'
#' @keywords plot
#'
#' @examples
#'
#' library(ggplot2)
#'
#' # (1) Load the data
#' data(zinc)
#' zinc <- survData(zinc)
#'
#' # (2) Plot survival data
#' plot(zinc, addlegend = TRUE)
#'
#' # (3) Plot survival data with a ggplot style
#' plot(zinc, style = "ggplot")
#'
#' # (4) To build a specific legend with a ggplot type
#' fu <- plot(zinc, style = "ggplot", addlegend = FALSE)
#' fu + theme(legend.position = "left") + scale_colour_hue("Replicate")
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot axis legend lines par points polygon
#' segments title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export
plot.survData <- function(x,
                          xlab,
                          ylab = "Number of surviving individuals",
                          main = NULL,
                          concentration = NULL,
                          style = "generic",
                          pool.replicate = FALSE,
                          log.scale = FALSE,
                          addlegend = FALSE,
                          remove.someLabels = FALSE, ...) {

  if(! is(x,"survData"))
    stop("plot.survData: object of class survData expected")

  if (pool.replicate) {
    # agregate by sum of replicate
    x <- cbind(aggregate(cbind(Nsurv, Ninit) ~ time + conc, x, sum),
               replicate = 1)
  }

  if (is.null(concentration)) {
    survDataPlotFull(x, xlab, ylab, style, remove.someLabels)
  }
  else {
    survDataPlotFixedConc(x, xlab, ylab, main, concentration,
                          style, addlegend, remove.someLabels)
  }
}


# [ReplicateIndex(data)] builds a list of indices, each one named after
# a replicate of [data], thus providing a dictionary from replicate names to
# integer keys.
ReplicateIndex <- function(data) {
  replicate <- unique(data$replicate)
  r <- as.list(seq(1, length(replicate)))
  names(r) <- as.character(replicate)
  return(r)
}


# [plotMatrixGeometry(n)] returns a vector [c(w,h)] such that a matrix of plots
# of dimension ([w], [h]) is big enough to display [n] plots in a pretty way.
# This will typically be used in [par(mfrow)] calls.
plotMatrixGeometry <- function(nblevels) {
  PlotPar <- c(c(2, 2), c(2, 3), c(2, 4), c(3, 3), c(2, 5), c(3, 4), c(3, 5),
               c(4, 4))
  NbPlotTheo <- matrix(ncol = 2, nrow = 8)
  NbPlotTheo[, 1] <- c(1, 3, 5, 7, 9, 11, 13, 15)
  NbPlotTheo[, 2] <- c(4, 6, 8, 9, 10, 12, 15, 16)
  if (nblevels < 15) {
    i <- NbPlotTheo[NbPlotTheo[, 2] - nblevels > 0, 1][1]
  } else {
    i <- 15
  }
  return(c(PlotPar[i], PlotPar[i + 1]))
}

# General full plot: one subplot for each concentration, and one color for
# each replicate (for generic graphics)
dataPlotFullGeneric <- function(data, xlab, ylab, resp) {
  replicate.index <- ReplicateIndex(data)

  # creation of a vector of colors
  colors <- rainbow(length(unique(data$replicate)))
  pchs <- as.numeric(unique(data$replicate))
  # split of the graphical window in subplots
  par(mfrow = plotMatrixGeometry(length(unique(data$conc))))

  by(data, data$conc, function(x) {
    # bakground
    plot(x$time, rep(0, length(x$time)),
         xlab = xlab,
         ylab = ylab,
         ylim = c(0, max(x[, resp])),
         type = "n",
         col = 'white',
         xaxt = "n",
         yaxt = "n")

    # axis
    axis(side = 1, at = sort(unique(x[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))

    # lines and points
    by(x, x$replicate, function(y) {
      index <- replicate.index[[y$replicate[1]]]
      lines(y$time, y[, resp],
            type = "l",
            col = colors[index])
      points(y$time, y[, resp],
             pch = pchs[index],
             col = colors[index])
    })

    # title
    title(paste("Conc: ", unique(x$conc), sep = ""))
  })
  
  par(mfrow = c(1, 1))
}

# general full plot (ggplot variant): one subplot for each concentration,
# and one color for each replicate
#' @import ggplot2
dataPlotFullGG <- function(data, xlab, ylab, resp, remove.someLabels) {
  
  time = NULL
  Nsurv = NULL
  
  data$response <- data[,resp]
  
  # create ggplot object Nsurv / time / replicate / conc
  fg <- ggplot(data, aes(time, response, colour = factor(replicate))) +
    geom_point() +
    geom_line() +
    labs(x = xlab, y = ylab) +
    facet_wrap(~conc, ncol = 2) +
    scale_x_continuous(breaks = unique(data$time),
                       labels = if (remove.someLabels) {
                         exclude_labels(unique(data$time))
                       } else {
                         unique(data$time)
                       }
    ) +
    scale_y_continuous(breaks = unique(round(pretty(c(0, max(data[, resp])))))) +
    expand_limits(x = 0, y = 0) +
    theme_minimal()
  
  fd <- fg + theme(legend.position = "none") # remove legend
  
  return(fd)
}

dataPlotFull <- function(data, xlab, ylab, resp, style = "generic",
                         remove.someLabels = FALSE) {

  if (missing(xlab)) xlab <- "Time"

  if (style == "generic")
    dataPlotFullGeneric(data, xlab, ylab, resp)
  else if (style == "ggplot")
    dataPlotFullGG(data, xlab, ylab, resp, remove.someLabels)
  else stop("Unknown plot style")
}

survDataPlotFull <- function(data, xlab, ylab, style = "generic",
                             remove.someLabels = FALSE) {
  dataPlotFull(data, xlab, ylab, "Nsurv", style, remove.someLabels)
}

dataPlotFixedConc <- function(x,
                              xlab,
                              ylab,
                              main,
                              resp,
                              concentration,
                              style = "generic",
                              addlegend = FALSE,
                              remove.someLabels = FALSE) {

  if (missing(xlab)) xlab <- "Time"

  legend.position <- ifelse(resp == "Nsurv", "bottomleft", "topleft")

  # check concentration value
  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select the concentration
  x <- filter(x, x$conc == concentration)

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(x$time, x[, resp],
         type = "n",
         xaxt = "n",
         yaxt = "n",
         main = main,
         xlim = range(x$time),
         ylim = c(0, max(x[, resp])),
         xlab = xlab,
         ylab = ylab)

    # one line by replicate
    by(x, list(x$replicate),
       function(x) {
         lines(x$time, x[,resp], # lines
               col = x$color)
         points(x$time, x[,resp], # points
                pch = 16,
                col = x$color)
       })

    # axis
    axis(side = 1, at = sort(unique(x[, "time"])))
    axis(side = 2, at = unique(round(pretty(c(0, max(x[, resp]))))))

    if (addlegend) {
      legend(legend.position, legend = unique(x$replicate) ,
             col = unique(x$color),
             pch = 16,
             lty = 1)
    }
  }
  else if (style == "ggplot") {
    x$response <- x[,resp]

    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = time, y = response))
    } else {
      df <- ggplot(x, aes(x = time, y = response,
                          color = factor(replicate),
                          group = replicate))
    }
    fd <- df + geom_line() + geom_point() + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_color_hue("Replicate") +
      scale_x_continuous(breaks = unique(x$time),
                         labels = if (remove.someLabels) {
                           exclude_labels(unique(x$time))
                         } else {
                           unique(x$time)
                         }) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(x$response)))))) +
      expand_limits(x = 0, y = 0)

    if (addlegend) {# only if pool.replicate == FALSE
      fd
    } else {
      fd + theme(legend.position = "none") # remove legend
    }
  }
  else stop("Unknown plot style")
}

survDataPlotFixedConc <- function(x,
                                  xlab,
                                  ylab,
                                  main,
                                  concentration,
                                  style = "generic",
                                  addlegend = FALSE,
                                  remove.someLabels = FALSE) {

  dataPlotFixedConc(x, xlab, ylab, main, "Nsurv", concentration,
                    style, addlegend, remove.someLabels)
}

