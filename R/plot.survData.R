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

#' @import ggplot2
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
    ylim(0, max(data$Nsurv)) + theme_minimal()

  # legend option
  if (addlegend){
    fd <- fg + scale_colour_hue(title.legend) # the default legend
  } else {
    fd <- fg + theme(legend.position = "none") # remove legend
  }
  return(fd)
}

survDataPlotFull <- function(data, style = "generic", addlegend = TRUE, ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Time"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Number of surviving individuals"

  if (style == "generic")
    survDataPlotFullGeneric(data, xlab, ylab, addlegend)
  else if (style == "ggplot")
    survDataPlotFullGG(data, xlab, ylab, addlegend)
  else stop("Unknown plot style")
}

#' @import ggplot2
#' @importFrom dplyr %>% filter
survDataPlotTargetTime <- function(x, target.time, style, addlegend, ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Concentration"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Number of surviving individuals"

  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")

  # select the target.time
  x <- filter(x, x$time == target.time)

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(x$conc, seq(0, max(x$Nsurv), length.out = length(x$conc)),
         type = "n",
         xaxt = "n",
         xlab = xlab,
         ylab = ylab)

    axis(side = 1, at = unique(x$conc),
         labels = unique(x$conc))

    # points
    if (length(unique(x$replicate)) == 1) {
      # points
      points(x$conc, x$Nsurv,
             pch = 16)
    } else {
      tt <- xyTable(x$conc, x$Nsurv)
      points(tt$x, tt$y,
             cex = (tt$number) / 3,
             pch = 16)
      if (addlegend) {
        legend("bottomleft",
               legend = sort(unique(tt$number)),
               pt.cex = sort(unique((tt$number) / 3)),
               title = "Overplotted replicates",
               pch = 16,
               bty = "n")
      }
    }
  }
  if (style == "ggplot") {
    if (length(unique(x$replicate)) == 1) {
      df <- ggplot(x, aes(x = conc, y = Nsurv))
    } else {
      df <- ggplot(x, aes(x = conc, y = Nsurv)) +
        stat_sum(aes(size = factor(..n..)))
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


survDataPlotFixedConc <- function(x,
                                  concentration,
                                  style = "generic",
                                  addlegend = TRUE,
                                  ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Time"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Number of surviving individuals"

  # check concentration value
  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select the concentration
  x <- filter(x, x$conc == concentration)

  # vector color
  x$color <- as.numeric(as.factor(x$replicate))

  if (style == "generic") {
    plot(x$time, x$Nsurv,
         type = "n",
         xlab = xlab,
         ylab = ylab)

    # one line by replicate
    by(x, list(x$replicate),
       function(x) {
         lines(x$time, x$Nsurv, # lines
               col = x$color)
         points(x$time, x$Nsurv, # points
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
      df <- ggplot(x, aes(x = time, y = Nsurv))
    } else {
      df <- ggplot(x, aes(x = time, y = Nsurv,
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


#' @importFrom dplyr %>% filter
survDataPlotReplicates <- function(x,
                                   target.time,
                                   concentration,
                                   style,
                                   addlegend,
                                   ...) {

  opt_args <- list(...)
  xlab <- if("xlab" %in% names(opt_args)) opt_args[["xlab"]] else "Replicate"
  ylab <- if("ylab" %in% names(opt_args)) opt_args[["ylab"]] else "Number of surviving individuals"

  # check [target.time] and [concentration]
  if (!target.time %in% x$time)
    stop("The argument [target.time] should correspond to one of the observed time points")

  if (!concentration %in% x$conc)
    stop("The argument [concentration] should correspond to one of the tested concentrations")

  # select for concentration and target.time
  x <- filter(x, conc == concentration & time == target.time)

  if (style == "generic") {
    plot(factor(x$replicate), x$Nsurv,
         type = "n",
         xlab = xlab,
         ylab = ylab)
  }

  if (style == "ggplot") {
    df <- ggplot(x, aes(x = replicate, y = Nsurv))
    df + geom_point() + labs(x = xlab, y = ylab) + theme_minimal()
  }
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
#' @param style graphical backend, can be \code{'generic'} or
#' \code{'ggplot'}
#' @param addlegend if \code{TRUE}, a default legend is added to the plot
#' @param pool.replicate If \code{TRUE}, the datapoints of each replicate are
#' summed for a same concentration
#' @param \dots further arguments to be passed to generic methods (xlab, ylab, ...).
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot2}} and returns an object of class \code{ggplot}.
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
#' # (3) Plot the survival data with a ggplot type
#' plot(zinc, style = "ggplot", addlegend = FALSE)
#'
#' # (4) To build a specific legend with a ggplot type
#' fu <- plot(zinc, style = "ggplot", addlegend = FALSE)
#' fu + theme(legend.position = "left") + scale_colour_hue("Replicate")
#'
#' # (5) Plot the survival rate in function of replicates for one concentration at
#' # one target.time with a generic type
#' plot(zinc, style = "generic", target.time = 21, concentration = 0.66)
#'
#' # (6) Plot the survival rate in function of replicates for one concentration at
#' # one target.time with a ggplot type
#' plot(zinc, style = "ggplot", target.time = 21, concentration = 0.66)

#' @export
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
plot.survData <- function(x,
                          target.time = NULL,
                          concentration = NULL,
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
    survDataPlotFull(x, style, addlegend, ...)
  else if (! is.null(target.time) && is.null(concentration))
    survDataPlotTargetTime(x, target.time, style, addlegend, ...)
  else if (is.null(target.time) && ! is.null(concentration))
    survDataPlotFixedConc(x, concentration, style, addlegend, ...)
  else
    survDataPlotReplicates(x, target.time, concentration, style, addlegend, ...)
}
