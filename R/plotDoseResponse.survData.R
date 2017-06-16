#' Plotting method for \code{survData} objects
#'
#' This is the generic \code{plotDoseResponse} S3 method for the \code{survData}
#' class. It plots the survival rate as a function of concentration (for a given
#' target time).
#' 
#' The function plots the observed values of the survival rate for a given time
#' as a function of concentration. The 95 \% binomial confidence interval is added
#' to each survival rate. It is calculated using function
#' \code{\link[stats]{binom.test}} from package \code{stats}.
#' Replicates are systematically pooled in this plot.
#'
#' @param x an object of class \code{survData}
#' @param xlab a title for the \eqn{x}-axis (optional)
#' @param ylab a label for the \eqn{y}-axis
#' @param main main title for the plot
#' @param target.time a numeric value corresponding to some observed time in \code{data}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param log.scale if \code{TRUE}, displays \eqn{x}-axis in log scale
#' @param remove.someLabels if \code{TRUE}, removes 3/4 of X-axis labels in
#' \code{'ggplot'} style to avoid the label overlap
#' @param addlegend if \code{TRUE}, adds a default legend to the plot
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @note When \code{style = "ggplot"}, the function calls function
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
#' 
#' @seealso \code{\link[stats]{binom.test}}
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
#' # (2) Plot dose-response
#' plotDoseResponse(zinc)
#'
#' # (3) Plot dose-response with a ggplot style
#' plotDoseResponse(zinc, style = "ggplot")
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom grid arrow unit
#' @importFrom graphics plot axis lines points segments title
#' @importFrom methods is
#' @importFrom stats aggregate
#'
#' @export
plotDoseResponse.survDataCstC <- function(x,
                                      xlab = "Concentration",
                                      ylab = "Survival rate",
                                      main = NULL,
                                      target.time = NULL,
                                      style = "generic",
                                      log.scale = FALSE,
                                      remove.someLabels = FALSE,
                                      addlegend = TRUE,
                                      ...) {
  if (is.null(target.time)) target.time <- max(x$time)
  
  if (!target.time %in% x$time || target.time == 0)
    stop("[target.time] is not one of the possible time !")
  
  if (style == "generic" && remove.someLabels)
    warning("'remove.someLabels' argument is valid only in 'ggplot' style.",
            call. = FALSE)
  
  # Create a new column named profile
  x$profile = as.character(x$conc)
  
  # agregate by sum of profile
  x <- x %>%
    dplyr::group_by(profile,conc, time) %>%
    dplyr::summarise(Nsurv = sum(Nsurv)) %>%
    ungroup()
  
  # create column Ninit
  x <- x %>%
    dplyr::group_by(profile) %>%
    # survDataCheck checked that Nsurv was decreasing in each profile and the present of time == 0
    dplyr::mutate(Ninit = max(Nsurv))
  
  
  x$resp <- x$Nsurv / x$Ninit
  # select the target.time
  xf <- filter(x, x$time == target.time)
  
  conf.int <- survConfInt(xf, log.scale)
  
  # Selection of datapoints that can be displayed given the type of scale
  sel <- if(log.scale) xf$conc > 0 else TRUE
  x <- xf[sel, ]
  transf_data_conc <- optLogTransform(log.scale, x$conc)
  
  # Concentration values used for display in linear scale
  display.conc <- (function() {
    x <- optLogTransform(log.scale, x$conc)
    if(log.scale) exp(x) else x
  })()
  
  # vector color
  # x$color <- as.numeric(as.factor(x$profile))
  
  if (style == "generic") {
    plot(transf_data_conc, seq(0, max(conf.int["qsup95",]),
                               length.out = length(transf_data_conc)),
         type = "n",
         xaxt = "n",
         yaxt = "n",
         main = main,
         xlab = xlab,
         ylab = ylab)
    
    axis(side = 1, at = transf_data_conc,
         labels = display.conc)
    axis(side = 2, at = unique(round(pretty(c(0, max(x$resp))))),
         labels = unique(round(pretty(c(0, max(x$resp))))))
    
    # points
    points(transf_data_conc, x$resp,
           pch = 20)
    
    # segment CI
    
    segments(transf_data_conc, x$resp,
             transf_data_conc, conf.int["qsup95", ])
    
    segments(transf_data_conc, x$resp,
             transf_data_conc, conf.int["qinf95", ])
    
    # add legend
    if (addlegend) {
      legend("bottomleft", pch = c(20, NA),
             lty = c(NA, 1),
             lwd = c(NA, 1),
             col = c(1, 1),
             legend = c("Observed values", "Confidence intervals"),
             bty = "n")
    }
  }
  else if (style == "ggplot") {
    # colors
    valCols <- fCols(x, fitcol = NA, cicol = NA)
    
    df <- data.frame(x,
                     transf_data_conc,
                     display.conc,
                     Points = "Observed values")
    dfCI <- data.frame(conc = transf_data_conc,
                       qinf95 = conf.int["qinf95",],
                       qsup95 = conf.int["qsup95",],
                       Conf.Int = "Confidence intervals")
    
    fd <- ggplot(df) +
      geom_point(aes(x = transf_data_conc, y = resp, fill = Points),
                 data = df, col = valCols$cols1) +
      geom_segment(aes(x = conc, xend = conc, y = qinf95,
                       yend = qsup95,
                       linetype = Conf.Int),
                   dfCI, col = valCols$cols3) +
      scale_fill_hue("") +
      scale_linetype(name = "") +
      expand_limits(x = 0, y = 0) + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = unique(df$transf_data_conc),
                         labels = if (remove.someLabels) {
                           exclude_labels(unique(df$display.conc))
                         } else {
                           unique(df$display.conc)
                         }
      ) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(df$resp)))))) +
      expand_limits(x = 0, y = 0)
    
if (addlegend) {
  fd
} else {
  fd + theme(legend.position = "none") # remove legend
}
  }
  else stop("Unknown plot style")
}

#' @importFrom stats aggregate binom.test
survConfInt <- function(x, log.scale) {
  # compute confidente interval on observed data
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # - log.scale : boolean
  # OUTPUT:
  
  # - ci : confidente interval
  ci <- apply(x, 1, function(x) {
    binom.test(x["Nsurv"], x["Ninit"])$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$conc
  
  if (log.scale) ci <- ci[ ,colnames(ci) != 0]
  
  return(ci)
}
