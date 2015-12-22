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


#' Plotting method for \code{survData} objects
#'
#' Plots the survival rate as a function of concentration (for a given target
#' time).
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
#' @param \dots Further arguments to be passed to generic methods.
#' @note When \code{style = "ggplot"}, the function calls package
#' \code{\link[ggplot2]{ggplot}} and returns an object of class \code{ggplot}.
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
plotDoseResponse.survData <- function(x,
                                      xlab = "Concentration",
                                      ylab = "Survival rate",
                                      main = NULL,
                                      target.time = NULL,
                                      style = "generic",
                                      log.scale = FALSE,
                                      remove.someLabels = FALSE,
                                      ...) {
  if (is.null(target.time)) target.time <- max(x$time)
  
  if (!target.time %in% x$time)
    stop("[target.time] is not one of the possible time !")
  
  # agregate by sum of replicate
  x <- cbind(aggregate(cbind(Nsurv, Ninit) ~ time + conc, x, sum),
             replicate = 1)
  
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
  x$color <- as.numeric(as.factor(x$replicate))
  
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
           pch = 16)
    
    # segment CI
    
    segments(transf_data_conc, x$resp,
             transf_data_conc, conf.int["qsup95", ])
    
    Bond <- if (log.scale) {
      0.03 * (max(transf_data_conc) - min(transf_data_conc))
    } else {
      0.03 * (max(transf_data_conc) - min(transf_data_conc[which(transf_data_conc != 0)]))
    }
    
    segments(transf_data_conc - Bond,
             conf.int["qsup95", ],
             transf_data_conc + Bond,
             conf.int["qsup95", ])
    
    segments(transf_data_conc, x$resp,
             transf_data_conc, conf.int["qinf95", ])
    
    segments(transf_data_conc - Bond,
             conf.int["qinf95", ],
             transf_data_conc + Bond,
             conf.int["qinf95", ])
  }
  else if (style == "ggplot") {
    df <- data.frame(x,
                     transf_data_conc,
                     display.conc)
    dfCI <- data.frame(conc = transf_data_conc,
                       qinf95 = conf.int["qinf95",],
                       qsup95 = conf.int["qsup95",],
                       Conf.Int = "Confidence interval")
    
    gp <- ggplot(df, aes(x = transf_data_conc, y = resp)) +
      geom_segment(aes(x = conc, xend = conc, y = qinf95,
                       yend = qsup95,
                       linetype = Conf.Int),
                   arrow = arrow(length = unit(0.25 , "cm"), angle = 90,
                                 ends = "both"), dfCI) +
      expand_limits(x = 0, y = 0)
    
    fd <- gp + geom_point() + ggtitle(main) +
      theme_minimal() +
      labs(x = xlab,
           y = ylab) +
      scale_x_continuous(breaks = df$transf_data_conc,
                         labels = if (remove.someLabels) {
                           exclude_labels(df$display.conc)
                         } else {
                           df$display.conc
                         }
      ) +
      scale_y_continuous(breaks = unique(round(pretty(c(0, max(df$resp)))))) +
      expand_limits(x = 0, y = 0)
    
    fd + theme(legend.position = "none") # remove legend
    
  }
  else stop("Unknown plot style")
}

