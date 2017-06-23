#' Plotting method for \code{predictFit} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFit}.  It plots the fits obtained for each
#' concentration of pollutant in the original dataset.
#'
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration.
#' The black dots depict the \strong{observed survival
#' rate} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95 \% credible intervals for the estimated survival
#' rate (by default the red area around the fitted curve) and 95 \% confidence
#' intervals for the observed survival rate (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (2 \% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{gm_urvFitTKTD}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param mainlab A mainlab title for the plot.
#' @param concentration A numeric value corresponding to some concentration in
#' \code{data}. If \code{concentration = NULL}, draws a plot for each concentration.
#' @param spaghetti if \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one per concentration.
#' @param adddata if \code{TRUE}, adds the observed data to the plot
#' with (frequentist) confidence intervals
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#'
#' @export
#'
#' @importFrom deSolve 
#'
plot.predictFit <- function(x){
 
  data_plt <- x$predict
  
    plt <- ggplot(data = data_plt) +
      theme_minimal() +
      # theme(legend.position ="top",
      #       strip.background = element_rect(fill="grey10"),
      #       strip.text = element_text(colour = "orange", face= "bold")) +
      scale_x_continuous(name = "time")+
      scale_y_continuous(name = "survival rate",
                         limits = c(0,1)) +
      theme(legend.position = "top") +
      geom_ribbon(aes(x = time, ymin = qinf95, ymax = qsup95), fill = "pink", alpha = 0.4) +
      geom_line(aes(x = time, y = q50), col = "red", size = 1)
    
    return(plt)
  }
