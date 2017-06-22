#' Plotting method for \code{survFit} objects
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
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#'
plot.survFit <- function(x,
                         xlab = "Time",
                         ylab = NULL,
                         mainlab = NULL,
                         data_type = "probability",
                         facetting = TRUE,
                         facet.label = "replicate") {
  
  
  ### compute posteriors median and 95 CI
  modelData = x$modelData
  
  postData = posteriorData(x$mcmc, model_type = x$model_type)
  
  if(data_type == "probability"){
    
    ylab = "Survival rate"
    
    df_plt = data_frame(Nsurv = modelData$Nsurv,
                        time = modelData$time,
                        replicate = modelData$replicate)
    
    if(!is.null(modelData$conc)){
      df_plt$conc = modelData$conc
    }
    
    df_plt = df_plt %>%
      group_by(replicate) %>%
      mutate(Ninit = max(Nsurv)) %>%
      ungroup() %>%
      mutate(Y = Nsurv / Ninit,
             timelag = ifelse(time == 0, time, lag(time)),
             Y_q50 = apply(postData$df_psurv, 2, quantile, probs = 0.5, na.rm = TRUE),
             Y_qinf95 = apply(postData$df_psurv, 2, quantile, probs = 0.025, na.rm = TRUE),
             Y_qsup95 = apply(postData$df_psurv, 2, quantile, probs = 0.975, na.rm = TRUE))
    
    
  } else if(data_type == "number"){
    
    ylab = "Number of survivors"
    
    df_plt = data_frame(Y = modelData$Nsurv,
                        time = modelData$time,
                        replicate = modelData$replicate,
                        Y_q50 = apply(postData$df_sim, 2, quantile, probs = 0.5, na.rm = TRUE),
                        Y_qinf95 = apply(postData$df_sim, 2, quantile, probs = 0.025, na.rm = TRUE),
                        Y_qsup95 = apply(postData$df_sim, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
      mutate(timelag = ifelse(time == 0, time, lag(time)))
    
    if(!is.null(modelData$conc)){
      df_plt$conc = modelData$conc
    }
    
  } else stop("type must be 'probability' (i.e., probability of survival) or 'number' (i.e., number of survivors)")
  
  
  plt_fit <- df_plt %>%
    ggplot() + theme_minimal() +
    theme(legend.position="none") +
    expand_limits(x = 0, y = 0) +
    labs(title = mainlab,
         x = xlab,
         y = ylab,
         colour = "Concentration" # legend title
    ) +
    scale_colour_gradient(
      name="Concentration",
      low="grey20", high="orange"
    ) +
    scale_fill_gradient(
      name="Concentration",
      low="grey20", high="orange"
    ) +
    # geom_ribbon(aes(x = time,
    #                 ymin = Y_qinf95,
    #                 ymax = Y_qsup95 ,
    #                 group = replicate), fill = "orange", alpha = 0.4, color = "grey90") +
    geom_rect(aes(xmin = timelag, xmax = time,
                  ymin = Y_qinf95,  ymax = Y_qsup95 ,
                  group = replicate), fill = "pink", alpha = 0.4) +
    geom_step(aes(x = time,
                  y = Y_q50,
                  group = replicate ), direction = "vh", color="red") +
    geom_point(aes(x = time,
                   y = Y,
                   group = replicate ))
  
  ##
  ## facetting
  ##
  if(facetting == TRUE){
    if(facet.label == "replicate"){
      OUT_fit <- plt_fit + facet_wrap(~ replicate)
    } else if (facet.label == "conc"){
      OUT_fit <- plt_fit + facet_wrap(~ conc)
    } else stop("'facet.label' is either 'replicate' (default) or 'conc'.")
  }
  
  return(OUT_fit)
}
