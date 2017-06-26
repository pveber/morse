#' Plotting method for \code{survFit} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFit}.  It plots the fits obtained for each
#' concentration of pollutant in the original dataset.
#'
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration (if \code{ data_type = "rate"}), or th
#' \strong{estimated number of survivors} as a function
#' of time for each concentration (if \code{ data_type = "number"})
#' The black dots depict the \strong{observed survival
#' rate} at each time point (if \code{adddata = TRUE}). Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#'
#' @param x An object of class \code{survFit}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param mainlab A mainlab title for the plot.
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one per concentration.
#' @param adddata if \code{TRUE}, adds the observed data to the plot
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param data_type A label \code{'rate'} or \code{'number'}. If \code{'rate'},
#' the fitted curves represent the \strong{estimated survival rate}, if \code{'number'}
#' thef itted curves represent the \strong{estimated number of survivors} as a function
#' of time for each concentration.
#' 
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#'
#' @export
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom dplyr contains
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' @importFrom tibble data_frame
#' @importFrom tibble as_data_frame
#'
plot.survFit <- function(x,
                         xlab = "Time",
                         ylab = NULL,
                         mainlab = NULL,
                         one.plot = FALSE,
                         adddata = TRUE,
                         addlegend = FALSE,
                         data_type = "rate",
                         # facetting = TRUE,
                         #facet.label = "replicate",
                         ...) {
  
  
  ### compute posteriors median and 95 CI
  modelData <- x$modelData
  
  postData <- posteriorData(x$mcmc, model_type = x$model_type)
  
  if(data_type == "rate"){
    
    ylab = "Survival rate"
    
    df_plt = data_frame(Nsurv = modelData$Nsurv,
                        time = modelData$time,
                        replicate = modelData$replicate)
    
    if(!is.null(modelData$conc)){
      df_plt$conc <- modelData$conc
    }
    
    df_plt <- df_plt %>%
      group_by(replicate) %>%
      mutate(Ninit = max(Nsurv)) %>%
      ungroup() %>%
      mutate(Y = Nsurv / Ninit,
             timelag = ifelse(time == 0, time, lag(time)),
             Y_q50 = apply(postData$df_psurv, 2, quantile, probs = 0.5, na.rm = TRUE),
             Y_qinf95 = apply(postData$df_psurv, 2, quantile, probs = 0.025, na.rm = TRUE),
             Y_qsup95 = apply(postData$df_psurv, 2, quantile, probs = 0.975, na.rm = TRUE))
    
    
  } else if(data_type == "number"){
    
    ylab <- "Number of survivors"
    
    df_plt <- data_frame(Y = modelData$Nsurv,
                        time = modelData$time,
                        replicate = modelData$replicate,
                        Y_q50 = apply(postData$df_sim, 2, quantile, probs = 0.5, na.rm = TRUE),
                        Y_qinf95 = apply(postData$df_sim, 2, quantile, probs = 0.025, na.rm = TRUE),
                        Y_qsup95 = apply(postData$df_sim, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
      mutate(timelag = ifelse(time == 0, time, lag(time)))
    
    if(!is.null(modelData$conc)){
      df_plt$conc <- modelData$conc
    }
    
  } else stop("type must be 'rate' (i.e., rate of survival) or 'number' (i.e., number of survivors)")
  
  
  plt_fit <- df_plt %>%
    ggplot() + theme_minimal() +
    expand_limits(x = 0, y = 0) +
    labs(title = mainlab,
         x = xlab,
         y = ylab,
         colour = "Concentration" # legend title
    ) +
    scale_colour_gradient(
      name = "Concentration",
      low = "grey20", high="orange"
    ) +
    scale_fill_gradient(
      name = "Concentration",
      low = "grey20",
      high = "orange"
    )
  
  
  if(data_type == "rate"){
   plt_fit <- plt_fit +
     geom_ribbon(aes(x = time,
                     ymin = Y_qinf95,
                     ymax = Y_qsup95 ,
                     group = replicate), fill = "lightgrey") +
    geom_line(aes(x = time,
                  y = Y_q50,
                  group = replicate ), color="orange")
  } else if(data_type == "number"){
    plt_fit <- plt_fit +
      geom_rect(aes(xmin = timelag, xmax = time,
                    ymin = Y_qinf95,  ymax = Y_qsup95 ,
                    group = replicate), fill = "lightgrey") +
      geom_step(aes(x = time,
                    y = Y_q50,
                    group = replicate ), direction = "vh", color="orange")
  
  } else stop("type must be 'rate' (i.e., rate of survival) or 'number' (i.e., number of survivors)")
  
  ##
  ## adddata
  ##
  if(!is.logical(adddata)) stop ("'adddata' must be a logical.")  
  if(adddata == TRUE){
    plt_fit <- plt_fit +
      geom_point(aes(x = time,
                     y = Y,
                     group = replicate ))
    
  }
  ##
  ## addlegend
  ##
  if(!is.logical(addlegend)) stop ("'addlegend' must be a logical.")
  if(addlegend == TRUE){
    plt_fit <- plt_fit +
      theme(legend.position="top")
  } else {
    plt_fit <- plt_fit +
      theme(legend.position="none")
  }
  
  ##
  ## facetting
  ##
  if(!is.logical(one.plot)) stop ("'one.plot' must be a logical.")  
  if(one.plot != TRUE){
    if(is.null(modelData$conc)){
      plt_fit <- plt_fit + facet_wrap(~ replicate)
    } else{
      plt_fit <- plt_fit + facet_wrap(~ conc)
    }
    # if(facet.label == "replicate"){
    #   plt_fit <- plt_fit + facet_wrap(~ replicate)
    # } else if (facet.label == "conc"){
    #   plt_fit <- plt_fit + facet_wrap(~ conc)
    # } else stop("'facet.label' is either 'replicate' (default) or 'conc'.")
  }
  
  return(plt_fit)
}
