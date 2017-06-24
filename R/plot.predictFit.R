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
#' 
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#'
#' @export
#'
#'
plot.predictFit <- function(x,
                            xlab = "Time",
                            ylab = "Survival rate",
                            mainlab = NULL,
                            one.plot = FALSE,
                            levels_vector = NULL,
                            adddata = FALSE,
                            addlegend = FALSE){
 
  
  predict_df <- x$predict_data
  
  if(!is.null(levels_vector)){
    predict_df$replicate = factor(predict_df$replicate, levels = levels_vector)
  }
  
  plt_fit <- predict_df %>%
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
    ) +
    geom_ribbon(aes(x = time, ymin = qinf95, ymax = qsup95), fill="pink")+
    geom_line(aes(x = time, y = q50), color= "red")
  
  ##
  ## adddata
  ##
  
  observ_df <- x$observ_data %>%
    group_by(replicate) %>%
    mutate(Ninit = max(Nsurv)) %>%
    ungroup() %>%
    mutate(psurv_obs = Nsurv / Ninit)
  
  if(!is.logical(adddata)) stop ("'adddata' must be a logical.")  
  if(adddata == TRUE){
    plt_fit <- plt_fit +
      geom_point(data = observ_df,
                 aes(x = time,
                     y = psurv_obs,
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
      plt_fit <- plt_fit + facet_wrap(~ replicate)
  }

  return(plt_fit)
}

