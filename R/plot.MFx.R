#' Plotting method for \code{MFx} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \\code{MFx} class. It plots the survival rate as a function of concentration.
#'
#'
#' @param x An object of class \code{MFx}.
#' @param x_variable A character to define the variable for the \eqn{X}-axis,
#'  either \code{"MFx"} or \code{"Time"}. The default is \code{"MFx"}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{NULL} and depend on the
#' argument \code{x_variable}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate median and 95 CI}.
#' @param main A main title for the plot.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#' 
#' @examples 
#' 
#' # (1) Load the data
#' data("propiconazole")
#' 
#' # (2) Create an object of class 'survData'
#' dataset <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFit function with model_type SD (or IT)
#' out_SD <- survFit(dataset, model_type = "SD")
#' 
#' # (4) estimate MF50 at time 4
#' MFx_SD <- MFx(out_SD, x_MFx = 50, time_MFx = 4)
#' 
#' # (5) plot the object of class 'MFx'
#' plot(MFx_SD)
#' }
#'
#' @export
#'
#'
#'
plot.MFx <- function(x,
                     spaghetti = FALSE,
                     x_variable = "MFx", # other option is "Time"
                     xlab = NULL,
                     ylab = "Survival rate \n median and 95 CI",
                     main = NULL,
                     log_scale = FALSE,
                     ...){
  
  # definition of xlab and check x_variable
  if(is.null(xlab)){
    if(x_variable == "MFx"){
      xlab = "Multiplication Factor"
    } else if(x_variable == "Time"){
      xlab = "Time"
    } else stop("the argument 'x_variable' must be 'MFx' or 'Time'. The default is 'MFx'.")
  }
  
  if(is.null(main)){
    main <- paste("Multiplication Factor MF is", 
                  round(x$optimal_MFx$MFx, digits = 3))
  } 
  
  # Load data frame from x
  if(x_variable == "MFx"){
    
    df_MFx <- x$df_MFx
    
    legend.point = data.frame(
      x.pts = rep(x$optimal_MFx$MFx, 3),
      y.pts = c(x$optimal_MFx$q50, x$optimal_MFx$qinf95, x$optimal_MFx$qsup95),
      pts.leg = c(paste("median: ", round(x$optimal_MFx$q50,digits = 2)),
                  paste("quantile 2.5%: ", round(x$optimal_MFx$qinf95, digits = 2)),
                  paste("quantile 97.5%: ", round(x$optimal_MFx$qsup95, digits = 2)))
    )
    
    
    MFx_plt <- ggplot() + theme_bw() +
      labs(title = main) +
      scale_x_continuous(name = xlab) +
      scale_y_continuous(name = ylab, limits = c(0,1)) +
      # Initial prediction:
      geom_ribbon(data = df_MFx,
                  aes(x = MFx, ymin = qinf95, ymax = qsup95),
                  fill = "lightgrey") +
      geom_line(data = df_MFx,
                aes(x = MFx, y = q50), color = "orange") +
      geom_vline(xintercept = x$optimal_MFx$MFx, col="grey70", linetype=2) +
      geom_point(data = legend.point,
                 aes(x = x.pts, y=y.pts,
                     color= pts.leg)) +
      scale_color_manual(values=c("orange", "black", "black"))
    
    if(log_scale == TRUE){
      MFx_plt <- MFx_plt +
        scale_x_log10()
    }
    
  }
  
  if(x_variable == "Time"){
    
    initial_prediction <- dplyr::filter(x$df_ls_predict, MFx == 1)
    final_prediction <- dplyr::filter(x$df_ls_predict, MFx == x$optimal_MFx$MFx)
    
    y_arrow <- dplyr::filter(initial_prediction, time == x$time_MFx)$q50
    yend_arrow <- dplyr::filter(final_prediction, time == x$time_MFx)$q50
    
    # Plot
    MFx_plt <- ggplot() + theme_bw() +
      labs(title = main) +
      scale_x_continuous(name = xlab) +
      scale_y_continuous(name = ylab, limits = c(0,1)) +
      # Initial prediction:
      geom_ribbon(data = initial_prediction,
                  aes(x = time, ymin = qinf95, ymax = qsup95),
                  fill = "grey80", alpha = 0.4) +
      geom_line(data = initial_prediction,
                aes(x = time, y = q50),
                col="orange", size = 0.3) +
      # Final prediction:
      geom_ribbon(data = final_prediction,
                  aes(x = time, ymin = qinf95, ymax = qsup95),
                  fill = "grey30", alpha = 0.4) +
      geom_line(data = final_prediction,
                aes(x = time, y = q50),
                col="orange", size = 1) +
      geom_segment(aes(x = x$time_MFx, y = y_arrow,
                       xend = x$time_MFx, yend =  yend_arrow),
                   arrow = arrow(length = unit(0.2,"cm")))
    
  }
  
  
  return(MFx_plt)
  
}

