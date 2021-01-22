#' Plotting method for \code{MFx} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{MFx} class. It plots the survival probability as a function of
#' the multiplication factor applied or as a function of time.
#'
#'
#' @param x An object of class \code{MFx}.
#' @param x_variable A character to define the variable for the \eqn{X}-axis,
#'  either \code{"MFx"} or \code{"Time"}. The default is \code{"MFx"}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{NULL} and depend on the
#' argument \code{x_variable}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival probability median and 95 CI}.
#' @param main A main title for the plot.
#' @param log_scale If \code{TRUE}, the x-axis is log-scaled. Default is \code{FALSE}.
#' @param ncol An interger for the number of columns when several panels are plotted.
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
#'# (4) data to predict
#' data_4prediction <- data.frame(time = 1:10, conc = c(0,0.5,3,3,0,0,0.5,3,1.5,0))
#' 
#' # (5) estimate MF for 30% reduction of survival at time 4
#' MFx_SD_30.4 <- MFx(out_SD, data_predict = data_4prediction , X = 30, time_MFx = 4)
#' 
#' # (6) plot the object of class 'MFx'
#' plot(MFx_SD_30.4)
#' 
#' # (6bis) plot with log-scale of x-axis
#' plot(MFx_SD_30.4, log_scale = TRUE)
#' 
#' # (6ter) plot with "Time" as the x-axis
#' plot(MFx_SD_30.4, x_variable = "Time") 
#' 
#' # (7) plot when X = NULL and along a MFx_range from 5 to 10:
#' MFx_SD_range <- MFx(out_SD, data_predict = data_4prediction ,
#'                     X = NULL, time_MFx = 4, MFx_range = seq(5, 10, length.out = 50))
#' plot(MFx_SD_range)
#' plot(MFx_SD_range, x_variable = "Time", ncol = 10)
#' }
#'
#' @export
#'
#'
#'
plot.MFx <- function(x,
                     x_variable = "MFx", # other option is "Time"
                     xlab = NULL,
                     ylab = "Survival probability \n median and 95 CI",
                     main = NULL,
                     log_scale = FALSE,
                     ncol = 3,
                     ...){
  
  # definition of xlab and check x_variable
  if(is.null(xlab)){
    if(x_variable == "MFx"){
      xlab = "Multiplication Factor"
    } else if(x_variable == "Time"){
      xlab = "Time"
    } else stop("the argument 'x_variable' must be 'MFx' or 'Time'. The default is 'MFx'.")
  }
  
  
  MFx_plt <- ggplot() + theme_minimal() +
    theme(legend.position = "top",
          legend.title = element_blank())+
    scale_y_continuous(limits = c(0,1)) 
  
  if(x_variable == "MFx"){
    
    if(is.null(main)){
      main <- paste("Multiplication Factor for MF",  x$X_prop_provided*100, "% at time", x$time_MFx)
    } 
    if(is.null(x$X_prop))  main <- paste("Survival over [", min(x$MFx_tested), ",", max(x$MFx_tested), "] MF range at time", x$time_MFx) 
    
    MFx_plt <- MFx_plt +
      scale_color_manual(values=c("orange", "black", "black")) +
      labs(title = main,
           x = xlab,
           y = ylab) +
      geom_ribbon(data = x$df_dose,
                  aes(x = MFx, ymin = qinf95, ymax = qsup95),
                     fill = "lightgrey") +
      geom_point(data = x$df_dose,
                 aes(x = MFx, y = q50), color = "orange", shape = 4) +
      geom_line(data = x$df_dose,
                aes(x = MFx, y = q50), color = "orange")
    
    if(!is.null(x$X_prop)){
      
      MFx_plt <- MFx_plt +
        geom_point(data = dplyr::filter(x$df_dose, id == "q50"),
                   aes(x = MFx, y = q50), color = "orange", shape = 4)  +
        geom_point(data = dplyr::filter(x$df_dose, id == "qinf95"),
                   aes(x = MFx, y = qinf95), color = "grey", shape = 4)  +
        geom_point(data = dplyr::filter(x$df_dose, id == "qsup95"),
                   aes(x = MFx, y = qsup95), color = "grey", shape = 4)
      
      legend.point = data.frame(
        x.pts = x$df_MFx$MFx,
        y.pts = rep(x$X_prop, 3),
        pts.leg = c(paste("median: ", round(x$df_MFx$MFx[1],digits = 2)),
                    paste("quantile 2.5%: ", round(x$df_MFx$MFx[2],digits = 2)),
                    paste("quantile 97.5%: ", round(x$df_MFx$MFx[3],digits = 2)))
      )
      
      
      MFx_plt <- MFx_plt +
        geom_hline(yintercept = x$X_prop, col="grey70", linetype=2) +
        geom_point(data = legend.point,
                   aes(x = x.pts, y = y.pts, color = pts.leg)) 
      
      warning("This is not an error message:
Just take into account that MFx as been estimated with a binary
search using the 'accuracy' argument. Cross point indicate the
position of evaluated time series. To improve the shape of the curve, you 
can use X = NULL, and compute time series around the median MFx, with the
          vector `MFx_range`.")
      }
    if(log_scale == TRUE){
      MFx_plt <- MFx_plt +
        scale_x_log10()
    }
    
  }
  
  if(x_variable == "Time"){
    
    # Plot
    if(is.null(main))  main <- paste("Survival over time. Multiplication Factor of", x$X_prop_provided*100, "percent") 
    if(is.null(x$X_prop))  main <- paste("Survival over time") 
    
    MFx = x$MFx_tested
    
    k <- 1:length(x$ls_predict)
    #
    # Make a dataframe with quantile of all generated time series
    #
    ls_predict_quantile <- lapply(k, function(kit){
      df_quantile <- x$ls_predict[[kit]]$df_quantile
      df_quantile$MFx <- rep(MFx[[kit]], nrow(x$ls_predict[[kit]]$df_quantile))
      return(df_quantile)
    })
    predict_MFx_quantile <- do.call("rbind", ls_predict_quantile)
    
    
    if(!is.null(x$X_prop)){
      
      initial_predict <- dplyr::filter(predict_MFx_quantile, MFx == 1)
      final_predict <- dplyr::filter(predict_MFx_quantile, MFx == x$df_MFx$MFx[1])
      
      y_arrow <- as.numeric(dplyr::filter(initial_predict, time == x$time_MFx)$q50)
      yend_arrow <- as.numeric(dplyr::filter(final_predict, time == x$time_MFx)$q50)
      
      MFx_plt <- MFx_plt +
        labs(title = main,
             x = xlab,
             y = ylab) +
        # final predict
        geom_ribbon(data = final_predict,
                    aes(x = time, ymin = qinf95, ymax = qsup95),
                    fill = "grey30", alpha = 0.4) +
        geom_line(data = final_predict,
                  aes(x = time, y = q50),
                  col="orange", size = 1) +
        # initial predict
        geom_ribbon(data = initial_predict,
                    aes(x = time, ymin = qinf95, ymax = qsup95),
                    fill = "grey60", alpha = 0.4) +
        geom_line(data = initial_predict,
                  aes(x = time, y = q50),
                  col="orange", size = 0.3) +
        # arrow
        geom_segment(aes(x = x$time_MFx, y = y_arrow,
                         xend = x$time_MFx, yend =  yend_arrow),
                     arrow = arrow(length = unit(0.2,"cm")))
      
      }
    if(is.null(x$X_prop)){
      
      MFx_plt <- MFx_plt +
        labs(title = main,
             x = xlab,
             y = ylab) +
        geom_ribbon(data = predict_MFx_quantile,
                    aes(x = time, ymin = qinf95, ymax = qsup95),
                    fill = "grey30", alpha = 0.4) +
        geom_line(data = predict_MFx_quantile,
                  aes(x = time, y = q50),
                  col="orange", size = 1) +
        facet_wrap( ~ round(MFx, digits = 2), ncol = ncol)
    }
  }
  
  
  return(MFx_plt)
  
}
