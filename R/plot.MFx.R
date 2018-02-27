#' Plotting method for \code{MFx} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \\code{MFx} class. It plots the survival rate as a function of concentration.
#'
#'
#' @param x An object of class \code{MFx}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
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
                     xlab = "Time",
                     ylab = "Survival rate \n median and 95 CI",
                     main = NULL, ...){
  
  Multiplication_factor = x$Multiplication_factor
  time_MFx = x$time_MFx
  initial_x_MFx = x$initial_x_MFx
  theoretical_x_MFx = x$theoretical_x_MFx
  final_x_MFx = x$final_x_MFx
  initial_prediction = x$initial_prediction
  final_prediction = x$final_prediction
  
  if(is.null(main)){
    main <- paste("Multiplication Factor-Response curve: MF is", 
                  round(Multiplication_factor, digits = 3),
                  "(surv. rate=", round(final_x_MFx, digits = 2), 
                  "; t=", round(time_MFx, digits = 2), ")")
  } 
  
  # Plot
  
  MFx_plt <- ggplot() + theme_bw() +
    labs(title = main) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab, limits = c(0,1)) +
    # Initial prediction:
    geom_ribbon(data = initial_prediction$df_quantile,
                aes(x = time, ymin = qinf95, ymax = qsup95),
                fill = "grey80", alpha = 0.4) +
    geom_line(data = initial_prediction$df_quantile,
              aes(x = time, y = q50),
              col="orange", size = 0.3) +
    # Final prediction:
    geom_ribbon(data = final_prediction$df_quantile,
                aes(x = time, ymin = qinf95, ymax = qsup95),
                fill = "grey30", alpha = 0.4) +
    geom_line(data = final_prediction$df_quantile,
              aes(x = time, y = q50),
              col="orange", size = 1) +
    geom_segment(aes(x = time_MFx, y = initial_x_MFx,
                     xend = time_MFx, yend = final_x_MFx),
                 arrow = arrow(length = unit(0.2,"cm")))
  
  MFx_plt
  return(MFx_plt)
  
}

