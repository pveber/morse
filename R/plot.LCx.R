#' Plotting method for \code{LCx} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{plot} of Lethal Concentration for \code{x} percent of the population.
#'
#'
#' @param x An object of class \code{LCx}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
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
#' dat <- survData("propiconazole")
#' 
#' ## Not run
#' # (3) Run the survFit function with model_type SD (or IT)
#' out_SD <- survFit(dat, model_type = "SD")
#' 
#' # (4) estimate LC50 at time 4
#' LCx_SD <- LCx(out_SD, X = 50, time_LCx = 4)
#' 
#' # (5) plot the object of class 'LCx'
#' plot(LCx_SD)
#'
#' @export
#'
#'
#'
plot.LCx <- function(x,
                     xlab = "Concentration",
                     ylab = "Survival rate \n median and 95CI",
                     main = NULL, ...){
  
  df_dose <- x$df_dose
  df_LCx <- x$df_LCx
  X_prop <- x$X_prop
  time_LCx <- x$time_LCx
  
  # loc.x.legend = max(df_dose$concentration)-0.4*max(df_dose$concentration)
  
  legend.point = data.frame(
    x.pts = df_LCx$LCx,
    y.pts = rep(X_prop,3),
    pts.leg = c(paste("median: ", round(df_LCx$LCx[1],digits = 2)),
                paste("quantile 2.5%: ", round(df_LCx$LCx[2],digits = 2)),
                paste("quantile 97.5%: ", round(df_LCx$LCx[3],digits = 2)))
  )
  
  if(is.null(main)){
    main <- paste("Dose response curve: LC for", X_prop*100, "% at time", time_LCx)
  } 
  
  
  LCx_plt <- ggplot() + theme_minimal() +
    scale_y_continuous(limits = c(0,1)) +
    labs(title = main,
         x = xlab,
         y = ylab) +
    geom_ribbon(data = df_dose,
                aes(x = concentration, ymin = qinf95, ymax = qsup95), fill = "lightgrey") + 
    geom_line(data = df_dose,
              aes( x = concentration, y = q50), color ="orange") +
    # LCx points
    geom_hline(yintercept = X_prop, col="grey70", linetype=2) +
    
    geom_point(data = legend.point,
               aes(x = x.pts, y=y.pts,
                   color= pts.leg))+
    theme(legend.position = "top",
          legend.title = element_blank())+
    scale_color_manual(values=c("orange", "black", "black"))
    
  LCx_plt 
  
  
  return(LCx_plt)
  
}
