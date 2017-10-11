#' Plotting method for \code{survFitPredict} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFitPredict}.  It plots the fit obtained for each
#' concentration of chemical compound in the provided dataset.
#'
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration.
#' The black dots depict the \strong{observed survival
#' rate} at each time point.
#' The function plots both 95\% binomial credible intervals for the estimated survival
#' rate.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted lines limiting the credible band, and a spaghetti plot is added to this band.
#' This spaghetti plot consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (10\% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{survFitPredict}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param spaghetti if \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one plot per concentration.
#' @param mcmc_size A numerical value refering by default to the size of the mcmc in object \code{survFitPredict}.
#'  This option is specific to \code{survFitPredict} objects for which computing time may be long.
#'  \code{mcmc_size} can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#'  
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#' 
#' @examples 
#'
#' # (1) Load the survival data
#' data("propiconazole_pulse_exposure")
#'
#' # (2) Create an object of class "survData"
#' dataset <- survData(propiconazole_pulse_exposure)
#'
#' \dontrun{
#' # (3) Run the survFit function
#' out <- survFit(dataset , model_type = "SD")
#'
#' # (4) Create a new data table for prediction
#' data_4prediction <- data.frame(time = 1:10, conc = c(0,5,5,5,0,0,5,5,5,5), replicate= rep("predict", 10))
#'
#' # (5) Predict on a new data set
#' predict_out <- predict(x = out, data_predict = data_4prediction, spaghetti = TRUE)
#'
#' # (6) Plot the predicted curve
#' plot(predict_out)
#' plot(predict_out, spaghetti = TRUE)
#' }
#' 
#' @export
#' 
#' @importFrom tidyr gather
#'
plot.survFitPredict <- function(x,
                               xlab = "Time",
                               ylab = "Survival rate",
                               main = NULL,
                               spaghetti = FALSE,
                               one.plot = FALSE,
                               mcmc_size = NULL,
                               ...) {

  df_prediction <-  x$df_quantile
  df_spaghetti <-  x$df_spaghetti
  
  # Plot
  
  plt <- ggplot() +
    theme_minimal() +
    theme(legend.position ="top",
          strip.background = element_rect(fill="grey90", color = "white"),
          strip.text = element_text(colour = "black")) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab,
                       limits = c(0,1)) +
    theme(legend.position = "top")
  
  # spaghetti
  if(spaghetti == TRUE){
    
    df_spaghetti_gather <- df_spaghetti %>%
      tidyr::gather(survRate_key, survRate_value, -c(time,conc,replicate))
    
    plt <- plt +
      geom_line(data = df_spaghetti_gather,
                aes(x = time, y = survRate_value, group = interaction(survRate_key, replicate)),
                alpha = 0.02) +
      geom_line(data = df_prediction,
                aes(x = time, y= qinf95, group = replicate),
                color = "orange", linetype = 2) +
      geom_line(data = df_prediction,
                aes(x = time, y = qsup95, group = replicate),
                color = "orange", linetype = 2)
  }
  if(spaghetti != TRUE){
    plt <- plt + 
      geom_ribbon(data = df_prediction,
                  aes(x = time, ymin = qinf95,ymax = qsup95, group = replicate),
                  fill = "grey30", alpha = 0.4)
  }
  
  # Prediction
  plt <- plt +
    geom_line(data = df_prediction,
              aes(x = time, y = q50, group = replicate),
              col="orange", size = 1)
  
  # facetting
  if(one.plot == FALSE){
    plt <- plt + facet_wrap(~ replicate)
  }  
  
  return(plt)
}

