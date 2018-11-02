#' Plotting method for \code{survFit} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFit}.  It plots the fit obtained for each
#' concentration of chemical compound in the original dataset.
#'
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration.
#' The black dots depict the \strong{observed survival
#' rate} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95\% binomial credible intervals for the estimated survival
#' rate (by default the grey area around the fitted curve) and 95\% binomial confidence
#' intervals for the observed survival rate (as black segments if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two types of  intervals.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted lines limiting the credible band, and a spaghetti plot is added to this band.
#' This spaghetti plot consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (10\% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{survFit}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param spaghetti if \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one plot per concentration.
#' @param adddata if \code{TRUE}, adds the observed data to the plot
#' with (frequentist binomial) confidence intervals.
#' @param mcmc_size A numerical value refering by default to the size of the mcmc in object \code{survFit}.
#'  This option is specific to \code{survFitVarExp} objects for which computing time may be long.
#'  \code{mcmc_size} can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#'@param scales Shape the scale of axis. Default is \code{"fixed"}, but can be \code{"free"}, or free
#' in only one dimension \code{"free_x"}, \code{"free_y"}. (See \code{ggplot2} documentation
#'  for more details.)
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
#' # (4) Summary look the estimated values (parameters)
#' summary(out)
#'
#' # (5) Plot the fitted curve
#' plot(out, adddata = FALSE)
#'
#' # (6) Plot the fitted curve with ggplot style and CI as spaghetti
#' plot(out, spaghetti = TRUE)
#' }
#' 
#' @export
#' 
#' @importFrom stats predict
#' @importFrom tidyr gather
#'
plot.survFitVarExp <- function(x,
                               xlab = "Time",
                               ylab = "Survival rate",
                               main = NULL,
                               spaghetti = FALSE,
                               one.plot = FALSE,
                               adddata = TRUE,
                               mcmc_size = NULL,
                               scales = "fixed",
                               ...) {
  
  
  df_predictTotal <- predict(x, spaghetti = spaghetti, mcmc_size = mcmc_size)
  
  df_prediction <-  df_predictTotal$df_quantile
  df_spaghetti <-  df_predictTotal$df_spaghetti
  
  df_observation <- filter(x$original.data, !is.na(Nsurv))
  
  # Plot
  
  plt <- ggplot() +
      theme_minimal() +
      scale_x_continuous(name = xlab) +
      scale_y_continuous(name = ylab,
                         limits = c(0,1)) +
      theme(legend.position = "top")
  
  # Observation
  if(adddata == TRUE){
    plt <- plt +
      geom_point(data = df_observation,
               aes(x = time, y = Nsurv/Ninit, group = replicate))
  }
  
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
                  fill = "lightgrey", alpha = 0.4)
  }
  
  # Prediction
  plt <- plt +
    geom_line(data = df_prediction,
              aes(x = time, y = q50, group = replicate),
              col="orange", size = 1)
    
  # facetting
  if(one.plot == FALSE){
    plt <- plt + facet_wrap(~ replicate, scales = scales)
  }  
      
   return(plt)
  }

