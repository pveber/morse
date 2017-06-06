#' Plotting method for \code{gm_survFitTKTD} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{gm_survFitTKTD}.  It plots the fits obtained for each
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
#' @param main A main title for the plot.
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
plot.gm_survFitTKTD <- function(x,
                             lab.x = "Time",
                             lab.y = NULL,
                             lab.main = NULL,
                             data_type = "probability") {

  ### check class
  if (!is(x, "gm_survFitTKTD")){
    stop("x is not of class 'gm_survFitTKTD'!")
  }

  ### compute posteriors median and 95 CI
  modelData = x$modelData

  postData = posteriorData(x$mcmc, model_type = x$model_type)

  if(data_type == "probability"){

    lab.y = "Survival rate"

    df_plt = data_frame(Nsurv = modelData$Nsurv,
                        time = modelData$time,
                        profile = modelData$profile) %>%
      group_by(profile) %>%
      mutate(Ninit = max(Nsurv)) %>%
      ungroup() %>%
      mutate(Y = Nsurv / Ninit,
             Y_q50 = apply(postData$df_psurv, 2, quantile, probs = 0.5, na.rm = TRUE),
             Y_qinf95 = apply(postData$df_psurv, 2, quantile, probs = 0.025, na.rm = TRUE),
             Y_qsup95 = apply(postData$df_psurv, 2, quantile, probs = 0.975, na.rm = TRUE))


  } else if(data_type == "number"){

    lab.y = "Number of survivors"

    df_plt = data_frame(Y = modelData$Nsurv,
                        time = modelData$time,
                        profile = modelData$profile,
                        Y_q50 = apply(postData$df_sim, 2, quantile, probs = 0.5, na.rm = TRUE),
                        Y_qinf95 = apply(postData$df_sim, 2, quantile, probs = 0.025, na.rm = TRUE),
                        Y_qsup95 = apply(postData$df_sim, 2, quantile, probs = 0.975, na.rm = TRUE))

  } else stop("type must be 'probability' (i.e., probability of survival) or 'number' (i.e., number of survivors)")


  plt_fit = df_plt %>%
    ggplot() + theme_bw() +
    theme(legend.position="none") +
    expand_limits(x = 0, y = 0) +
    labs(title = lab.main,
         x = lab.x,
         y = lab.y,
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
    geom_ribbon(aes(x = time,
                  ymin = Y_qinf95,
                  ymax = Y_qsup95 ,
                  group = profile), fill = "orange", alpha = 0.4, color = "grey90") +
    geom_line(aes(x = time,
                   y = Y_q50,
                   group = profile )) +
    geom_point(aes(x = time,
                  y = Y,
                  group = profile )) +
    facet_wrap(~ profile)

  return(plt_fit)
}

