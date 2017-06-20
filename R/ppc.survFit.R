#' Posterior predictive check plot for \code{survFit} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitTKTD} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{gm_survFitTKTD} objects.
#'
#' The black points show the observed number of survivors (pooled
#' replicates, on \eqn{X}-axis) against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise. Samples with equal observed
#' value are shifted on the \eqn{X}-axis. For that reason, the
#' bisecting line (y = x), is represented by steps when observed
#' values are low. That way we ensure green intervals do intersect the
#' bisecting line.
#'
#' @param x An object of class \code{survFit}
#' @param style graphical backend
#' @param \dots Further arguments to be passed to generic methods
#'
#' @examples
#'
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
#' @export
ppc.survFit <- function(x, ...) {
  
  ### compute posteriors median and 95 CI
  modelData = x$modelData
  postData = posteriorData(x$mcmc, model_type = x$model_type)
  
  df_plt = data_frame(Nsurv = modelData$Nsurv,
                      time = modelData$time,
                      replicate = modelData$replicate,
                      Nsurv_q50 = apply(postData$df_ppc, 2, quantile, probs = 0.5, na.rm = TRUE),
                      Nsurv_qinf95 = apply(postData$df_ppc, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Nsurv_qsup95 = apply(postData$df_ppc, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
    mutate(col_range = ifelse(Nsurv_qinf95 > Nsurv | Nsurv_qsup95 < Nsurv, "out", "in"))
  
  
  ppc_plt = df_plt %>%
    ggplot() + theme_bw() +
    theme(legend.position="none") +
    # expand_limits(x = 0, y = 0) +
    scale_colour_manual(values = c("green", "red")) +
    scale_x_continuous(name="Observed nb of survivors") +
    scale_y_continuous(name="Predicted nb of survivors") +
    geom_abline(slope = 1) +
    geom_linerange(aes(x = Nsurv,
                       ymin = Nsurv_qinf95,
                       ymax = Nsurv_qsup95 ,
                       group = replicate,
                       color = col_range),
                   position = position_dodge(width=0.5)) +
    geom_point(aes(x = Nsurv,
                   y = Nsurv_q50,
                   group = replicate ),
               position = position_dodge(width=0.5))
  
  return(ppc_plt)
  
}

