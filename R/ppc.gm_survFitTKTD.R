#' Posterior predictive check plot for \code{gm_survFitTKTD} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{gm_survFitTKTD} class. It
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
#' @param x An object of class \code{survFitTKTD}
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
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
ppc.gm_survFitTKTD <- function(x, ...) {

  ### check class
  if (!is(x, "gm_survFitTKTD")){
    stop("x is not of class 'gm_survFitTKTD'!")
  }

  ### compute posteriors median and 95 CI
  modelData = x$modelData
  postData = posteriorData(x$mcmc, model_type = x$model_type)

  df_plt = data_frame(Nsurv = modelData$Nsurv,
                      time = modelData$time,
                      profile = modelData$profile,
                      Nsurv_q50 = apply(postData$df_ppc, 2, quantile, probs = 0.5, na.rm = TRUE),
                      Nsurv_qinf95 = apply(postData$df_ppc, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Nsurv_qsup95 = apply(postData$df_ppc, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
    mutate(col_range = ifelse(Nsurv_qinf95 > Nsurv | Nsurv_qsup95 < Nsurv, "out", "in"))


  ppc_plt = df_plt %>%
    ggplot() + theme_bw() +
    theme(legend.position="none") +
    expand_limits(x = 0, y = 0) +
    scale_colour_manual(values = c("green", "red")) +
    scale_x_continuous(name="Observed nb of survivors") +
    scale_y_continuous(name="Predicted nb of survivors") +
    geom_abline(slope = 1) +
    geom_linerange(aes(x = Nsurv,
                       ymin = Nsurv_qinf95,
                       ymax = Nsurv_qsup95 ,
                       group = profile,
                       color = col_range),
                   position = position_dodge(width=0.5)) +
    geom_point(aes(x = Nsurv,
                   y = Nsurv_q50,
                   group = profile ),
               position = position_dodge(width=0.5))

  return(ppc_plt)

}

# To retrive posterior object from a gm_survFitTKTD object
#

posteriorData = function(mcmc, model_type){
  mctot <- do.call("rbind", mcmc)

  df_mctot = as_data_frame(mctot)

  df_psurv = select(df_mctot, contains("psurv"))
  df_ppc = select(df_mctot, contains("Nsurv_ppc"))
  df_sim = select(df_mctot, contains("Nsurv_sim"))

  if(model_type == "SD"){

    df_param = select(df_mctot, c(kd_log10, hb_log10, kk_log10, z_log10))

  } else if(model_type == "IT"){

    df_param = select(df_mctot, c(kd_log10, hb_log10, alpha_log10, beta_log10))
  }

  ls_posteriorData = list( df_psurv = df_psurv,
                           df_ppc = df_ppc,
                           df_sim = df_sim,
                           df_param = df_param)

  return(ls_posteriorData)
}

