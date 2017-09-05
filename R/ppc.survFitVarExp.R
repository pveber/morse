#' Posterior predictive check plot for \code{survFitVarExp} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitVarExp} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{survFit} objects.
#'
#' The black points show the observed number of survivors (on \eqn{X}-axis)
#'  against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise.
#'
#' @param x An object of class \code{survFitVarExp}
#' @param style graphical backend
#' @param \dots Further arguments to be passed to generic methods
#'
#' @export

ppc.survFitVarExp <- function(x,
                              xlab = "Observed nb of survivors",
                              ylab = "Predicted nb of survivors",
                              main = NULL, ...) {
  
  
  ### compute posteriors median and 95 CI
  jags.data <- x$jags.data
  df_ppc <- posteriorData(x)$df_ppc
  
  df_plt <- data_frame(Nsurv = jags.data$Nsurv,
                      time = jags.data$time,
                      replicate = jags.data$replicate,
                      Nsurv_q50 = apply(df_ppc, 2, quantile, probs = 0.5, na.rm = TRUE),
                      Nsurv_qinf95 = apply(df_ppc, 2, quantile, probs = 0.025, na.rm = TRUE),
                      Nsurv_qsup95 = apply(df_ppc, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
    mutate(col_range = ifelse(Nsurv_qinf95 > Nsurv | Nsurv_qsup95 < Nsurv, "out", "in"))
  
  
  ppc_plt <- df_plt %>%
    ggplot() + theme_bw() +
    theme(legend.position="none") +
    # expand_limits(x = 0, y = 0) +
    scale_colour_manual(values = c("green", "red")) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
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


#' @param x An object of class \code{survFitVarExp}

posteriorData = function(x){
  
  model_type = x$model_type
  mcmc = x$mcmc
  
  mctot <- do.call("rbind", mcmc)
  
  df_mctot = as_data_frame(mctot)
  
  df_psurv = select(df_mctot, contains("psurv"))
  df_ppc = select(df_mctot, contains("Nsurv_ppc"))
  
  ls_posteriorData = list( df_psurv = df_psurv,
                           df_ppc = df_ppc)
  
  return(ls_posteriorData)
}
