#' Posterior predictive check plot for \code{survFitPredict_Nsurv} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitPredict_Nsurv} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{survFitPredict_Nsurv} objects.
#' 
#' For \code{survFitPredict_Nsurv} object, PPC is based on times series simulated
#' for each replicate. In addition, the black points show the observed
#'  number of survivors (on \eqn{X}-axis)
#'  against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise.
#'
#' @param x An object of class \code{survFitPredict_Nsurv}
#' @param xlab A label for the \eqn{X}-axis, by default \code{Observed nb of survivors}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Predicted nb of survivors}.
#' @param main A main title for the plot.
#' @param \dots Further arguments to be passed to generic methods
#'
#' @export

ppc.survFitPredict_Nsurv <- function(x,
                              xlab = "Observed nb of survivors",
                              ylab = "Predicted nb of survivors",
                              main = NULL,
                              ...) {
  
  
  ### add color range "in" or "out"
  df_plt <- x$df_quantile %>%
    mutate(col_range = ifelse(Nsurv_qinf95_valid > Nsurv | Nsurv_qsup95_valid < Nsurv, "out", "in"))

  ppc_plt <- df_plt %>%
    ggplot() + theme_bw() +
    theme(legend.position="none") +
    expand_limits(x = 0, y = 0) +
    scale_colour_manual(values = c("green", "red")) +
    scale_x_continuous(name = xlab) +
    scale_y_continuous(name = ylab) +
    geom_abline(slope = 1) +
    geom_abline(slope = rep(1, 2 ), intercept = c(-12.5, 12.5), linetype = "dashed") +
    geom_abline(slope = rep(1, 2 ), intercept = c(-25, 25), linetype = "dotted") +
    geom_linerange(aes(x = Nsurv,
                       ymin = Nsurv_qinf95_valid,
                       ymax = Nsurv_qsup95_valid ,
                       group = replicate,
                       color = col_range),
                   position = position_dodge(width=0.5)) +
    geom_point(aes(x = Nsurv,
                   y = Nsurv_q50_valid,
                   group = replicate ),
               position = position_dodge(width=0.5))
  
  return(ppc_plt)
}


# @param x An object of class \code{survFitVarExp}

posteriorData <- function(x){
  
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
