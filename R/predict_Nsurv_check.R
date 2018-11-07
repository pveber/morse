#' Checking goodness-of-fit method for \code{survFitPredict} and
#'  \code{survFitPredict_Nsurv} objects
#' 
#' It returns measures of goodness-of-fit for prediction
#' 
#' @param object an object used to select a method \code{ppc_percent}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
#'
predict_Nsurv_check <- function(object, ...){
  UseMethod("predict_Nsurv_check")
}

#' Compute criteria to check model performance
#'
#' Provide various criteria for assessment of the model performance:
#' (i) percentage of observation within the 95% credible
#' interval of the Posterior Prediction Check (PPC), the Normalised Root Mean
#' Square Error (NRMSE) and the Survival Probability Prediction Error (SPPE)
#'
#' @param object an object of class \code{survFitPredict_Nsurv}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @return The function return a list with three items:
#' \item{Percent_PPC}{The criterion compares the predicted median numbers 
#' of survivors associated to their uncertainty limits with the observed numbers
#' of survivors. Based on experience, PPC resulting in less than \eqn{50\%} of the
#' observations within the uncertainty limits indicate poor model performance. A fit of
#' \eqn{100\%} may hide too large uncertainties of prediction (so covering all data).}
#' \item{NRMSE}{The criterion is based on the classical root-mean-square error (RMSE),
#'  used to aggregate the magnitudes of the errors in predictions for various time-points
#'  into a single measure of predictive power. In order to provide a criterion expressed
#'   as a percentage, NRMSE is the normalised RMSE by the mean of the observations.}
#' \item{SPPE}{The SPPE indicator is negative (between \eqn{0} and \eqn{-100\%} for an 
#' underestimation of effects, and positive (between \eqn{0} and \eqn{100}) for an 
#' overestimation of effects. An SPPE value of \eqn{0} means an exact prediction
#'  of the observed survival probability at the end of the experiment.}
#' 
#' @export
#'
predict_Nsurv_check.survFitPredict_Nsurv <- function(object, ...){
  
  df_ppc <- object$df_quantile %>%
    mutate(ppc_matching_check = ifelse(Nsurv_qinf95_check > Nsurv | Nsurv_qsup95_check < Nsurv, 0, 1),
           ppc_matching_valid = ifelse(Nsurv_qinf95_valid > Nsurv | Nsurv_qsup95_valid < Nsurv, 0, 1),
           rmse_id = (Nsurv - Nsurv_q50_valid)^2)
  
  percent_ppc_graphic <- sum(df_ppc$ppc_matching_check) / nrow(df_ppc) * 100
  percent_ppc_timeserie <- sum(df_ppc$ppc_matching_valid) / nrow(df_ppc) * 100
  
  # --- NRMSE
  nrmse <- mean(df_ppc$rmse_id) / mean(df_ppc$Nsurv) * 100
  ## version with distribution
  # rmse <- (jags.data$Nsurv - t(df_sim))^2
  # y_mean <- mean(jags.data$Nsurv)
  # nrmse <- rmse / y_mean 
  
  # --- SPPE
  
  df_sppe <- df_ppc %>%
    select(Nsurv, time, Nsurv_q50_valid, Nsurv_q50_check, replicate, time) %>%
    group_by(replicate) %>%
    arrange(replicate,time) %>%
    summarise(SPPE = (last(Nsurv) - last(Nsurv_q50_valid)) / first(Nsurv) * 100 )
    # summarise(sppe_check = (last(Nsurv) - last(Nsurv_q50_check)) / first(Nsurv) * 100,
    #           sppe_valid = (last(Nsurv) - last(Nsurv_q50_valid)) / first(Nsurv) * 100 )

  
  return( list(Percent_PPC = percent_ppc_timeserie,
               NRMSE = nrmse,
               SPPE = as.data.frame(df_sppe))
  )
  
}