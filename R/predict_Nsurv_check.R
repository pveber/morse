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
#' (i) percentage percentage of observation within the 95% credible
#' interval of the posterior prediction check (PPC), the Normalised Root Mean
#' Square Error (NRMSE) and the survival-probability prediction error (SPPE)
#'
#' @param object an object of class \code{survFitPredict_Nsurv}
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
  nrmse <- mean(df_ppc$rmse_id) / mean(df_ppc$Nsurv)
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