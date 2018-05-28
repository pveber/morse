#' Predict x\% Multiplication Factor at any specified time point for 
#' a \code{survFit} object.
#' 
#' The function \code{MFx}, \eqn{x}\% Multiplication Factor at time \eqn{t}, (\eqn{MF(x,t)}),
#' is used to compute the multiplication factor
#' applied to the concentration exposure profile in order to
#' reduce by \eqn{x}\% (argument \code{X}) the survival rate at a
#'  specified test duration \eqn{t} (argument \code{time_MFx}) (default is the maximum
#'  time point of the experiment).
#'  
#'  Mathematical definition of \eqn{x}\% Multiplication Factor at time \eqn{t}
#'  (at the end of a time series \eqn{T = \{0, \dots, t\}}),
#'  denoted \eqn{MF(x,t)}, is given by:
#'  
#'  \eqn{S(MF(x,t) * C_w(\tau \in T), t) = S( C_w(\tau \in T), t)*(1- x/100)},
#'  
#'  where \eqn{C_w(\tau \in T)} is the initial exposure profile without
#'  multiplication factor. And so the expression \eqn{S(MF(x,t)* C_w(\tau \in T), t)}
#'   is the survival rate after an exposure profile
#'    \eqn{MF(x,t)* C_w(\tau \in T)} at time \eqn{t}.
#'   
#' 
#' @param object An object of class \code{survFit}
#' @param data_predict A dataframe with two columns \code{time} and \code{conc}.
#' @param X Percentage of survival change (e.g., \eqn{50} for survival decrease of 50\%
#'  , or \eqn{-50} for survival increase of 50\%).The default is 50. 
#'  Only time series computed during the adaptation using a binary search in
#'  \eqn{O(log(n))} are returned. However, if \code{NULL}, all time series
#'  computed from the vector \code{MFx_range} are returned.
#' @param time_MFx A number giving the time at which  \eqn{MF(x,t)} has to be estimated. 
#' If NULL, the latest time point of the profile is used.
#' @param MFx_range A vector from which lower and upper bound of the range of the
#'  multiplication factor \code{MFx} are generated. The default is a vector \code{c(0, 1000)}.
#' If argument \code{X} is \code{NULL}, then all the time series generated with
#' \code{MFx_range} are returned.
#' @param mcmc_size Can be used to reduce the number of MCMC samples in order to speed up
#'  the computation. The default is 1000.
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into
#'  account from the posterior.
#' If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param accuracy Accuracy of the multiplication factor. The default is 0.01.
#' @param quiet If \code{FALSE}, print the evolution of accuracy.
#' @param threshold_iter Threshold number of iteration.
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns an object of class \code{MFx}, which is a list
#'  with the following information:
#'  \item{X_prop}{Survival rate for \code{X} percent of reduction of the initial median 
#' survival rate at time \code{time_MFx}.}
#' \item{X_prop_provided}{A number giving the proportion of reduction in survival.}
#' \item{time_MFx}{A number giving the time at which  \eqn{MF(x,t)} has to be
#'  estimated as provided in arguments or if NULL, the latest time point of the
#'   profile is used.}
#' \item{df_MFx}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
#'  of \eqn{MF(x,t)} at time \eqn{t}, \code{time_MFx}, for \eqn{x}\% of survival reduction.}
#' \item{df_dose}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
#'  of survival rate along the computed multiplication factor and at time \code{time_MFx}.}
#' \item{MFx_tested}{A vector of all multiplication factors computed.} 
#' \item{ls_predict}{A list of all object of class \code{survFitPredict} obtained
#' from computing survival rate for every profiles build from the vector of
#' multiplication factors \code{MFx_tested}.}
#' 
#'    
#' @examples 
#' 
#' # (1) Load the data
#' data("propiconazole")
#' 
#' # (2) Create an object of class 'survData'
#' dataset <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFit function with model_type SD (or IT)
#' out_SD <- survFit(dataset, model_type = "SD")
#' 
#' # (4) data to predict
#' data_4prediction <- data.frame(time = 1:10, conc = c(0,0.5,3,3,0,0,0.5,3,1.5,0))
#' 
#' # (5) estimate MF for 30% reduction of survival at time 4
#' MFx_SD_10.4 <- MFx(out_SD, data_predict = data_4prediction , X = 30, time_MFx = 4)
#' }
#' 
#' 
#' @export
#' 
MFx.survFit <- function(object,
                        data_predict,
                        X = 50,
                        time_MFx = NULL,
                        MFx_range = c(0,1000),
                        mcmc_size = 1000,
                        hb_value = TRUE,
                        spaghetti = FALSE,
                        accuracy = 0.01,
                        quiet = FALSE,
                        threshold_iter = 100,
                        ...){
  
  ## Analyse data_predict data.frame
  if(!all(colnames(data_predict) %in% c("conc", "time")) || ncol(data_predict) != 2){
    stop("The argument 'data_predict' is a dataframe with two columns 'time' and 'conc'.")
  }
  
  ## Check time_MFx
  if(is.null(time_MFx))  time_MFx = max(data_predict$time)

  if(!(time_MFx %in% data_predict$time)){
    stop("Please provide a 'time_MFx' corresponding to a time-point at which concentration is provided.
            Interpolation of concentration is too specific to be automatized.")
  }
  
  ls_data_predict <- list()
  ls_predict <- list()

  ls_data_predict[[1]] <- data_predict
  ls_data_predict[[1]]$replicate <- rep("predict_MFx_1", nrow(data_predict))
  
  ls_predict[[1]] <- predict( object = object,
                              data_predict = ls_data_predict[[1]],
                              spaghetti = spaghetti,
                              mcmc_size = mcmc_size,
                              hb_value = hb_value)

  filter_time_MFx = dplyr::filter(ls_predict[[1]]$df_quantile, time == time_MFx)

  median_Mortality_test <- filter_time_MFx$q50
  theoretical_X <- (100 - X) / 100 * filter_time_MFx$q50 # Necessary to compared with accuracy

  if(!is.null(X)){
    binarySearch_MFx_q50 <- binarySearch_MFx(object = object,
                                             spaghetti = spaghetti,
                                             mcmc_size = mcmc_size,
                                             hb_value = hb_value,
                                             MFx_range = MFx_range,
                                             time_MFx = time_MFx,
                                             theoretical_X = theoretical_X,
                                             value_mortality_test = median_Mortality_test,
                                             accuracy = accuracy,
                                             data_predict = data_predict,
                                             ls_data_predict = ls_data_predict,
                                             ls_predict = ls_predict,
                                             quiet = quiet,
                                             quantile = "q50",
                                             threshold_iter = threshold_iter) # "q50", "qinf95", "qsup95"
    binarySearch_MFx_qinf95 <- binarySearch_MFx(object = object,
                                                spaghetti = spaghetti,
                                                mcmc_size = mcmc_size,
                                                hb_value = hb_value,
                                                MFx_range = MFx_range,
                                                time_MFx = time_MFx,
                                               theoretical_X = theoretical_X,
                                               value_mortality_test = filter_time_MFx$qinf95,
                                               accuracy = accuracy,
                                               data_predict = data_predict,
                                               ls_data_predict = ls_data_predict,
                                               ls_predict = ls_predict,
                                               quiet = quiet,
                                               quantile = "qinf95",
                                               threshold_iter = threshold_iter) # "q50", "qinf95", "qsup95"
    binarySearch_MFx_qsup95 <- binarySearch_MFx(object = object,
                                                spaghetti = spaghetti,
                                                mcmc_size = mcmc_size,
                                                hb_value = hb_value,
                                                MFx_range = MFx_range,
                                                time_MFx = time_MFx,
                                               theoretical_X = theoretical_X,
                                               value_mortality_test = filter_time_MFx$qsup95,
                                               accuracy = accuracy,
                                               data_predict = data_predict,
                                               ls_data_predict = ls_data_predict,
                                               ls_predict = ls_predict,
                                               quiet = quiet,
                                               quantile = "qsup95",
                                               threshold_iter = threshold_iter) # "q50", "qinf95", "qsup95"
    
    #
    # Make a dataframe with quantile of all generated time series
    #
    
    ls_predict_quantile_q50 <- lapply(binarySearch_MFx_q50$k, function(kit){
      df_quantile <- binarySearch_MFx_q50$ls_predict[[kit]]$df_quantile
      df_quantile$MFx <- rep(binarySearch_MFx_q50$MFx[[kit]], nrow(binarySearch_MFx_q50$ls_predict[[kit]]$df_quantile))
      return(df_quantile)
    })
    ls_predict_quantile_qinf95 <- lapply(binarySearch_MFx_qinf95$k, function(kit){
      df_quantile <- binarySearch_MFx_qinf95$ls_predict[[kit]]$df_quantile
      df_quantile$MFx <- rep(binarySearch_MFx_qinf95$MFx[[kit]], nrow(binarySearch_MFx_qinf95$ls_predict[[kit]]$df_quantile))
      return(df_quantile)
    })
    ls_predict_quantile_qsup95 <- lapply(binarySearch_MFx_qsup95$k, function(kit){
      df_quantile <- binarySearch_MFx_qsup95$ls_predict[[kit]]$df_quantile
      df_quantile$MFx <- rep(binarySearch_MFx_qsup95$MFx[[kit]], nrow(binarySearch_MFx_qsup95$ls_predict[[kit]]$df_quantile))
      return(df_quantile)
    })
    
    predict_MFx_quantile_q50 <- do.call("rbind", ls_predict_quantile_q50)
    predict_MFx_quantile_qinf95 <- do.call("rbind", ls_predict_quantile_qinf95)
    predict_MFx_quantile_qsup95 <- do.call("rbind", ls_predict_quantile_qsup95)
    #
    # doseResponse dataframe at specific time_MFx
    #
    df_dose_q50 <- dplyr::filter(predict_MFx_quantile_q50, time == time_MFx)
    df_dose_q50$id = rep("q50", nrow(df_dose_q50))
    df_dose_qinf95 <- dplyr::filter(predict_MFx_quantile_qinf95, time == time_MFx)
    df_dose_qinf95$id = rep("qinf95", nrow(df_dose_qinf95))
    df_dose_qsup95 <- dplyr::filter(predict_MFx_quantile_qsup95, time == time_MFx)
    df_dose_qsup95$id = rep("qsup95", nrow(df_dose_qsup95))
    
    ## Additional element to return
    df_dose <- do.call("rbind", list(df_dose_q50, df_dose_qinf95, df_dose_qsup95))
    MFx <- binarySearch_MFx_q50$MFx
    ls_predict <- binarySearch_MFx_q50$ls_predict
    
  }
  if(is.null(X)){
    theoretical_X = NULL # to return in the final object
    
    MFx <- MFx_range
    
    k <- 1:length(MFx_range)

    ls_data_predict <- lapply(k, function(kit){
      profil_test <- data_predict
      profil_test$conc <- MFx[kit] * data_predict$conc
      profil_test$replicate <- rep(paste0("predict_MFx_", MFx[kit]), nrow(data_predict))
      return(profil_test)
    })
    
    ls_predict <- lapply(k, function(kit){
      predict(object = object,
              data_predict = ls_data_predict[[kit]],
              spaghetti = spaghetti,
              mcmc_size = mcmc_size,
              hb_value = hb_value)
    })
    
    #
    # Make a dataframe with quantile of all generated time series
    #
    
    ls_predict_quantile <- lapply(k, function(kit){
      df_quantile <- ls_predict[[kit]]$df_quantile
      df_quantile$MFx <- rep(MFx[[kit]], nrow(ls_predict[[kit]]$df_quantile))
      return(df_quantile)
    })
    predict_MFx_quantile <- do.call("rbind", ls_predict_quantile)
    
    df_dose <- dplyr::filter(predict_MFx_quantile, time == time_MFx)
    
  }
  
  #
  # Compute table with the optimal MFx obtained if X != NULL
  #
  if(!is.null(X)){
    
    MFx_q50 = df_dose_q50$MFx[nrow(df_dose_q50)]
    MFx_qinf95 = df_dose_qinf95$MFx[nrow(df_dose_qinf95)]
    MFx_qsup95 = df_dose_qsup95$MFx[nrow(df_dose_qsup95)]
    # Compute MFx q95:
    # pts_MFx <- pointsMFx(df_dose, median_Mortality_test)
    # MFx_qinf95 <- pts_MFx$MFx_qinf95
    # MFx_qsup95 <- pts_MFx$MFx_qsup95
    # 
    # Return dataframe of quantiles MFx
    
    df_MFx <- data.frame(quantile = c("median", "quantile 2.5%", "quantile 97.5%"),
                         MFx = c(MFx_q50, MFx_qinf95, MFx_qsup95))
  } else{
    df_MFx <- data.frame(quantile = c("median", "quantile 2.5%", "quantile 97.5%"),
                         MFx = c(NA, NA, NA))
  }
  
# warning("This is not an error message:
# Just take into account that MFx as been estimated with a binary
# search using the 'accuracy' argument. To improve the shape of the curve, you
# can use X = NULL, and computed time series around the median MFx, with the
#           vector `MFx_range`.")
  
  ls_out = list(X_prop = theoretical_X,
                X_prop_provided = X/100,
                time_MFx = time_MFx,
                df_MFx = df_MFx,
                df_dose = df_dose, # return MFx at specific time
                MFx_tested = MFx,
                ls_predict = ls_predict)
  
  class(ls_out) = c("list", "MFx")
  
  return(ls_out)
}


# points for LCx
# 
# 
# pointsMFx <- function(df_dose, X_prop){
#   
#   if(min(df_dose$qinf95) < X_prop & X_prop < max(df_dose$qinf95)){
#     df.qinf95 <- select(df_dose, c(MFx, qinf95))%>%
#       dplyr::add_row(qinf95 = X_prop)%>%
#       dplyr::arrange(qinf95)%>%
#       dplyr::mutate(MFx = na.approx(MFx, qinf95, na.rm = FALSE))%>%
#       dplyr::filter(qinf95 == X_prop)
#     
#     MFx_qinf95 <- df.qinf95$MFx
#     
#   } else {
#     MFx_qinf95 <- NA
#     
#     warning(paste("No 95%inf for survival rate of", X_prop ,
#                   " in the range of multiplication factors under consideration: [",
#                   min(df_dose$MFx), ";", max(df_dose$MFx), "]"))
#   }
#   
#   if(min(df_dose$qsup95) < X_prop & X_prop < max(df_dose$qsup95)){
#     df.qsup95 <- select(df_dose, c(MFx,qsup95)) %>%
#       add_row(qsup95 = X_prop) %>%
#       arrange(qsup95) %>%
#       mutate(MFx = na.approx(MFx,qsup95, na.rm = FALSE)) %>%
#       filter(qsup95 == X_prop)
#     
#     MFx_qsup95 <- df.qsup95$MFx
#     
#   } else {
#     
#     MFx_qsup95 <- NA
#     warning(paste("No 95%sup for survival rate of", X_prop,
#                   " in the range of multiplication factors under consideration: [",
#                   min(df_dose$MFx), ";", max(df_dose$MFx), "]"))
#   }
#   
#   return(list(MFx_qinf95 = MFx_qinf95,
#               MFx_qsup95 = MFx_qsup95))
# }



##########################
#
#
#

#
# binary search of MFx in O(log n)
#

binarySearch_MFx <- function(object,
                             spaghetti,
                             mcmc_size,
                             hb_value,
                             MFx_range,
                             time_MFx,
                             theoretical_X,
                             value_mortality_test,
                             accuracy,
                             data_predict,
                             ls_data_predict,
                             ls_predict,
                             quiet,
                             quantile, # "q50", "qinf95", "qsup95"
                             threshold_iter
                             ){
    #
    # binary search of MFx in O(log n)
    #
    i = 1
    MFx = 1
    MFx_min = min(MFx_range)
    MFx_max = max(MFx_range)
    MFx_test = max(MFx_range)
    
    while(abs(theoretical_X - value_mortality_test) > accuracy){
      
      MFx = c(MFx, MFx_test)
      
      ls_data_predict[[i+1]] <- data_predict
      ls_data_predict[[i+1]]$conc <- MFx_test * data_predict$conc
      ls_data_predict[[i+1]]$replicate <- rep(paste0("predict_MFx_", MFx_test), nrow(data_predict))
      
      ls_predict[[i+1]] <- predict(object = object,
                                   data_predict = ls_data_predict[[i+1]],
                                   spaghetti = spaghetti,
                                   mcmc_size = mcmc_size,
                                   hb_value = hb_value)
      
      filter_time_MFx = dplyr::filter(ls_predict[[i+1]]$df_quantile, time == time_MFx)
      if(quantile == "q50"){ value_mortality_test = filter_time_MFx$q50 }
      if(quantile == "qinf95"){ value_mortality_test = filter_time_MFx$qinf95 }
      if(quantile == "qsup95"){ value_mortality_test = filter_time_MFx$qsup95 }
      
      if(quiet == FALSE){
        cat(quantile, i,"accuracy:", abs(theoretical_X - value_mortality_test), " with multiplication factor:",  MFx_test, "\n")
      }
      
      i = i + 1
      if(theoretical_X - value_mortality_test < 0){
        MFx_min = MFx_test
        MFx_test = MFx_test + (MFx_max - MFx_min)/2
      }
      if(theoretical_X - value_mortality_test > 0){
        MFx_max = MFx_test
        MFx_test = MFx_test - (MFx_max - MFx_min)/2
      }
      if(MFx_test == max(MFx_range)){
        MFx_test <- NULL
        warning(paste("For", quantile, "The multiplication factor is over the bound of", max(MFx_range)))
        break
      }
      if(i > threshold_iter){
        MFx_test <- NULL
        warning(paste("For", quantile, "the number of iteration reached the threshold number of iteration."))
        break
      }
    }
    k <- 1:length(MFx)
    
    return(list(k = k,
                MFx = MFx,
                ls_predict = ls_predict,
                ls_data_predict = ls_data_predict))
  }
  
  
  
  
