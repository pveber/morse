#' Predict \eqn{x}\% Multiplication Factor at any specified time point for 
#' a \code{survFit} object.
#' 
#' The function \code{MFx}, \eqn{x}\% Multiplication Factor (\eqn{MF_x}), is use to compute
#'  the factor of multiplication required to kill \eqn{x}\% of the members of a tested population
#'  after a specified test duration (\code{time_MFx}) (default is the maximum
#'  time point of the experiment).
#'  
#'  Mathematical definition of \eqn{x}\% Multiplication Factor at time \eqn{t},
#'  denoted \eqn{MF(x,t)}, is:
#'  
#'  \eqn{S(MF(x,t) * C_w(\tau), t) = S( C_w(\tau), t)*(1- x/100)},
#'  
#'  where \eqn{S(MF(x,t)* C_w(\tau), t)} is the survival rate at concentration
#'  \eqn{MF(x,t)* C_w(\tau)} at time \eqn{t}.
#'   
#' 
#' @param object An object of class \code{survFit}
#' @param profile A dataframe with two columns \code{time} and \code{conc}.
#' @param x_MFx Percentage of survival reduction (e.g., \eqn{50} for \eqn{MF_{50}},
#'  \eqn{10} for \eqn{MF_{10}}, ...). The default is 50.
#' @param time_MFx A number giving the time at which  \eqn{MF_{x}} has to be estimated. 
#' If NULL, the latest time point of the experiment is used.
#' @param mcmc_size Can be used to reduce the number of mcmc samples in order to speed up
#'  the computation. The default is 1000.
#' @param MFx_max Upper bound of the multiplication factor. The default is 1000.
#' @param accuracy Accuracy of the multiplication factor. The default is 0.01.
#' @param quiet If \code{FALSE}, print the evolution of accuracy.
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns an object of class \code{MFx}, which is a list
#'  with the following information:
#' \item{Multiplication_factor}{The multiplication factor to obtain a median 
#' survial rate at time \code{time_MFx} of \eqn{X}\% lower than the median mortality
#' obtain from the input exposure profile, argument \code{profile} (\eqn{C_w(t)}).}
#' \item{time_MFx}{A number giving the time at which  \eqn{MF_{x}} has to be
#'  estimated as provided in arguments or if NULL, the latest time point of the
#'   experiment is used.}
#' \item{initial_x_MFx}{The median background mortality (i.e. \eqn{S(C_w(\tau), t)})}
#' \item{theoretical_x_MFx}{The theoretical mortality to reach after a \code{x_MFx}\% increase
#' of the mortality (i.e. \eqn{S(C_w(\tau), t)*(1- x/100)}).}
#' \item{final_x_MFx}{The mortality computed with the accuracy specified in the argument.}
#' \item{initial_prediction}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
#'  of the survival rate along time using the exposure profile without multiplication factor.}
#' \item{final_prediction}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
#'  of survival rate along time, leading to a increase of \code{x_MFx}\%  of mortality at
#'   time \code{time_MFx} compared to the provided exposure profile.}
#' 
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
#' # (5) estimate MF10 at time 4
#' MFx(out_SD, x_MFx = 50, time_MFx = 4)
#' }
#' 
#' 
#' @export
#' 
MFx.survFit <- function(object,
                        profile,
                        x_MFx = 50,
                        time_MFx = NULL,
                        mcmc_size = 1000,
                        MFx_max = 1000,
                        accuracy = 0.01,
                        quiet = FALSE, ...){
  
  ## Analyse profile data.frame
  if(!all(colnames(profile) %in% c("conc", "time")) || ncol(profile) != 2){
    stop("The argument 'profile' is a dataframe with two columns 'time' and 'conc'.")
  }
  
  ## Check time_MFx
  if(is.null(time_MFx))  time_MFx = max(profile$time)

  if(!(time_MFx %in% profile$time)){
    stop("Please provide a 'time_MFx' corresponding to a time-point at which concentration is provided.
            Interpolation of concentration is too specific to be automatized.")
  }
  
  profile$replicate <- rep("predict_MFx", nrow(profile))
  
  predict_init <- predict(object = object, data_predict = profile, spaghetti = FALSE, mcmc_size = mcmc_size)
  
  median_backgroundMortality_Conc0 = dplyr::filter(predict_init$df_quantile, time == time_MFx)$q50
  x_MFx_out = (100 - x_MFx) / 100 * median_backgroundMortality_Conc0
  
  profile_test = profile
  profile_test$conc = MFx_max * profile$conc
  predict_test <- predict(object = object, data_predict = profile_test, spaghetti = FALSE, mcmc_size = mcmc_size)
  median_Mortality_test = dplyr::filter(predict_test$df_quantile, time == time_MFx)$q50
  
  i_round = 0
  MFx_min = 0 ; MFx_test = MFx_max
  while(abs(x_MFx_out - median_Mortality_test) > accuracy){
    if(x_MFx_out - median_Mortality_test < 0){
      MFx_min = MFx_test
      MFx_test = MFx_test + (MFx_max - MFx_min)/2
    }
    if(x_MFx_out-median_Mortality_test > 0){
      MFx_max = MFx_test
      MFx_test = MFx_test - (MFx_max - MFx_min)/2
    }
    profile_test$conc = MFx_test * profile$conc
    predict_test <- predict(object = object, data_predict = profile_test, spaghetti = FALSE, mcmc_size = mcmc_size)
    median_Mortality_test = dplyr::filter(predict_test$df_quantile, time == time_MFx)$q50
    
    if(quiet == FALSE){
      i_round = i_round +1
      cat(i_round,"accuracy:", abs(x_MFx_out - median_Mortality_test), " with multiplication factor:",  MFx_test, "\n")
    }
  }
  
  ls_out = list(Multiplication_factor = MFx_test,
                x_MFx = x_MFx,
                time_MFx = time_MFx,
                initial_x_MFx = median_backgroundMortality_Conc0,
                theoretical_x_MFx = x_MFx_out,
                final_x_MFx = median_Mortality_test,
                initial_prediction = predict_init,
                final_prediction = predict_test)
  
  class(ls_out) = c("MFx", "list")
  
  return(ls_out)
}
