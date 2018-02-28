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
#' @param MFx_range A vector of length 2 with lower and upper bound of the 
#' range of multiplication factor \code{MFx}. The default is a vector \code{c(1, 1000)}.
#' @param npoints Number of points to compute within \code{MFx_range}. If \code{NULL}, only point computed 
#' during the adaptation using a binary search in \eqn{O(log n)} are taken into account.
#' @param mcmc_size Can be used to reduce the number of mcmc samples in order to speed up
#'  the computation. The default is 1000.
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into account from the posterior.
#' If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
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
                        MFx_range = c(1,1000),
                        npoints = NULL,
                        mcmc_size = 1000,
                        hb_value = TRUE,
                        spaghetti = FALSE,
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
  
  ls_profile <- list()
  ls_predict <- list()
  
  
  ls_profile[[1]] <- profile
  ls_predict[[1]] <- predict( object = object,
                              data_predict = ls_profile[[1]],
                              spaghetti = spaghetti,
                              mcmc_size = mcmc_size,
                              hb_value = hb_value)
  
  filter_time_MFx = dplyr::filter(ls_predict[[1]]$df_quantile, time == time_MFx)

  median_Mortality_test <- filter_time_MFx$q50
  theoretical_x_MFx <- (100 - x_MFx) / 100 * filter_time_MFx$q50 # Necessary to compared with accuracy
  
  if(is.null(npoints)){
    #
    # binary search of MFx in O(log n)
    #
    i = 1
    MFx = 1
    MFx_min = min(MFx_range)
    MFx_test = max(MFx_range)
    
    while(abs(theoretical_x_MFx - median_Mortality_test) > accuracy){
    
      MFx = c(MFx, MFx_test)
        
      ls_profile[[i+1]] <- profile
      ls_profile[[i+1]]$conc <- MFx_test * profile$conc
      
      ls_predict[[i+1]] <- predict(object = object,
                                  data_predict = ls_profile[[i+1]],
                                  spaghetti = spaghetti,
                                  mcmc_size = mcmc_size,
                                  hb_value = hb_value)
      
      filter_time_MFx = dplyr::filter(ls_predict[[i+1]]$df_quantile, time == time_MFx)
      median_Mortality_test = filter_time_MFx$q50
      df_MFx <- dplyr::add_row(df_MFx,
                               MFx = MFx_test,
                               q50 = filter_time_MFx$q50,
                               qinf95 = filter_time_MFx$qinf95,
                               qsup95 = filter_time_MFx$qsup95 )
      
      
      
      if(quiet == FALSE){
        cat(i,"accuracy:", abs(theoretical_x_MFx - median_Mortality_test), " with multiplication factor:",  MFx_test, "\n")
      }
      
      i = i + 1
      if(theoretical_x_MFx - median_Mortality_test < 0){
        MFx_min = MFx_test
        MFx_test = MFx_test + (MFx_max - MFx_min)/2
      }
      if(theoretical_x_MFx - median_Mortality_test > 0){
        MFx_max = MFx_test
        MFx_test = MFx_test - (MFx_max - MFx_min)/2
      }
    }
    k <- 1:length(MFx)
  }
  if(!is.null(npoints)){
    
    MFx = seq(min(MFx_range), max(MFx_range), length.out = npoints)
    
    k <- 1:npoints

    ls_profile <- lapply(k, function(kit){
      profil_test <- profile
      profil_test$conc <- MFx[kit] * profile$conc
      return(profil_test)
    })
    
    ls_predict <- lapply(k, function(kit){
      predict(object = object,
              data_predict = ls_profile[[kit]],
              spaghetti = spaghetti,
              mcmc_size = mcmc_size,
              hb_value = hb_value)
    })
    
  }
  
  ls_predict_quantile <- lapply(k, function(kit){
    df_quantile <- ls_predict[[kit]]$df_quantile
    df_quantile$MFx <- rep(MFx[[kit]], nrow(ls_predict[[kit]]$df_quantile))
    return(df_quantile)
  })
  
  df_ls_predict <- do.call("rbind", ls_predict_quantile)
  
  df_MFx <- dplyr::filter(df_ls_predict, time == time_MFx) %>%
    mutate(accuracy = abs(theoretical_x_MFx - q50)) %>%
    mutate(optimal_MFx = ifelse(accuracy == min(accuracy), MFx, NA))

  optimal_MFx <- dplyr::filter(df_MFx, !is.na(optimal_MFx))
  
  ls_out = list(optimal_MFx = optimal_MFx,
                x_MFx = x_MFx,
                time_MFx = time_MFx,
                initial_x_MFx = df_MFx$q50,
                theoretical_x_MFx = theoretical_x_MFx,
                df_MFx = df_MFx, # MFx at specific time
                ls_predict = ls_predict,
                df_ls_predict = df_ls_predict,
                ls_profile = ls_profile)
  
  class(ls_out) = c("list", "MFx")
  
  return(ls_out)
}


