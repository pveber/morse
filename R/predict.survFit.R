#' Predict method for \code{survFit} objects
#'
#' This is the generic \code{predict} S3 method for the \code{survFit} class.
#' It provides simulation for "SD" or "IT" models under constant or time-variable exposure.
#'
#' @param x A object of class \code{survFit}
#' @param data_predict A dataframe with three columns \code{time}, \code{conc} and \code{replicate}
#'  used for prediction. If \code{NULL}, prediction is based on \code{x} object of 
#'  class \code{survFit} used for fitting.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param mcmc_size A numerical value refering by default to the size of the mcmc in object \code{survFit}.
#'  \code{mcmc_size} can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#' 
#' @examples 
#'
#' # (1) Load the survival data
#' data("propiconazole_pulse_exposure")
#'
#' # (2) Create an object of class "survData"
#' dataset <- survData(propiconazole_pulse_exposure)
#'
#' \dontrun{
#' # (3) Run the survFit function
#' out <- survFit(dataset , model_type = "SD")
#'
#' # (4) Create a new data table for prediction
#' data_4prediction <- data.frame(time = 1:10, conc = c(0,5,30,30,0,0,5,30,15,0), replicate= rep("predict", 10))
#'
#' # (5) Predict on a new data set
#' predict_out <- predict(x = out, data_predict = data_4prediction, spaghetti = TRUE)
#'
#' }
#' 
#' 
#' @export
#'
predict.survFit <- function(x,
                            data_predict = NULL,
                            spaghetti = FALSE,
                            mcmc_size = NULL){
  
  # Initialisation
  mcmc <- x$mcmc
  model_type <- x$model_type
  extend_time <- 100
  
  if(is.null(data_predict)){
    if("survFitVarExp" %in% class(x)){
      x_interpolate = data.frame(
        time = x$jags.data$time_long,
        conc = x$jags.data$conc_long,
        replicate = x$jags.data$replicate_long)
    } else{
      data_predict = data.frame(
        time = x$jags.data$time,
        conc = x$jags.data$conc,
        replicate = x$jags.data$replicate)
      
      x_interpolate <- predict_interpolate(data_predict,  extend_time = extend_time) %>%
        dplyr::arrange(replicate, time)
    }
  }
  if(!is.null(data_predict)){
    x_interpolate <- predict_interpolate(data_predict,  extend_time = extend_time) %>%
      dplyr::arrange(replicate, time)
  }
    
  df <- data.frame(
    time = x_interpolate$time,
    conc = x_interpolate$conc,
    replicate = x_interpolate$replicate)

  unique_replicate <- unique(df$replicate)
  
  ls_time <- list()
  ls_conc <- list()
  
  for(i in 1:length(unique_replicate)){
    
    ls_time[[i]] <- dplyr::filter(df, replicate == unique_replicate[i])$time
    ls_conc[[i]] <- dplyr::filter(df, replicate == unique_replicate[i])$conc
    
  }
  
  # ------- Computing
  
  mcmc.samples = mcmc
  
  if(!is.null(mcmc_size)){
    reduc_tab = lapply(mcmc.samples, "[", 
                       seq(1, nrow(mcmc.samples[[1]]), length = mcmc_size),
                       1:ncol(mcmc.samples[[1]]))
    mcmc.samples = reduc_tab
  }
  
  mctot = do.call("rbind", mcmc.samples)
  kd = 10^mctot[, "kd_log10"]
  hb <- 10^mctot[, "hb_log10"]
  
  k = 1:length(unique_replicate)
  
  
  if(model_type == "SD"){
    kk <- 10^mctot[, "kk_log10"]
    z <- 10^mctot[, "z_log10"]
    
    dtheo = lapply(k, function(kit) { # Pour chaque replicat
      Surv.SD_Cext(Cw = ls_conc[[kit]],
                   time = ls_time[[kit]],
                   kk=kk,
                   kd=kd,
                   hb=hb,
                   z=z)
    })
    
  }
  if(model_type == "IT"){
    
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    dtheo = lapply(k, function(kit) { # Pour chaque replicat
      Surv.IT_Cext (Cw = ls_conc[[kit]],
                    time = ls_time[[kit]],
                    kd = kd,
                    hb = hb,
                    alpha = alpha,
                    beta = beta)
    })
    
  }
  
  #transpose
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  df_quantile = dplyr::data_frame(
    time = df$time,
    conc = df$conc,
    replicate = df$replicate,
    q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
    qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE),
    qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  if(spaghetti == TRUE){
    random_column <- sample(1:ncol(dtheo), size = round(10/100 * ncol(dtheo)))
    df_spaghetti <- as_data_frame(dtheo[, random_column]) %>%
      mutate(time = df$time,
             conc = df$conc,
             replicate = df$replicate)
  } else df_spaghetti <- NULL
  
  
  return_object <- list(df_quantile = df_quantile,
                        df_spaghetti = df_spaghetti)
    
  class(return_object) <- c(class(return_object), "survFitPredict")

  return(return_object)
  
}

# Survival function for "IT" model with external concentration changing with time
#
# @param Cw A scalar of external concentration
# @param time A vector of time
# @param kk a vector of parameter
# @param kd a vector of parameter
# @param z a vector of parameter
# @param hb a vector of parameter
# 
#
# @return A matrix generate with coda.samples() function
#

Surv.SD_Cext <- function(Cw, time, kk, kd, z, hb){
  
  time.prec = dplyr::lag(time, 1) ; time.prec[1] = time[1] #   time[1] = tprec[1]
  diff.int = (exp(time %*% t(kd)) * Cw + exp(time.prec %*% t(kd)) * Cw )/2 * (time-time.prec) #OK time[1]-tprec[1] = 0
  
  D = kd * exp(-kd %*% t(time)) * t(apply(diff.int,2,cumsum))
  
  lambda = kk * pmax(D-z,0) + hb # the pmax function is important here for elementwise maximum with 0 and D[i,j]-z ATTENTION: pmax(0,D) != pmax(D,0)
  
  lambda.prec = dplyr::lag(lambda, 1 ) ; lambda.prec[1] = lambda[1]
  
  int.lambda =  t(t((lambda + lambda.prec)/2) * (time-time.prec))
  
  S <- exp(-t(apply(int.lambda,1,cumsum)))
  
  return(S)
}

# Survival function for "IT" model with external concentration changing with time
#
# @param Cw A scalar of external concentration
# @param time A vector of time
# @param kk a vector of parameter
# @param kd a vector of parameter
# @param z a vector of parameter
# @param hb a vector of parameter
# 
#
# @return A matrix generate with coda.samples() function
#

Surv.IT_Cext <- function(Cw, time, kd, hb, alpha, beta)
{
  time.prec = dplyr::lag(time, 1) ; time.prec[1] = time[1] #   time[1] = tprec[1]
  diff.int = (exp(time %*% t(kd)) * Cw + exp(time.prec %*% t(kd)) * Cw )/2 * (time-time.prec) #OK time[1]-tprec[1] = 0
  
  D <- kd * exp(-kd %*% t(time)) * t(apply(diff.int,2,cumsum))
  
  D.max <- t(apply(D,1,cummax))
  
  S <-  exp(-hb %*% t(time))*(1-plogis(log(D.max),location=log(alpha),scale=1/beta))
  
  return(S)
}


# Create a dataset for survival analysis when the replicate of concentration is variable
#
# @param x An object of class \code{survData}
#
# @return A dataframe
#
predict_interpolate <- function(x, extend_time = 100){
  
  ## data.frame with time
  
  df_MinMax <- x %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(min_time = min(time, na.rm = TRUE),
                     max_time = max(time, na.rm = TRUE)) %>%
    dplyr::group_by(replicate) %>%
    # dplyr::do(data.frame(replicate = .$replicate, time = seq(.$min_time, .$max_time, length = extend_time)))
    dplyr::do(data_frame(replicate = .$replicate, time = seq(.$min_time, .$max_time, length = extend_time)))
  
  x_interpolate <- dplyr::full_join(df_MinMax, x,
                                    by = c("replicate", "time")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>% # organize in replicate and time
    dplyr::mutate(conc = zoo::na.approx(conc, time, na.rm = FALSE)) %>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    dplyr::mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc) )
  
  return(x_interpolate)
}

