#' Predict method for \code{survFit} objects
#'
#' This is the generic \code{predict} S3 method for the \code{survFit} class.
#' It provides simulation for "SD" or "IT" models under constant or time-variable exposure.
#'
#' @param object An object of class \code{survFit}
#' @param data_predict A dataframe with three columns \code{time}, \code{conc} and \code{replicate}
#'  used for prediction. If \code{NULL}, prediction is based on \code{x} object of 
#'  class \code{survFit} used for fitting.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param mcmc_size Can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into account from the posterior.
#' If \code{FALSE}, parameter \code{hb} is set to 0. The default is \code{TRUE}.
#' @param \dots Further arguments to be passed to generic methods
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
#' data_4prediction <- data.frame(time = 1:10,
#'                                conc = c(0,5,30,30,0,0,5,30,15,0),
#'                                replicate= rep("predict", 10))
#'
#' # (5) Predict on a new dataset
#' predict_out <- predict(object = out, data_predict = data_4prediction, spaghetti = TRUE)
#'
#' }
#' 
#' @import deSolve
#' 
#' @export
#'
predict_ode.survFit <- function(object,
                                data_predict = NULL,
                                spaghetti = FALSE,
                                mcmc_size = NULL,
                                hb_value = TRUE,
                                interpolate_length = NULL,
                                interpolate_method = "linear",
                                na.rm = TRUE,
                                ...) {
  x <- object # Renaming to satisfy CRAN checks on S3 methods
  # arguments should be named the same when declaring a
  # method and its instantiations
  
  # Initialisation
  mcmc <- x$mcmc
  model_type <- x$model_type
  
  if(is.null(data_predict)){
    x_interpolate = data.frame(
        time = x$jags.data$time,
        conc = x$jags.data$conc,
        replicate = x$jags.data$replicate)
  }
  if(!is.null(data_predict)){
    x_interpolate <- data_predict
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
  if(is.null(mcmc_size)){
    mcmc_size = nrow(mctot)
  }
  
  kd = 10^mctot[, "kd_log10"]
  
  if(hb_value == TRUE){
    hb <- 10^mctot[, "hb_log10"]
  } else if(hb_value == FALSE){
    hb <- rep(0, nrow(mctot))
  }
  
  k = 1:length(unique_replicate)
  
  if(model_type == "SD"){
    kk <- 10^mctot[, "kk_log10"]
    z <- 10^mctot[, "z_log10"]
    
    dtheo = lapply(k, function(kit) { # For each replicate
      SurvSD_ode(Cw = ls_conc[[kit]],
                 time = ls_time[[kit]],
                 kk=kk,
                 kd=kd,
                 hb=hb,
                 z=z,
                 mcmc_size = mcmc_size,
                 interpolate_length = interpolate_length,
                 interpolate_method = interpolate_method)
    })
    
  }
  if(model_type == "IT"){
    
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    dtheo = lapply(k, function(kit) { # For each replicate
      SurvIT_ode(Cw = ls_conc[[kit]],
                 time = ls_time[[kit]],
                 kd = kd,
                 hb = hb,
                 alpha = alpha,
                 beta = beta,
                 mcmc_size = mcmc_size,
                 interpolate_length = interpolate_length,
                 interpolate_method = interpolate_method)
    })
    
  }
  
  # Transpose
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  df_quantile = dplyr::data_frame(
    time = df$time,
    conc = df$conc,
    replicate = df$replicate,
    q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = na.rm),
    qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = na.rm),
    qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = na.rm)
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

SurvSD_ode <- function(Cw, time, kk, kd, z, hb, mcmc_size = NULL, interpolate_length = NULL, interpolate_method = "linear"){
  
  ## external signal with several rectangle impulses
  signal <- data.frame(times = time, 
                       import = Cw)
  
  sigimp <- approxfun(signal$times, signal$import, method = interpolate_method, rule = 2)
  
  if(!is.null(interpolate_length)){
    times <- seq(min(time), max(time), length = interpolate_length)
  } else{
    times <- signal$times
  }
  
  ## The parameters
  parms  <- list( kd = kd,
                  kk = kk,
                  z = z,
                  mcmc_size = mcmc_size)
  
  ## Start values for steady state
  xstart <- c(D = rep(0, mcmc_size),
              H = rep(0, mcmc_size))
  
  ## Solve model
  out <- ode(y = xstart,
             times = times,
             func = model_SD,
             parms,
             input = sigimp)
  
  mat_4cast_H = out %>%
    as.data.frame() %>%
    select(contains("H"),-c(time,signal)) %>%
    as.matrix()
  
  S <- exp(-t(mat_4cast_H))
  dtheo = t(S)
  
  return(S)
}

model_SD <- function(t, State, parms, input)  {
  with(as.list(c(parms, State)), {
    
    import = input(t)
    C = rep(import, length = mcmc_size)
    
    D = State[1:mcmc_size]
    H = State[(mcmc_size+1):(2*mcmc_size)] 
    
    dD <- kd * (C - D)     # internal concentration
    dH <- kk * max(D - z,0) + hb # risk function
    
    res <- c(dD, dH)
    list(res, signal = import)
  })
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

SurvIT_ode <- function(Cw, time, kd, hb, alpha, beta, mcmc_size = NULL, interpolate_length = NULL, interpolate_method = "linear"){
  
  ## external signal with several rectangle impulses
  signal <- data.frame(times = time, 
                       import = Cw)
  
  sigimp <- approxfun(signal$times, signal$import, method = interpolate_method, rule = 2)
  
  if(!is.null(interpolate_length)){
    times <- seq(min(time), max(time), length = interpolate_length)
  } else{
    times <- signal$times
  }
  
  ## The parameters
  parms  <- list( kd = kd,
                  alpha = alpha,
                  beta = beta,
                  mcmc_size = mcmc_size)
  
  ## Start values for steady state
  xstart <- c(D = rep(0, mcmc_size),
              H = rep(0, mcmc_size))
  
  ## Solve model
  out <- ode(y = xstart,
             times = times,
             func = model_IT,
             parms,
             input = sigimp)
  
  mat_4cast <- out %>%
    as.data.frame() %>%
    select(-c(time,signal)) %>%
    as.matrix()
  
  D <- t(mat_4cast)
  D.max <- t(apply(D,1,cummax))
  
  S <- 1-plogis(log(D.max),location=log(parms$alpha),scale=1/parms$beta)
  dtheo <- t(S)
  
  
  return(S)
  
}

model_IT <- function(t, State, parms, input) {
  with(as.list(c(parms, State)), {
    
    import <- input(t)
    C = rep(import, length = mcmc_size)
    
    D = State[1:mcmc_size]
    
    dD <- kd*(C - D)    # internal concentration
    
    list(dD = dD, signal = import)
    
  })
}
