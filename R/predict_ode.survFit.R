#' Predict method for \code{survFit} objects
#' 
#' This is a \code{method} to replace function \code{predict} used on \code{survFit}
#' object when computing issues happen. \code{predict_ode} uses the \code{deSolve}
#' library to improve robustness. However, time to compute may be longer.
#' 
#' 
#' @param object an object used to select a method \code{ppc}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
predict_ode <- function(object, ...){
  UseMethod("predict_ode")
}

#' Predict method for \code{survFit} objects
#'
#' This is the generic \code{predict} S3 method for the \code{survFit} class.
#' It provides predicted survival rate for "SD" or "IT" models under constant or time-variable exposure.
#'
#' @param object An object of class \code{survFit}.
#' @param data_predict A dataframe with three columns \code{time}, \code{conc} and \code{replicate}
#'  used for prediction. If \code{NULL}, prediction is based on \code{x} object of 
#'  class \code{survFit} used for fitting.
#' @param spaghetti If \code{TRUE}, return a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param mcmc_size Can be used to reduce the number of mcmc samples in order to speed up
#'  the computation. \code{mcmc_size} is the number of selected iterations for one chain. Default
#'  is 1000. If all MCMC is wanted, set argument to \code{NULL}.
#' @param hb_value If \code{TRUE}, the background mortality \code{hb} is taken into account from the posterior.
#' If \code{FALSE}, parameter \code{hb} is set to a fixed value. The default is \code{TRUE}.
#' @param interpolate_length Length of the time sequence for which output is wanted.
#' @param interpolate_method The interpolation method for concentration. See package \code{deSolve} for details.
#' Default is \code{linear}.
#' @param  hb_valueFORCED If \code{hb_value} is \code{FALSE}, it fix \code{hb}.
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
#' # (5) Predict on a new data set
#' predict_out <- predict_ode(object = out, data_predict = data_4prediction,
#'                            mcmc_size = 1000, spaghetti = TRUE)
#'
#' }
#' 
#' @import deSolve
#' @importFrom stats approxfun
#' 
#' @export
#'
predict_ode.survFit <- function(object,
                                data_predict = NULL,
                                spaghetti = FALSE,
                                mcmc_size = 1000,
                                hb_value = TRUE,
                                interpolate_length = 100,
                                interpolate_method = "linear",
                                hb_valueFORCED = NA,
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
  #if(is.null(mcmc_size)){
  mcmc_size = nrow(mctot)
  #}
  
  kd = 10^mctot[, "kd_log10"]
  
  if(hb_value == TRUE){
    hb <- 10^mctot[, "hb_log10"]
  } else if(hb_value == FALSE){
    if(is.na(hb_valueFORCED)){
      if(is.na(x$hb_valueFIXED)){
        stop("Please provide value for `hb` using `hb_valueFORCED`.")
      } else{
        hb <- rep(x$hb_valueFIXED, nrow(mctot))
      } 
    } else{
      hb <- rep(hb_valueFORCED, nrow(mctot))
    }
  }
  
  k = 1:length(unique_replicate)
  
  if(model_type == "SD"){
    kk <- 10^mctot[, "kk_log10"]
    z <- 10^mctot[, "z_log10"]
    
    dtheo = lapply(k, function(kit) { # For each replicate
      SurvSD_ode(Cw = ls_conc[[kit]],
                 time = ls_time[[kit]],
                 replicate = unique_replicate[kit],
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
                 replicate = unique_replicate[kit],
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
  df_theo <- do.call("rbind", dtheo)
  
  df_quantile = select(df_theo, time, conc, replicate, q50, qinf95, qsup95)
  
  if(spaghetti == TRUE){
    df_spaghetti <- df_theo
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

SurvSD_ode <- function(Cw, time, replicate, kk, kd, z, hb, mcmc_size = 1000, interpolate_length = NULL, interpolate_method=c("linear","constant")) {
  interpolate_method <- match.arg(interpolate_method)
  
  ## external signal with several rectangle impulses
  signal <- data.frame(times=time,import=Cw)
  if(!is.null(interpolate_length)){
    times <- seq(min(time), max(time), length = interpolate_length)
  } else{
    times <- signal$times
  }
  
  xstart <- c(rep(c(D=0),mcmc_size),rep(c(H=0),mcmc_size))
  # ordering of parameters required by compiled function
  parms <- c(mcmc_size,kd,hb,z,kk)
  # solve model
  on.exit(.C("gutsredsd_free")) # clean up
  deSolve::ode(y=xstart,
               times=times,
               parms=parms,
               method="lsoda",
               dllname="morse",
               initfunc="gutsredsd_init",
               func="gutsredsd_func",
               initforc="gutsredsd_forc",
               forcings=signal,
               fcontrol=list(method=interpolate_method,rule=2,ties="ordered"),
               nout=1
  ) -> out
  
  dtheo <- exp(-out[,grep("H",colnames(out))])

  # Manage vector case
  if(mcmc_size == 1){
    q50 = dtheo
    qinf95 = dtheo
    qsup95 = dtheo
  } else{
    qs <- apply(as.matrix(dtheo), 1, quantile, probs=c(0.5,0.025,0.975), names=FALSE, na.rm=TRUE)
    q50 = qs[1,]
    qinf95 = qs[2,]
    qsup95 = qs[3,]
  }
  
  dtheo <- as.data.frame(dtheo)
  names(dtheo) <- paste0("H",seq(1,mcmc_size))
  dtheo <- dtheo %>%
    dplyr::mutate(time = times,
           conc = out[,ncol(out)],
           replicate = c(replicate),
           q50 = q50,
           qinf95 = qinf95,
           qsup95 = qsup95)
  
  return(dtheo)
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

SurvIT_ode <- function(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = NULL, interpolate_length = NULL, interpolate_method=c("linear","constant")){
  interpolate_method <- match.arg(interpolate_method)
  
  ## external signal with several rectangle impulses
  signal <- data.frame(times=time,import=Cw)
  if(!is.null(interpolate_length)){
    times <- seq(min(time), max(time), length = interpolate_length)
  } else{
    times <- signal$times
  }
  
  ## The parameters
  parms  <- c(mcmc_size,kd,hb)
  
  ## Start values for steady state
  xstart <- c(rep(c(D=0),mcmc_size),rep(c(H=0),mcmc_size))
  
  ## Solve model
  on.exit(.C("gutsredit_free")) # clean up
  deSolve::ode(y=xstart,
               times=times,
               parms=parms,
               method="lsoda",
               dllname="morse",
               initfunc="gutsredit_init",
               func="gutsredit_func",
               initforc="gutsredit_forc",
               forcings=signal,
               fcontrol=list(method=interpolate_method,rule=2,ties="ordered"),
               nout=1
  ) -> out

  D <- out[,grep("D",colnames(out))]
  cumMax_D <- if(is.null(dim(D))) cummax(D) else apply(D, 2, cummax)
  thresholdIT <- t(1 / (1 + (t(cumMax_D) / alpha)^(-beta)))
  
  dtheo <- (1 - thresholdIT) * exp(times %*% t(-hb))

  # Manage vector case
  if(mcmc_size == 1){
    q50 = dtheo
    qinf95 = dtheo
    qsup95 = dtheo
  } else{
    qs <- apply(as.matrix(dtheo), 1, quantile, probs=c(0.5,0.025,0.975), names=FALSE, na.rm=TRUE)
    q50 = qs[1,]
    qinf95 = qs[2,]
    qsup95 = qs[3,]
  }
  
  dtheo <- as.data.frame(dtheo)
  names(dtheo) <- paste0("H",seq(1,mcmc_size))
  dtheo <- dtheo %>%
    dplyr::mutate(time = out[, "time"],
            conc = out[,ncol(out)],
            replicate = c(replicate),
            q50 = q50,
            qinf95 = qinf95,
            qsup95 = qsup95)

  return(dtheo)
  
}
