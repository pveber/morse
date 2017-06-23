#' Return prediction from the results of model fitting with \code{survFit}
#' function.
#'
#' This function creates a \code{predictFit} object from a \code{survFit}
#' object.
#'
#' @aliases reproData
#'
#' @param object a dataframe of class \code{survFit}
#'
#' @return An object of class \code{reproData}.
#'
#' @keywords transformation
#'
#' @export
#' 
#' @importFrom deSolve ode
#' 
predict.survFit <- function(object,
                            time_predict,
                            conc_predict,
                            n_size,
                            interpolate_length = NULL,
                            interpolate_method = "linear"){
  
  ##
  ## 1. Recover parameters
  ##
  mctot = do.call("rbind", object$mcmc)
  mctot_samples = mctot[round(seq(1,nrow(mctot),length = n_size)),]
  
  parms  <- list( n_size = n_size,
                  kd = 10^mctot_samples[, "kd_log10"],
                  hb = 10^mctot_samples[, "hb_log10"])
  
  if(object$model_type == "IT"){
    parms$alpha <- 10^mctot_samples[, "alpha_log10"]
    parms$beta <- 10^mctot_samples[, "beta_log10"]
  }
  if (object$model_type == "SD"){
    parms$z <- 10^mctot_samples[, "z_log10"]
    parms$kk <- 10^mctot_samples[, "kk_log10"]
  } 
  
  
  ##
  ## 2. ODE solver 
  ##
  
  ### external signal with several rectangle impulses
  signal <- data.frame(times = time_predict, 
                       import = conc_predict)
  
  sigimp <- approxfun(signal$times, signal$import, method = interpolate_method, rule = 2)
  
  ### time vector if interpolate_length is specified
  if(!is.null(interpolate_length)){
    times = seq(min(time_predict),max(time_predict),length = interpolate_length)
  } else{
    times = signal$times
  }
  
  ### parameterization of initial condition + ODE model to use SD or IT

  if(object$model_type == "IT"){
    ## Start values for steady state
    xstart <- c(D = rep(0,n_size))
    ## model
    model_TKTD <- model_TKTD_IT
  }
  if (object$model_type == "SD"){
    xstart <- c(D = rep(0,n_size),
                H = rep(0,n_size))
    ## model
    model_TKTD <- model_TKTD_SD
  } 
  
  ### Solve  ODE system
  out <- ode(y = xstart,
             times = times,
             func = model_TKTD,
             parms,
             input = sigimp)
  
  ##
  ## 3. Construction of the object to return 
  ##
  
  if(object$model_type == "IT"){
    mat_predict <- out%>%
      as.data.frame()%>%
      select(-c(time,signal))%>%
      as.matrix()
    D <- t(mat_predict)
    D.max <- t(apply(D, 1, cummax))
    
    S <- 1 - plogis(log(D.max), location = log(parms$alpha), scale = 1 / parms$beta)
    dtheo <- t(S)
    
  } 
  if(object$model_type == "SD"){
    
    mat_predict_H <- out%>%
      as.data.frame()%>%
      select(contains("H"),-c(time,signal))%>%
      as.matrix()
    S <- exp(-t(mat_predict_H))
    dtheo <- t(S)
    
  }
  
  # -------- df.theo
  
  OUT_df <- data_frame( time = signal$times,
                        conc = signal$import,
                        q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
                        qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE),
                        qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE))
  
  
  OUT_ls <- list(predict = OUT_df,
                n_size = parms$n_size)
  
  class(OUT_ls) <- c("predictFit")
  
  return(OUT_ls)
  
}


#' Forecasting ith ODE solver


# ----------------------------

model_TKTD_SD <- function(t, State, parms, input)  {
  with(as.list(c(parms, State)), {
    
    import <- input(t)
    C = rep(import, length = n_size)
    
    D = State[1:n_size]
    H = State[(n_size+1):(2*n_size)] 
    
    dD <- kd*(C - D)        # internal concentration
    dH <- kk*max(D-z,0)     # risk function
    
    res <- c(dD, dH)
    list(res, signal = import)
    
  })
}

model_TKTD_IT <- function(t, State, parms, input) {
  with(as.list(c(parms, State)), {
    
    import <- input(t)
    C = rep(import, length = n_size)
    
    D = State[1:n_size]
    
    dD <- kd*(C - D)    # internal concentration
    
    list(dD = dD, signal = import)
    
  })
}
