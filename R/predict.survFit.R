#' \code{predict} is a generic function for predictions from \code{survFiT} model. 
#' The function invokes particular methods which depend on the class of the first argument.
#'
#' This function creates a \code{predictFit} object from a \code{survFit}
#' object.
#'
#' @param object An object of class \code{survFit}
#' @param newdata An optional data frame in which to look for variables with which to predict.
#' If omitted, the fitted values are used.
#'
#' @return An object of class \code{reproData}.
#'
#' @keywords transformation
#'
#' @export
#' 
#' @importFrom deSolve ode
#' @importFrom dplyr bind_rows
#' 
#' 
#' 
predict.survFit <- function(object,
                            newdata = NULL,
                            only.newdata = FALSE,
                            n_size = 100, # number of iteration selected
                            interpolate_length = 100,
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
  
  ### 1. select object and newdata
  if(!is.null(object$modelData$conc)){
    x <- data.frame(time = object$modelData$time,
                    conc = object$modelData$conc,
                    replicate = object$modelData$replicate)
  } else{
    x <- data.frame(time = object$modelData$time_long,
                    conc = object$modelData$conc_long,
                    replicate = object$modelData$replicate_ID_long)
  }
  
  if(!is.null(newdata)){
    if(only.newdata == TRUE){
      x <- newdata
    } else{
      x <- rbind(x, newdata)
    }
  } 
  
  ### 2. Check number of replicate
  
  replicate_unique <- unique(x$replicate)
  replicate_uniqueLength <- length(replicate_unique)
  
  ### 3. run the ode
  
  ls_return <- lapply(1:replicate_uniqueLength,
                      function(replic) ode_predict(x, parms, replicate_unique, replic, n_size,
                                                   interpolate_method, interpolate_length, model_type = object$model_type))
  
  ### 3. produce the output data.frame
  
  #---- Predict
  
  predict_df <- dplyr::bind_rows(ls_return) %>%
    as_data_frame()
  
  predict_df$replicate <- rep(replicate_unique, rep(interpolate_length, replicate_uniqueLength))
  
  #---- Observation
  
  if(!is.null(object$modelData$conc)){
    observ_df <-  data.frame(time = object$modelData$time,
                              conc = object$modelData$conc,
                              Nsurv = object$modelData$Nsurv,
                              replicate = object$modelData$replicate)
  } else{
    observ_df <-  data.frame(time = object$modelData$time,
                             Nsurv = object$modelData$Nsurv,
                             replicate_name = object$modelData$replicate,
                             replicate = object$modelData$replicate_ID)
  }
  
    
    OUT_ls <- list(predict_data = predict_df,
                   observ_data = observ_df)
  

  class(OUT_ls) <- "predictFit"

  return(OUT_ls)
  
}


#' ode_predict
#' 
ode_predict <- function(x, parms, replicate_unique, replic, n_size, interpolate_method, interpolate_length, model_type){
  
  filter_x <- filter(x, replicate == replicate_unique[replic])
  
  if(is.null(filter_x$conc)){
    signal <- data.frame(times = filter_x$time_long, 
                         import = filter_x$conc_long)
  } else{
    signal <- data.frame(times = filter_x$time, 
                         import = filter_x$conc)
  }
  
  sigimp <- approxfun(signal$times, signal$import, method = interpolate_method, rule = 2)
  
  ### time vector if interpolate_length is specified
  if(!is.null(interpolate_length)){
    ode_times = seq(min(filter_x$time),max(filter_x$time),length = interpolate_length)
  } else{
    ode_times = signal$times
  }
  
  ### parameterization of initial condition + ODE model to use SD or IT
  
  if(model_type == "IT"){
    ## Start values for steady state
    xstart <- c(D = rep(0,n_size))
    ## model
    model_TKTD <- model_TKTD_IT
  }
  if (model_type == "SD"){
    xstart <- c(D = rep(0,n_size),
                H = rep(0,n_size))
    ## model
    model_TKTD <- model_TKTD_SD
  }
  
  ### Solve  ODE system
  out <- ode(y = xstart,
             times = ode_times,
             func = model_TKTD,
             parms,
             input = sigimp)
  
  out_df = as.data.frame(out)
  
  ##
  ## 3. Construction of the object to return 
  ##
  
  if(model_type == "IT"){
    mat_predict <- out_df %>%
      select(-c(time,signal)) %>%
      as.matrix()
    D <- t(mat_predict)
    D.max <- t(apply(D, 1, cummax))
    
    S <- exp(-parms$hb %*% t(out_df$time))*(1 - plogis(log(D.max), location = log(parms$alpha), scale = 1 / parms$beta))
    dtheo <- t(S)
    
  } 
  if(model_type == "SD"){
    
    mat_predict_H <- out_df %>%
      select(contains("H"),-c(time,signal)) %>%
      as.matrix()
    S <- exp(-t(mat_predict_H))
    dtheo <- t(S)
    
  }
  
  OUT_ode <- data_frame( time = out_df$time,
                        conc = out_df$signal,
                        q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
                        qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE),
                        qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE))
  
  
  return(OUT_ode)
  
}


# ----------------------------

model_TKTD_SD <- function(t, State, parms, input)  {
  with(as.list(c(parms, State)), {
    
    import <- input(t)
    C = rep(import, length = n_size)
    
    D = State[1:n_size]
    H = State[(n_size+1):(2*n_size)] 
    
    dD <- kd*(C - D)        # internal concentration
    dH <- kk*max(D-z,0) + hb    # risk function
    
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

