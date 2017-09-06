#' Lethal Concentration for \code{x} percent of the population
#' 
#' 
#' @param object an object of class \code{survFit}
#' @param X Percentage of population dying: 50 for LC50, 10 for LC10, ...
#' @param time_LCx a number giving the time at which LCx has to be applied. If NULL, latest time of experiment is used.
#' @param conc_range a vector of length 2 with minimal and maximal value of the range of concentration. If NULL, range is
#' define between 0 and the maximal concentration used in the experiment.
#' 
#' @import zoo
#' 
#' @export
#' 
LCx.survFit = function(object, X, time_LCx = NULL, conc_range = NULL, npoints = 100){
  
  if(is.null(conc_range)){
    conc_range = seq(0, max(object$jags.data$conc), length.out = npoints)
  } else{
    if(length(conc_range) != 2){
      stop('conc_range must a vector of length 2 with minimal and maximal value of the range of concentration')
    }
    conc_range = seq(conc_range[1], conc_range[2], length.out = npoints)
  }
  
  if(is.null(time_LCx)){
    time_LCx = max(object$jags.data$time)
  }
  
  df_dose <- doseResponse_survFitCstExp(x = object, time_LCx = time_LCx, conc_range, npoints)
  
  X_prop <- X/100
  
  df_LCx <- pointsLCx(df_dose, X_prop)
  
  object_LCx <- list(X_prop = X/100, time_LCx = time_LCx, df_LCx = df_LCx , df_dose = df_dose)
  class(object_LCx) <- c("LCx", "list")
  
  return(object_LCx)
}

#' dose response curve
#' 
doseResponse_survFitCstExp <- function(x, time_LCx,
                                       conc_range, npoints){
  
  model_type = x$model_type
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  kd <- 10^mctot[, "kd_log10"]
  hb <- 10^mctot[, "hb_log10"]
  
  # all theorical
  k <- 1:length(conc_range)
  j <- 1:npoints
  
  if(model_type == "SD"){
    
    if(is.null(time_LCx)){
      time_LCx = max(x$jags.data$time)
    }
    
    z <- 10^mctot[, "z_log10"]
    kk <- 10^mctot[, "kk_log10"]
    
    dtheo <- lapply(k, function(kit) { # conc
        Surv_SD(Cw = conc_range[kit],
                time = time_LCx,
                kk = kk,
                kd = kd,
                z = z,
                hb = hb)
    })
  }
  if(model_type == "IT"){
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    if(is.null(time_LCx)){
      time_LCx = max(x$jags.data$time)
    }
    
    dtheo <- lapply(k, function(kit) { # concentration pour chaque concentration
      Surv_IT_LCx(Cw = conc_range[kit],
              time = time_LCx,
              kd = kd,
              hb = hb,
              alpha = alpha,
              beta = beta)
    })
  }
  
  # transpose dtheo
  dtheo <- do.call("rbind", lapply(dtheo, t))
  
  # quantile
  qinf95 <- apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE)
  qsup95 <- apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  q50 <- apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE)
  
  df_dose_Resp = data.frame(concentration = conc_range,
                            q50 = q50,
                            qinf95 = qinf95,
                            qsup95 = qsup95)
}

Surv_IT_LCx <- function(Cw, time, kd, hb, alpha, beta)
{
  D <- Cw*(1-exp(-kd * time))
  S <- exp(-hb * time) * ( 1- 1/(1 + (D/alpha)^(- beta))) 
  return(S)
}

#' points for LCx
#' 

pointsLCx <- function(df_dose, X_prop){
  
  if(min(df_dose$q50) < X_prop & X_prop < max(df_dose$q50)){
    df.q50 = select(df_dose, c(concentration, q50)) %>%
      add_row(q50 = X_prop) %>%
      arrange(q50) %>%
      mutate(concentration = zoo::na.approx(concentration,q50, na.rm=FALSE))%>%
      filter(q50 == X_prop)
    
    LCX_q50 = df.q50$concentration
    
  } else {
    LCX_q50 = NA
    
    warning(paste("No median for LC at", X_prop * 100,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  if(min(df_dose$qinf95) < X_prop & X_prop < max(df_dose$qinf95)){
    df.qinf95=select(df_dose, c(concentration,qinf95))%>%
      add_row(qinf95=X_prop)%>%
      arrange(qinf95)%>%
      mutate(concentration = na.approx(concentration,qinf95, na.rm=FALSE))%>%
      filter(qinf95==X_prop)
    
    LCX_qinf95 = df.qinf95$concentration
    
  } else{
    LCX_qinf95 = NA
    
    warning(paste("No 95%inf for LC at", X_prop * 100,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  
  if(min(df_dose$qsup95) < X_prop & X_prop < max(df_dose$qsup95)){
    df.qsup95=select(df_dose, c(concentration,qsup95))%>%
      add_row(qsup95=X_prop)%>%
      arrange(qsup95)%>%
      mutate(concentration = na.approx(concentration,qsup95, na.rm=FALSE))%>%
      filter(qsup95==X_prop)
    
    LCX_qsup95 = df.qsup95$concentration
    
  } else{
    
    LCX_qsup95 = NA
    warning(paste("No 95%sup for LC at", X_prop * 100,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  df_LCx <- data.frame(quantile = c("median", "quantile 2.5%", "quantile 97.5%"),
                       LCx = c(LCX_q50, LCX_qinf95, LCX_qsup95))
  
  return(df_LCx)
  
}


