#' Lethal Concentration for \code{x} percent of the population
#' 
#' 
#' 
#' @param object an object of class \code{survFitCstExp}
#' @param X Percentage of population dying: 50 for LC50, 10 for LC10, ...
#' @param time_LCx a number giving the time at which LCx has to be applied. If NULL, lastest time of experiment is used.
#' @param conc_range range of concentrations for 
#' 
#' @export
#' 
LCx.survFitCstExp = function(object, X, time_LCx = NULL, conc_range = NULL, npoints = 100){
  
  df_dose <- doseResponse_survFitCstExp(x = object, time_LCx, conc_range, npoints)
  
  X_prop <- X/100
  
  if(min(df_dose$q50) < X_prop & X_prop < max(df_dose$q50)){
    df.q50 = select(df_dose, c(concentration,q50)) %>%
      add_row(q50 = X_prop) %>%
      arrange(q50) %>%
      mutate(concentration = na.approx(concentration,q50, na.rm=FALSE))%>%
      filter(q50 == X_prop)
    
    LCX_q50 = df.q50$concentration
    
  } else {
    LCX_q50 = NA
    
    warning(paste("No median for LC at", X,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  if(min(df_dose$qinf) < X_prop & X_prop < max(df_dose$qinf)){
    df.qinf=select(df_dose, c(concentration,qinf))%>%
      add_row(qinf=X_prop)%>%
      arrange(qinf)%>%
      mutate(concentration = na.approx(concentration,qinf, na.rm=FALSE))%>%
      filter(qinf==X_prop)
    
    LCX_qinf = df.qinf$concentration
    
  } else{
    LCX_qinf = NA
    
    warning(paste("No 95%inf for LC at", X,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  
  if(min(df_dose$qsup) < X_prop & X_prop < max(df_dose$qsup)){
    df.qsup=select(df_dose, c(concentration,qsup))%>%
      add_row(qsup=X_prop)%>%
      arrange(qsup)%>%
      mutate(concentration = na.approx(concentration,qsup, na.rm=FALSE))%>%
      filter(qsup==X_prop)
    
    LCX_qsup = df.qsup$concentration
    
  } else{
    
    LCX_qsup = NA
    warning(paste("No 95%sup for LC at", X,
                  "% in the range of concentration considered: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  return(c(
    LCX.q50 = LCX_q50,
    LCX.qinf = LCX_qinf,
    LCX.qsup = LCX_qsup
  ))
}

#' dose response curve
#' 
#' 
doseResponse_survFitCstExp <- function(x, time_LCx = NULL,
                                       conc_range = NULL, npoints = 100){
  
  if(is.null(conc_range)){
    conc_range = seq(0, max(x$jags.data$conc), length.out = npoints)
  }
  
  if(is.null(time_LCx)){
    time_LCx = max(x$jags.data$time)
  }
  
  model_type = x$model_type
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  kd <- 10^mctot[, "kd_log10"]
  hb <- 10^mctot[, "hb_log10"]
  
  # all theorical
  k <- 1:length(conc_range)
  j <- 1:npoints
  
  model_type = x$model_type
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
      time_LCx = seq(0, max(x$jags.data$time), 100)
    } else{
      time_LCx = seq()
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
                            qinf = qinf95,
                            qsup = qsup95)
}

Surv_IT_LCx <- function(Cw, time, kd, hb, alpha, beta)
{
  D <- Cw*(1-exp(-kd * time))
  S <- exp(-hb * time) * ( 1- 1/(1 + (D/alpha)^(- beta))) 
  return(S)
}