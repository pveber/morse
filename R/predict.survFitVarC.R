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
#' @keywords transformation
#'
#' @export
#' 
#' @importFrom dplyr bind_rows
#' @importFrom dplyr count
#' 
#' 
predict.survFitVarC <- function(object,
                                newdata = NULL,
                                only.newdata = FALSE,
                                interpolate_length = 100,
                                n.pts.mcmc.thin = NULL){
  
  ##
  ## 1. Creation of the dataset
  ##
  
  x <-  data.frame(time = object$modelData$time_long,
                   conc = object$modelData$conc_long,
                   replicate = object$modelData$replicate_ID_long)
  
  if(!is.null(newdata)){
    if(only.newdata == TRUE){
      x <- newdata
    } else{
      x <- rbind(x, newdata)
    }
  } 
  
  mat_survData = x[x != 0, ] %>%
    arrange(replicate,time) # ARRANGE To make proper data matrix
  
  count_data = count(mat_survData, replicate)
  max.count_data = max(count_data$n)
  
  ls_out = split( mat_survData , f = mat_survData$replicat )
  
  mat_time = ls_out%>%
    lapply("[[", "time")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)

  mat_conc = ls_out%>%
    lapply("[[", "conc")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  mat_replicate = ls_out%>%
    lapply("[[", "replicate")%>%
    lapply(`length<-`, max.count_data)%>%
    unlist()%>%
    matrix(ncol = max.count_data, byrow = TRUE)
  
  ##
  ## 2. Simulations
  ##

  mcmc.samples = object$mcmc
  
  # ---------REDUCTION OF MCMC SIZE FOR PLOTTING EFFICIENCY
  if(!is.null(n.pts.mcmc.thin)){
    reduc_tab = lapply(mcmc.samples,"[",
                       #(nrow(mcmc.samples[[1]])-5e4):nrow(mcmc.samples[[1]]),
                       seq(1,nrow(mcmc.samples[[1]]),nrow(mcmc.samples[[1]])/n.pts.mcmc.thin),
                       1:4)
    
    mcmc.samples = reduc_tab
  }
  # -------------------------------------------------
  
  mctot <- do.call("rbind", mcmc.samples)
  kd <- 10^mctot[, "kd_log10"]
  hb <- 10^mctot[, "hb_log10"]
  
  k = 1:nrow(mat_time)
  
  if(object$model_type == "SD"){
    kk <- 10^mctot[, "kk_log10"]
    z <- 10^mctot[, "z_log10"]
    
    dtheo_ls  <- lapply(k, function(kit) { # For every group
      Surv.SD_Cext (Cw = mat_conc[kit,],
                    time = mat_time[kit,],
                    kd = kd,
                    hb = hb,
                    kk = kk,
                    z = z)
    })
  }
  if(object$model_type == "IT"){
    
    alpha <- 10^mctot[, "alpha_log10"]
    beta <- 10^mctot[, "beta_log10"]
    
    dtheo_ls  <- lapply(k, function(kit) { # Pour chaque replicat
      Surv.IT_Cext (Cw = mat_conc[kit,],
                    time = mat_time[kit,],
                    kd = kd,
                    hb = hb,
                    alpha = alpha,
                    beta = beta)
    })
  }
  # transpose dtheo
  dtheo <- do.call("rbind", lapply(dtheo_ls, t))
  
  ##
  ## 3. Building of return object
  ##
  
  # -------- df.theo
  
  predict_df = data.frame(
    conc = as.vector(t(mat_conc)),
    time = as.vector(t(mat_time)),
    replicate = as.vector(t(mat_replicate)),
    #
    q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
    qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE),
    qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE))

  predict_df = filter(predict_df, !is.na(replicate)) # resulting from the matrix building
  
  #---- Observation
  observ_df <-  data.frame(time = object$modelData$time,
                           Nsurv = object$modelData$Nsurv,
                           replicate_name = object$modelData$replicate,
                           replicate = object$modelData$replicate_ID)
  
  
  OUT_ls <- list(predict_data = predict_df,
                 observ_data = observ_df)
  
  
  class(OUT_ls) <- c("predictFitCst","predictFit")
  
  return(OUT_ls)
  
}


#' Survival function for "IT" model with external concentration changing with time
#'
#' @param Cw A scalar of external concentration
#' @param time A vector of time
#' @param kk a vector of parameter
#' @param kd a vector of parameter
#' @param z a vector of parameter
#' @param hb a vector of parameter
#' 
#'
#' @return A matrix generate with coda.samples() function
#'
#' @export
#'
#'

Surv.SD_Cext <- function(Cw, time, kk, kd, z, hb){
  
  time.prec = dplyr::lag(time, 1) ; time.prec[1] = time[1] #   time[1] = tprec[1]
  diff.int = (exp(time %*% t(kd)) * Cw + exp(time.prec %*% t(kd)) * Cw )/2 * (time-time.prec) #OK time[1]-tprec[1] = 0
  
  D = kd * exp(-kd %*% t(time)) * t(apply(diff.int,2,cumsum))
  
  lambda = kk * pmax(D-z,0) + hb # the pmax function is important here for elementwise maximum with 0 and D[i,j]-z ATTENTION: pmax(0,D) != pmax(D,0)
  
  lambda.prec = dplyr::lag(lambda, 1 ) ; lambda.prec[1] = lambda[1]
  
  int.lambda =  t(t((lambda + lambda.prec)/2) * (time-time.prec))
  
  psurv <- exp(-t(apply(int.lambda,1,cumsum)))
  
  return(psurv)
}

#' Survival function for "IT" model with external concentration changing with time
#'
#' @param Cw A scalar of external concentration
#' @param time A vector of time
#' @param kk a vector of parameter
#' @param kd a vector of parameter
#' @param z a vector of parameter
#' @param hb a vector of parameter
#' 
#'
#' @return A matrix generate with coda.samples() function
#'
#' @export

Surv.IT_Cext <- function(Cw, time, kd, hb, alpha, beta)
{
  time.prec = dplyr::lag(time, 1) ; time.prec[1] = time[1] #   time[1] = tprec[1]
  diff.int = (exp(time %*% t(kd)) * Cw + exp(time.prec %*% t(kd)) * Cw )/2 * (time-time.prec) #OK time[1]-tprec[1] = 0
  
  D = kd * exp(-kd %*% t(time)) * t(apply(diff.int,2,cumsum))
  
  D.max=t(apply(D,1,cummax))
  S <-  exp(-hb %*% t(time))*(1-plogis(log(D.max),location=log(alpha),scale=1/beta))
  return(S)
}
