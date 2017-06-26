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
#' @importFrom deSolve ode
#' @importFrom dplyr bind_rows
#' 
#' 
#' 
predict.survFitCstC <- function(object,
                                newdata = NULL,
                                only.newdata = FALSE,
                                interpolate_length = 100){
  
  
    ##
    ## 1. Creation of the dataset
    ##

    x <-  data.frame(time = object$modelData$time,
                     conc = object$modelData$conc,
                     replicate = object$modelData$replicate)
    
    if(!is.null(newdata)){
      if(only.newdata == TRUE){
        x <- newdata
      } else{
        x <- rbind(x, newdata)
      }
    } 
    
    conc.obs = unique(x$conc)
    replicate.obs = unique(x$replicate)
    
    tfin <- seq(0, max(x$time), length.out = interpolate_length)
    
    ##
    ## 2. Simulations
    ##
    
    mcmc.samples = object$mcmc
    
    mctot <- do.call("rbind", mcmc.samples)
    kd <- 10^mctot[, "kd_log10"]
    hb <- 10^mctot[, "hb_log10"]
   
    # all theorical
    k <- 1:length(conc.obs)
    j <- 1:interpolate_length
    
    if(object$model_type == "SD"){
      # parameters
      kk <- 10^mctot[, "kk_log10"]
      z <- 10^mctot[, "z_log10"]
      
      dtheo_ls <- lapply(k, function(kit) { # concentration pour chaque concentration
        sapply(j, function(n.t) { # time
          Surv.SD(Cw = conc.obs[kit],
                  time = tfin[n.t],
                  kk = kk,
                  kd = kd,
                  z = z,
                  hb = hb)
        })
      })
    } 
    if(object$model_type=="IT"){
      # parameters
      alpha <- 10^mctot[, "alpha_log10"]
      beta <- 10^mctot[, "beta_log10"]
      
      dtheo_ls <- lapply(k, function(kit) { # concentration pour chaque concentration
        Surv.IT(Cw = conc.obs[kit],
                time = tfin,
                kd = kd,
                hb = hb,
                alpha = alpha,
                beta = beta)
      })
    }
    # transpose dtheo
    dtheo = do.call("rbind", lapply(dtheo_ls, t))  

    ##
    ## 3. Building of return object
    ##
    
    predict_df = data.frame(conc = rep(conc.obs, rep(interpolate_length, length(conc.obs))),
                            replicate = rep(replicate.obs, rep(interpolate_length, length(conc.obs))),
                            time = rep(tfin, length(conc.obs)),
                            q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
                            qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE) ,
                            qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE))  
    
    
    #---- Observation
      observ_df <-  data.frame(time = object$modelData$time,
                               conc = object$modelData$conc,
                               Nsurv = object$modelData$Nsurv,
                               replicate = object$modelData$replicate)
    
    
    OUT_ls <- list(predict_data = predict_df,
                   observ_data = observ_df)
    
    
    class(OUT_ls) <- c("predictFitCst","predictFit")
    
    return(OUT_ls)
    
  }


#' Survival function for "SD" model
#'
#' @param Cw A scalar of external concentration
#' @param time A scalar of time
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

Surv.SD <- function(Cw, time, kk, kd, z, hb)
{
  S <- exp(-hb*time)
  
  x <- ifelse(Cw > z, 1 - z/Cw, NA)
  tNEC <- -(1/kd)*log(x)
  
  y <- ifelse(time > tNEC,
              exp( kk/kd*Cw*(exp(-kd*tNEC) -exp(-kd*time))
                   - kk*(Cw-z)*(time - tNEC)),
              NA)
  
  return(ifelse(!is.na(x) & !is.na(y), S * y, S))
}

#' Survival function for "IT" model
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

Surv.IT <- function(Cw, time, kd, hb, alpha, beta)
{
  D = Cw*(1-exp(-kd %*% t(time)))
  D.max = t(apply(D,1,cummax))
  S = exp(-hb %*% t(time))*(1-plogis(log(D.max),location=log(alpha),scale=1/beta))
  return(S)
}