#' Predict \eqn{x}\% Lethal Concentration at any specified time point for 
#' a \code{survFit} object.
#' 
#' The function \code{LCx}, \eqn{x}\% Lethal Concentration (\eqn{LC_x}), is use to compute
#'  the dose required to kill \eqn{x}\% of the members of a tested population
#'  after a specified test duration (\code{time_LCx}) (default is the maximum
#'  time point of the experiment).
#'  
#'  Mathematical definition of \eqn{x}\% Lethal Concentration at time \eqn{t},
#'  denoted \eqn{LC(x,t)}, is:
#'  
#'  \eqn{S(LC(x,t), t) = S(0, t)*(1- x/100)},
#'  
#'  where \eqn{S(LC(x,t), t)} is the survival probability at concentration
#'  \eqn{LC(x,t)} at time \eqn{t}, and \eqn{S(0,t)} is the survival probability at
#'  no concentration (i.e. concentration is \eqn{0}) at time \eqn{t} which
#'  reflect the background mortality \eqn{h_b}:
#'  
#'  \eqn{S(0, t) = exp(-hb* t)}.
#'   
#'  In the function \code{LCx}, we use the median of \eqn{S(0,t)} to rescale the
#'  \eqn{x}\% Lethal Concentration at time \eqn{t}.
#'  
#' @rdname LCX
#' 
#' @param object An object of class \code{survFit}
#' @param X Percentage of individuals dying (e.g., \eqn{50} for \eqn{LC_{50}}, \eqn{10} for \eqn{LC_{10}}, ...)
#' @param time_LCx A number giving the time at which  \eqn{LC_{x}} has to be estimated. 
#' If NULL, the latest time point of the experiment is used.
#' @param conc_range A vector of length 2 with minimal and maximal value of the 
#' range of concentration. If NULL, the range is
#' define between 0 and the highest tested concentration of the experiment.
#' @param npoints Number of time point in \code{conc_range} between 0 and the maximal concentration. 100 by default.
#' @param \dots Further arguments to be passed to generic methods
#'
#' @return The function returns an object of class \code{LCx}, which is a list
#'  with the following information:
#' \item{X_prop}{Survival probability of individuals surviving considering the median
#'  of the background mortality (i.e. \eqn{S(0, t)*(1- x/100)})}
#' \item{X_prop_provided}{Survival probability of individuals surviving as provided in arguments (i.e. \eqn{(100-X)/100)}}
#' \item{time_LCx}{A number giving the time at which  \eqn{LC_{x}} has to be
#'  estimated as provided in arguments or if NULL, the latest time point of the
#'   experiment is used.}
#' \item{df_LCx}{A \code{data.frame} with quantiles (median, 2.5\% and 97.5\%)
#'  of \eqn{LC_{X}} at time \code{time_LCx} for \eqn{X}\% of individuals}
#' \item{df_dose}{A \code{data.frame} with four columns: \code{concentration}, and median \code{q50} and 95\% credible interval
#'  (\code{qinf95} and \code{qsup95}) of the survival probability at time \code{time_LCx}}
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
#' # (4) estimate LC50 at time 4
#' LCx(out_SD, X = 50, time_LCx = 4)
#' }
#' 
#' @import zoo
#' @importFrom stats approx
#' 
#' @export
#' 
LCx.survFit <- function(object,
                        X,
                        time_LCx = NULL,
                        conc_range = NULL,
                        npoints = 100,
                        ...){
  
  if(is.null(conc_range)){
    conc_range = seq(0, max(object$jags.data$conc), length.out = npoints)
  } else{
    if(length(conc_range) != 2){
      stop('conc_range must a vector of length 2 with minimal and maximal value of the range of concentration')
    }
    conc_range = seq(conc_range[1], conc_range[2], length.out = npoints)
  }
  
  if(min(conc_range) != 0){
    stop("Minimal value of 'conc_range' must be 0.")
  }
  
  if(is.null(time_LCx)){
    time_LCx = max(object$jags.data$time)
  }
  
  df_dose <- doseResponse_survFitCstExp(x = object, time_LCx = time_LCx, conc_range, npoints)
  
  median_backgroundMortality_Conc0 = dplyr::filter(df_dose, concentration == 0)$q50
  
  X_prop_provided <- (100-X)/100
  
  X_prop <- (100-X)/100*median_backgroundMortality_Conc0
  
  df_LCx <- pointsLCx(df_dose, X_prop)
  
  object_LCx <- list(X_prop = X_prop,
                     X_prop_provided = X_prop_provided,
                     time_LCx = time_LCx,
                     df_LCx = df_LCx,
                     df_dose = df_dose)
  class(object_LCx) <- c("LCx", "list")
  
  return(object_LCx)
}

# dose response curve
# 
doseResponse_survFitCstExp <- function(x, time_LCx,
                                       conc_range, npoints){
  
  model_type = x$model_type
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  kd <- 10^mctot[, "kd_log10"]
  # "hb" is not in survFit object of morse <v3.2.0
  if("hb" %in% colnames(mctot)){
    hb <- mctot[, "hb"]  
  } else{ hb <- 10^mctot[, "hb_log10"] }
  
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

# points for LCx
# 
pointsLCx <- function(df_dose, X_prop){
  
  if(min(df_dose$q50) < X_prop & X_prop < max(df_dose$q50)){
    LCX_q50 = approx( df_dose$q50, df_dose$concentration, xout = X_prop, ties = mean)$y
  } else {
    LCX_q50 = NA
    warning(paste("No median for survival probability of", X_prop,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  if(min(df_dose$qinf95) < X_prop & X_prop < max(df_dose$qinf95)){
    LCX_qinf95 = approx( df_dose$qinf95, df_dose$concentration, xout = X_prop, ties = mean)$y
  } else{
    LCX_qinf95 = NA
    warning(paste("No 95%inf for survival probability of", X_prop ,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  if(min(df_dose$qsup95) < X_prop & X_prop < max(df_dose$qsup95)){
    LCX_qsup95 = approx( df_dose$qsup95, df_dose$concentration, xout = X_prop, ties = mean)$y
  } else{
    LCX_qsup95 = NA
    warning(paste("No 95%sup for survival probability of", X_prop,
                  " in the range of concentrations under consideration: [",
                  min(df_dose$concentration), ";", max(df_dose$concentration), "]"))
  }
  
  df_LCx <- data.frame(quantile = c("median", "quantile 2.5%", "quantile 97.5%"),
                       LCx = as.numeric(c(LCX_q50, LCX_qinf95, LCX_qsup95)))
    # as.numeric is needed here because if all values are NA, LCx has type logical
  return(df_LCx)
}


