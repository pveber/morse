#' Plotting method for \code{survFit} objects
#'
#' This is the generic \code{plot} S3 method for the
#' \code{survFit}.  It plots the fits obtained for each
#' concentration of pollutant in the original dataset.
#'
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration.
#' The black dots depict the \strong{observed survival
#' rate} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95 \% credible intervals for the estimated survival
#' rate (by default the grey area around the fitted curve) and 95 \% confidence
#' intervals for the observed survival rate (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (2 \% of the MCMC chains are randomly
#' taken for this sample).
#'
#' @param x An object of class \code{survFit}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param spaghetti if \code{TRUE}, draws a set of survival curves using
#' parameters drawn from the posterior distribution.
#' @param one.plot if \code{TRUE}, draws all the estimated curves in
#' one plot instead of one per concentration.
#' @param adddata if \code{TRUE}, adds the observed data to the plot
#' with (frequentist) confidence intervals.
#' @param mcmc_size A scalar refering to the size of the MCMC used for
#'  sampling parameters, by default, all the MCMC is used.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @keywords plot
#'
#' @export
#' 
#'
#'
plot.survFitVarExp <- function(x,
                               xlab = "Time",
                               ylab = "Survival rate",
                               main = NULL,
                               spaghetti = FALSE,
                               one.plot = FALSE,
                               adddata = TRUE,
                               mcmc_size = NULL,
                               ...) {
  
  
  df_predictTotal <- predict.survFitVarExp(x, spaghetti, mcmc_size)
  
  df_prediction <-  df_predictTotal$df_quantile
  
  df_observation <- filter(x$original.data, !is.na(Nsurv))
  
  # Plot
  
  plt <- ggplot() +
      theme_minimal() +
      theme(legend.position ="top",
            strip.background = element_rect(fill="grey90", color = "white"),
            strip.text = element_text(colour = "black")) +
      scale_x_continuous(name = xlab) +
      scale_y_continuous(name = ylab,
                         limits = c(0,1)) +
      theme(legend.position = "top") +
    # Prediction
    geom_ribbon(data = df_prediction,
                aes(x = time, ymin = qinf95,ymax = qsup95, group = replicate),
                fill = "grey30", alpha = 0.4) +
    geom_line(data = df_prediction,
              aes(x = time, y = q50, group = replicate),
              col="orange", size = 1)
  
  # Observation
  if(adddata == TRUE){
    plt <- plt +
      geom_point(data = df_observation,
               aes(x = time, y = Nsurv/Ninit, group = replicate))
  }
  
  # # spaghetti
  # if(spaghetti == TRUE){
  #   
  #   df_spaghetti <- predictTotal$df_spaghetti[, 1:1000]
  #     
  #   plt <- plt +
  #     geom_lines(data = df_observation,
  #                aes(x = time, y = ??? , group = replicate))
  # }
    
  # facetting
  if(one.plot == FALSE){
    plt <- plt + facet_wrap(~ replicate)
  }  
      
   return(plt)
  }

#' Theoretical simulation for "SD" or "IT" model with time varying exposure
#'
#' @param x A object of class \code{survFitVarExp}
#' @param mcmc_size \code{mcmc_size} is by default the size of the mcmc in object \code{x}.
#'  \code{mcmc_size} can be used to reduce the number of mcmc samples in order to speed up
#'  the computation.
#'  
#'    
#' @return A data.frame with theoretical simulations
#'
#' @export
#'
predict.survFitVarExp <- function(x,
                                  spaghetti = FALSE,
                                  mcmc_size = NULL){
  
  # Initialisation
  
  mcmc <- x$mcmc
  
  model_type <- x$model_type
  
  df <- data.frame(
    time = x$jags.data$time_long,
    conc = x$jags.data$conc_long,
    replicate = x$jags.data$replicate_long
  )
  
  unique_replicate <- unique(df$replicate)
  
  ls_time <- list()
  ls_conc <- list()
  
  for(i in 1:length(unique_replicate)){
    
    ls_time[[i]] <- filter(df, replicate == unique_replicate[i])$time
    ls_conc[[i]] <- filter(df, replicate == unique_replicate[i])$conc
    
  }
  
  # ------- Computing
  
  mcmc.samples = mcmc
  
  if(!is.null(mcmc_size)){
    reduc_tab = lapply(mcmc.samples,"[",
                       #(nrow(mcmc.samples[[1]])-5e4):nrow(mcmc.samples[[1]]),
                       seq(1,nrow(mcmc.samples[[1]]), length = mcmc_size) , 1:4)
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
  
  df_quantile = data_frame(
    time = df$time,
    conc = df$conc,
    replicate = df$replicate,
    q50 = apply(dtheo, 1, quantile, probs = 0.5, na.rm = TRUE),
    qinf95 = apply(dtheo, 1, quantile, probs = 0.025, na.rm = TRUE),
    qsup95 = apply(dtheo, 1, quantile, probs = 0.975, na.rm = TRUE)
  )
  
  if(spaghetti == TRUE){
    df_spaghetti <- as_data_frame(dtheo) %>%
      mutate(time = df$time,
             conc = df$conc,
             replicate = df$replicate)
  } else df_spaghetti <- NULL
  
  return(list(df_quantile = df_quantile,
              df_spaghetti = df_spaghetti))
  
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