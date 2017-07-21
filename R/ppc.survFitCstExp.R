#' Posterior predictive check plot for \code{survFitCstExp} objects
#'
#' This is the generic \code{ppc} S3 method for the \code{survFitTKTD} class. It
#' plots the predicted values along with 95\% credible intervals
#' versus the observed values for \code{survFit} objects.
#'
#' The black points show the observed number of survivors (pooled
#' replicates, on \eqn{X}-axis) against the corresponding predicted
#' number (\eqn{Y}-axis). Predictions come along with 95\% prediction
#' intervals, which are depicted in green when they contain the
#' observed value and in red otherwise. Samples with equal observed
#' value are shifted on the \eqn{X}-axis. For that reason, the
#' bisecting line (y = x), is represented by steps when observed
#' values are low. That way we ensure green intervals do intersect the
#' bisecting line.
#'
#' @param x An object of class \code{survFit}
#' @param style graphical backend
#' @param \dots Further arguments to be passed to generic methods
#'
#' @import ggplot2
#' @import grDevices
#' @importFrom graphics plot
#'
#' @export
#' 
ppc.survFitCstExp <- function(x, style = "ggplot", ...) {

  xlab <- "Observed nb of survivors"
  ylab <- "Predicted nb of survivors"
  
  ppc_gen(EvalsurvTKTDPpc_CstExp(x), style, xlab, ylab)
}

#' @importFrom stats rbinom quantile
EvalsurvTKTDPpc_CstExp <- function(x) {
  tot.mcmc <- do.call("rbind", x$mcmc)

  model_type <- x$model_type
  
  kd <- 10^(tot.mcmc[, "kd_log10"])
  hb <- 10^(tot.mcmc[, "hb_log10"])
  if(model_type == "SD"){
    z <- 10^(tot.mcmc[, "z_log10"])
    kk <- 10^(tot.mcmc[, "kk_log10"])
  } else if (model_type == "IT"){
    alpha <- 10^(tot.mcmc[, "alpha_log10"])
    beta <- 10^(tot.mcmc[, "beta_log10"])
  } else{
    stop("'model_type must be 'SD' or 'IT'")
  }
  
  
  if(model_type == "SD"){
    
    niter <- nrow(tot.mcmc)
    #n <- x$jags.data$ndat
    n <- x$jags.data$n_data
    
    #xconc <- x$jags.data$x
    xconc <- x$jags.data$conc
    
    #t <- x$jags.data$t
    t <- x$jags.data$time
    
    tprec <- x$jags.data$tprec
    
    #NsurvObs <- x$jags.data$y
    NsurvObs <- x$jags.data$Nsurv
    
    Nprec <- x$jags.data$Nprec
    
    NsurvPred <- matrix(NA, nrow = niter, ncol = n)
    psurv = NULL
    
    # bigtime <- x$jags.data$bigtime
    bigtime <- max(x$jags.data$time) + 10
    
    for (i in 1:n) {
      for (j in 1:length(kd)) {
        xcor <- ifelse(xconc[i] > 0, xconc[i], 10)
        R <- ifelse(xconc[i] > z[j], z[j]/xcor, 0.1)
        tz <- ifelse(xconc[i] > z[j], -1 / kd[j] * log(1 - R), bigtime)
        tref <- max(tprec[i], tz)
        psurv[j] <- exp(-hb[j] * (t[i] - tprec[i]) +
                          if (t[i] > tz) {
                            -kk[j] * ((xconc[i] - z[j]) * (t[i] - tref) +
                                        xconc[i]/kd[j] * (exp(-kd[j] * t[i]) - exp(-kd[j] * tref)))
                          } else {
                            0
                          })
      }
      NsurvPred[, i] <- rbinom(niter, Nprec[i], psurv)
    }
    NsurvPred <- t(NsurvPred)
  }
  if(model_type == "IT"){
    # all theorical
    conc.obs = unique(x$jags.data$conc)
    time.obs = unique(x$jags.data$time)
    
    k <- 1:length(conc.obs)
    j <- 1:length(time.obs)
    
    dtheo.IT <- lapply(k, function(kit) { # concentration pour chaque concentration
      Surv_IT(Cw = conc.obs[kit],
              time = time.obs,
              kd = kd,
              hb = hb,
              alpha = alpha,
              beta = beta)
    })
    # transpose dtheo
    dtheo <- do.call("rbind", lapply(dtheo.IT, t))
    
    Nprec = x$jags.data$Nprec
    i_prec = x$jags.data$i_prec
    ncol_NsurvPred = ncol(dtheo)
    
    NsurvPred = matrix(NA, ncol = ncol_NsurvPred, nrow = nrow(dtheo))
    
    for(i in 1:nrow(dtheo)){
      NsurvPred[i, ] = rbinom(ncol_NsurvPred, size = Nprec[i], prob = as.numeric(dtheo[i,] / dtheo[i_prec[i],]))
    }
    
  }
  
  QNsurvPred <- t(apply(NsurvPred, 1, quantile,
                        probs = c(2.5, 50, 97.5) / 100, na.rm = TRUE))
  
  tab <- data.frame(QNsurvPred,
                    Nprec,
                    NsurvObs,
                    col = ifelse(QNsurvPred[,"2.5%"] > NsurvObs |
                                   QNsurvPred[,"97.5%"] < NsurvObs,
                                 "red", "green"))
  colnames(tab) <- c("P2.5", "P50", "P97.5", "Nprec", "Obs", "col")
  
  return(tab)
}

# ppc.survFit <- function(x, ...) {
#   
#   ### compute posteriors median and 95 CI
#   modelData = x$modelData
#   postData = posteriorData(x$mcmc, model_type = x$model_type)
#   
#   df_plt = data_frame(Nsurv = modelData$Nsurv,
#                       time = modelData$time,
#                       replicate = modelData$replicate,
#                       Nsurv_q50 = apply(postData$df_ppc, 2, quantile, probs = 0.5, na.rm = TRUE),
#                       Nsurv_qinf95 = apply(postData$df_ppc, 2, quantile, probs = 0.025, na.rm = TRUE),
#                       Nsurv_qsup95 = apply(postData$df_ppc, 2, quantile, probs = 0.975, na.rm = TRUE)) %>%
#     mutate(col_range = ifelse(Nsurv_qinf95 > Nsurv | Nsurv_qsup95 < Nsurv, "out", "in"))
#   
#   
#   ppc_plt = df_plt %>%
#     ggplot() + theme_bw() +
#     theme(legend.position="none") +
#     # expand_limits(x = 0, y = 0) +
#     scale_colour_manual(values = c("green", "red")) +
#     scale_x_continuous(name="Observed nb of survivors") +
#     scale_y_continuous(name="Predicted nb of survivors") +
#     geom_abline(slope = 1) +
#     geom_linerange(aes(x = Nsurv,
#                        ymin = Nsurv_qinf95,
#                        ymax = Nsurv_qsup95 ,
#                        group = replicate,
#                        color = col_range),
#                    position = position_dodge(width=0.5)) +
#     geom_point(aes(x = Nsurv,
#                    y = Nsurv_q50,
#                    group = replicate ),
#                position = position_dodge(width=0.5))
#   
#   return(ppc_plt)
#   
# }

