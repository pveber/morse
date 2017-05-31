#' Create a list of scalars giving priors to use in Bayesian modelling (JAGS or Stan)
#'
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return A list of scalar for parameterization of priors for Bayesian modeling
#'
#' @export
#'


gm_priors = function(gm_survData){
  
  data = filter(gm_survData, time != 0)

  # Parameter calculation of concentration min and max
  conc_min = min(data$conc[data$conc != 0], na.rm=TRUE) # to remove 0 and NA
  conc_max = max(data$conc, na.rm=TRUE)
  
  time_min = min(data$time)
  time_max = max(data$time)
  
  conc_unic = sort(unique(data$conc))
  conc_unicPrec = dplyr::lag(conc_unic)
  conc_minDelta = min(conc_unic - conc_unicPrec, na.rm=TRUE)   

  ##
  ## killing rate parameter: kk
  ##
  
  kk_max = -log(0.001) / (time_min * conc_minDelta)
  kk_min = -log(0.999) / (time_max * (conc_max - conc_min))
  
  kk_meanlog10 = (log10(kk_max) + log10(kk_min)) / 2
  kk_sdlog10 = (log10(kk_max) - log10(kk_min)) / 4
  kk_taulog10 = 1 / kk_sdlog10^2
  
  ##
  ## dominant rate constant: kd
  ##
  
  kd_max = -log(0.001) / time_min
  kd_min = -log(0.999) / time_max
  
  kd_meanlog10 = (log10(kd_max) + log10(kd_min)) / 2
  kd_sdlog10 = (log10(kd_max) - log10(kd_min)) / 4
  kd_taulog10 = 1 / kd_sdlog10^2
  
  ##
  ## background hazard rate
  ##
  
  hb_max = -log(0.5) / time_min
  hb_min = -log(0.999) / time_max
  
  hb_meanlog10 = (log10(hb_max) + log10(hb_min)) / 2
  hb_sdlog10 = (log10(hb_max) - log10(hb_min)) / 4
  hb_taulog10 = 1/ hb_sdlog10^2
  
  ##
  ## non effect threshold
  ##
  
  z_meanlog10 = (log10(conc_max) + log10(conc_min)) / 2
  z_sdlog10 = (log10(conc_max) - log10(conc_min)) / 4 
  z_taulog10 = 1/ z_sdlog10^2
  
  ##
  ## non effect threshold: scale parameter & median of a log-logistic distribution
  ##
  
  alpha_meanlog10 = (log10(conc_max) + log10(conc_min)) / 2
  alpha_sdlog10 = (log10(conc_max) - log10(conc_min)) / 4 
  alpha_taulog10 = 1/ alpha_sdlog10^2

  ##
  ## shape parameter of a log-logistic distribution
  ##
  
  beta_minlog10 = -2 
  beta_maxlog10 = 2
  
  return(list(
    ## kk
    kk_meanlog10 = kk_meanlog10,
    kk_sdlog10 = kk_sdlog10,
    kk_taulog10 = kk_taulog10,
    ## kd
    kd_meanlog10 = kd_meanlog10,
    kd_sdlog10 = kd_sdlog10,
    kd_taulog10 = kd_taulog10,
    ## hb
    hb_meanlog10 = hb_meanlog10,
    hb_sdlog10 = hb_sdlog10,
    hb_taulog10 = hb_taulog10,
    ## z
    z_meanlog10 = z_meanlog10,
    z_sdlog10 = z_sdlog10,
    z_taulog10.z = z_taulog10,
    ## alpha
    alpha_meanlog10 = alpha_meanlog10,
    alpha_sdlog10 = alpha_sdlog10,
    alpha_taulog10 = alpha_taulog10,
    ## beta
    beta_minlog10 = beta_minlog10,
    beta_maxlog10 = beta_maxlog10 
  ))
}
