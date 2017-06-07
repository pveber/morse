#' Create a list of scalars giving priors to use in Bayesian modelling (JAGS or Stan)
#'
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return A list of scalar for parameterization of priors for Bayesian modeling
#'
#' @export
#'


gm_priors = function(gm_survData, model_type = NULL){

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
  ## dominant rate constant: kd
  ##

  kd_max = -log(0.001) / time_min
  kd_min = -log(0.999) / time_max

  ##
  ## background hazard rate
  ##

  hb_max = -log(0.5) / time_min
  hb_min = -log(0.999) / time_max

  ##
  ## killing rate parameter: kk
  ##

  kk_max = -log(0.001) / (time_min * conc_minDelta)
  kk_min = -log(0.999) / (time_max * (conc_max - conc_min))

  ##
  ## beta
  ##

  beta_minlog10 = -2
  beta_maxlog10 = 2

  priorsMinMax= list(
    conc_min = conc_min,
    conc_max = conc_max,

    kd_min = kd_min,
    kd_max = kd_max,

    hb_min = hb_min,
    hb_max = hb_max

  )

  ##
  ## Construction of the list of priors
  ##

  priorsList =  list(
    ##
    ## dominant rate constant: kd
    ##
    kd_meanlog10 = (log10(kd_max) + log10(kd_min)) / 2 ,
    kd_sdlog10 = (log10(kd_max) - log10(kd_min)) / 4 ,
    ##
    ## background hazard rate
    ##
    hb_meanlog10 = (log10(hb_max) + log10(hb_min)) / 2 ,
    hb_sdlog10 = (log10(hb_max) - log10(hb_min)) / 4
    )

  if(model_type == "IT"){

    ## priorsMinMax
    priorsMinMax$beta_min = beta_minlog10
    priorsMinMax$beta_max = beta_maxlog10

    ## priorsList
    ### non effect threshold: scale parameter & median of a log-logistic distribution
    priorsList$alpha_meanlog10 = (log10(conc_max) + log10(conc_min)) / 2
    priorsList$alpha_sdlog10 = (log10(conc_max) - log10(conc_min)) / 4

    ### shape parameter of a log-logistic distribution
    priorsList$beta_minlog10 = beta_minlog10
    priorsList$beta_maxlog10 = beta_maxlog10

  } else if (model_type == "SD"){

    ## priorsMinMax
    priorsMinMax$kk_min = kk_min
    priorsMinMax$kk_max = kk_max

    ## priorsList
    ### killing rate parameter: kk
    priorsList$kk_meanlog10 = (log10(kk_max) + log10(kk_min)) / 2
    priorsList$kk_sdlog10 = (log10(kk_max) - log10(kk_min)) / 4
    ### non effect threshold: z
    priorsList$z_meanlog10 = (log10(conc_max) + log10(conc_min)) / 2
    priorsList$z_sdlog10 = (log10(conc_max) - log10(conc_min)) / 4
  } else stop("please, provide the 'model_type': 'IT' or 'SD'")


  return(list(priorsList = priorsList,
              priorsMinMax = priorsMinMax))
}
