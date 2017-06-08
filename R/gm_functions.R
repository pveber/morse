################################################################################
#
#  DATACHECK
#
################################################################################

#' Checks if an object can be used to perform survival analysis
#'
#' The \code{gm_survDataCheck} function can be used to check if an object
#' containing survival data is formatted according to the expectations of the
#' \code{gm_survData} function.
#'
#'
#' @aliases gm_survDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, the function may produce diagnosis plots
#'
#' @return
#' @seealso \code{\link{survData}}
#'
#' @examples
#' # Run the check data function
#' data(zinc)
#' gm_survDataCheck(zinc)
#'
#' # Now we insert an error in the dataset, by artificially increasing the
#' # number of survivors at some time point, in such a way that the number
#' # of indivuals increases in some replicate
#' zinc[25, "Nsurv"] <- as.integer(20)
#' gm_survDataCheck(zinc)
#'
#'
#' @export
gm_survDataCheck <- function(data_surv, data_conc) {

  ##
  ## 0. check if data.frame
  ##
  if(!("data.frame" %in% class(data_surv) & "data.frame" %in% class(data_conc))) {
    stop("'data_surv' and/or 'data_conc' are not of class 'data.frame'.")
  }

  ##
  ## 1. check column names
  ##

    ### 1.1 data_surv
    ref.namesSurv <- c("profile","replicate","time","Nsurv")

    if ( ncol(data_surv) != 4  # check same number of column
         | !(all(colnames(data_surv) %in% ref.namesSurv))  ) { # check if the name are the same, order doesn't matter
      stop( "'data_surv' must have 4 columns with names: 'profile','replicate','time' and 'Nsurv' (order does not matter)." )
    }

    ### 1.2 data_conc
    ref.namesConc <- c("profile","time","conc")

    if ( ncol(data_conc) != 3  # check same number of column
         | !(all(colnames(data_conc) %in% ref.namesConc))  ) { # check if the name are the same, order doesn't matter
      stop( "'data_conc' must have 3 columns with names: 'profile', 'time' and 'conc' (order does not matter)." )
    }

  ##
  ## 2. assert the first time point is zero for each (profile, replicate)
  ##
  df_check0Surv = data_surv %>%
    group_by(profile,replicate) %>%
    summarize(check_0 = 0 %in% time)
  if(all(df_check0Surv$check_0) != TRUE ){
    stop( "In each '(profile, replicate)' of 'data_surv', the first time point must be '0'." )
  }

  df_check0Conc = data_conc %>%
    group_by(profile) %>%
    summarize(check_0 = 0 %in% time)

  if(all(df_check0Conc$check_0) !=TRUE ){
    stop( "In each 'profile' of 'data_conc', the first time point must be '0'." )
  }

  ##
  ## 3. time are numeric, Nsurv are interger, assert concentrations are numeric
  ##
  ### 3.1 data_surv
  if (!is.integer(data_surv$Nsurv) | any(data_surv$Nsurv < 0)
      | !is.numeric(data_surv$time) | any(data_surv$time < 0)      ) {
    stop( "In table 'data_surv', column 'time' must contain only positive numeric values and column 'Nsurv' must contain only positive integer values. " )
  }
  ### 3.2 data_conc
  if (!is.numeric(data_conc$conc) | any(data_conc$conc < 0)
      | !is.numeric(data_conc$time) | any(data_conc$time < 0)      ) {
    stop( "In table 'data_conc', column 'time' must contain only positive numeric values and column 'conc' must contain only positive integer values. " )
  }


  ##
  ## 4. assert Nsurv != 0 at time 0
  ##

  df_checkNeq0 = data_surv %>%
    filter(time == 0) %>%
    mutate( check_NsurvNeq0 = (Nsurv != 0))

  if(all(df_checkNeq0$check_NsurvNeq0) != TRUE){
    stop( "In each '(profile, replicate)', 'Nsurv' must be greater than '0'." )
  }

  # ##
  # ## 5. assert there is the same name and/or number of replicates between tables
  # ##
  #
  # df_checkReplicate_surv = data_surv %>%
  #   group_by(profile) %>%
  #   summarize( check_Replicate = length(unique(replicate)))
  #
  # df_checkReplicate_conc = data_conc %>%
  #   group_by(profile) %>%
  #   summarize( check_Replicate = length(unique(replicate)))
  #
  # if( all( df_checkReplicate_surv$check_Replicate == df_checkReplicate_conc$check_Replicate ) != TRUE){
  #   stop( "At least one 'replicate' is missing or different between tables 'data_surv' and 'data_conc'." )
  # }

  ##
  ## 6. assert there is the same name and/or number of profile between tables
  ##

  df_checkProfile_surv = data_surv %>%
    summarize( check_Profile = length(unique(profile)))

  df_checkProfile_conc = data_conc %>%
    summarize( check_Profile = length(unique(profile)))

  if( all( df_checkProfile_surv$check_Profile == df_checkProfile_conc$check_Profile ) != TRUE){
    stop( "At least one 'profile' is missing or different between tables 'data_surv' and 'data_conc'." )
  }

  ##
  ## 7. assert Nsurv never increases with time
  ##

  df_checkSurvIncrease = data_surv %>%
    group_by(profile, replicate) %>%
    arrange(time) %>%
    mutate(Nprec = ifelse(time == 0, Nsurv, lag(Nsurv))) %>%
    mutate( check_SurvIncrease = Nsurv <= Nprec)

  if(all(df_checkSurvIncrease$check_SurvIncrease) != TRUE){
    stop( "'Nsurv' increases at some time points." )
  }


  # ##
  # ## 8. check if the profile of concentration are the same for every `replicate` of a same `profile`
  # ##
  #
  # df_checkProfile = data_conc %>%
  #   group_by(profile, time) %>%
  #   summarize(checkProfile = length(unique(conc)))
  #
  # if(any(df_checkProfile$checkProfile) != 1){
  #   stop( "Elements in 'conc' must be the same between 'replicate' of a same '(profile, time)'." )
  # }

  ##
  ## 9. WARNING - assert max(time in data_conc) >= max(time in data_surv)
  ##

  df_checkMaxTimeSurv = data_surv %>%
    group_by(profile) %>%
    filter(time == max(time))

  df_checkMaxTimeConc = data_conc %>%
    group_by(profile) %>%
    filter(time == max(time))

  if(!all(df_checkMaxTimeConc$time >= df_checkMaxTimeSurv$time) ){
    warning( "In each 'profile', maximum time in 'data_conc' should be greater or equal to maximum time in 'data_surv'.
             Otherwise, last concentration is taken to fill concentration profile until the maximum time in 'data_surv'." )
  }

  ##
  ## 10. for every profile, 'time' must be the same in 'data_surv'
  ##

  df_checkTimeSurv = data_surv %>%
    mutate(concate = paste0(profile,replicate)) %>%
    mutate(time_ref = time)%>%
    select(concate, time, time_ref) %>%
    spread(key = concate, value = time)

  if(length(unique(as.list(df_checkTimeSurv))) != 1) {
    stop( "In 'data_surv' all 'time' vector must be the same." )
  }
}

################################################################################
#
#  PRIORS
#
################################################################################

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


################################################################################
#
#  MODELDATA
#
################################################################################

#' Test if concentration of profile in 'gm_survData' is constant
#'
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return a boolean \code{TRUE} if concentration in profile is constant, or \code{FALSE} if the concentration in at least one of the profile is variable
#'
#' @export

test_cst_conc = function(gm_survData){
  df.test = gm_survData %>%
    group_by(profile)%>%
    summarise(count = n_distinct(conc))
  return(all(df.test$count==1))
}

#' Create a dataset for survival analysis when the profile of concentration is variable
#'
#' @param gm_survData An object of class \code{gm_survData}
#'
#' @return A dataframe of class \code{gm_survData}.
#'
#' @export
#'
#'

gm_survData_interpolate = function(gm_survData,
                                   extend.time = 100){

  ## data.frame with time

  df_MinMax = gm_survData %>%
    dplyr::group_by(profile) %>%
    dplyr::summarise(min_time = min(time, na.rm = TRUE),
              max_time = max(time, na.rm = TRUE)) %>%
    dplyr::group_by(profile) %>%
    # dplyr::do(data.frame(profile = .$profile, time = seq(.$min_time, .$max_time, length = extend.time)))
    dplyr::do(data_frame(profile = as.character(.$profile), time = seq(.$min_time, .$max_time, length = extend.time)))

  gm_survData.Interpolate = dplyr::full_join(df_MinMax,
                      gm_survData,
                      by = c("profile", "time")) %>%
    dplyr::group_by(profile) %>%
    dplyr::arrange(profile, time) %>% # organize in profile and time
    dplyr::mutate(conc.origin = conc) %>%
    dplyr::mutate(conc = na.approx(conc,time, na.rm=FALSE)) %>%
    dplyr::mutate(id_conc_interp = ifelse(is.na(Nsurv), NA, row_number()))%>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    dplyr::mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc)) %>%# FOR THE LAST VALUE
    # 'lag' function copy values lagged by 1 (see 'dplyr' package)
    dplyr::mutate( tprec_long = ifelse( time == 0, time, dplyr::lag(time) ) ) %>%
    dplyr::mutate( concprec_long = ifelse( time == 0, conc, dplyr::lag(conc) ) ) %>%
    dplyr::group_by(profile) %>%
    dplyr::mutate(time_ID_long = row_number()) %>%
    dplyr::mutate(tprec_ID_long = ifelse(time_ID_long==1, time_ID_long,  dplyr::lag(time_ID_long))) %>%
    dplyr::ungroup() %>%
    # Group by profile to profile an indice of profile:
    dplyr::mutate(profile_ID_long = group_indices_(., .dots="profile"))


  class(gm_survData.Interpolate) = c("gm_survData", "data.frame")
  return(gm_survData.Interpolate)
}

#' Create a list of scalars giving priors to use in Bayesian modelling (JAGS or Stan)
#'
#' @param gm_survData An object of class 'gm_survData'
#'
#' @return A list of scalar for parameterization of priors for Bayesian modeling
#'
#' @export
#'


gm_modelData = function(gm_survData,
                        model_type = NULL,
                        extend_time = 100){

  cst_conc = test_cst_conc(gm_survData)

  # return priors for model
  priorsData = gm_priors(gm_survData = gm_survData, model_type = model_type)

  if(cst_conc){

    gm_survData = gm_survData %>%
      # Group by profile to profile an indice of profile:
      dplyr::mutate(profile_ID = group_indices_(., .dots="profile")) %>%
      dplyr::group_by(profile) %>%
      dplyr::arrange(profile, time) %>%
      dplyr::mutate(tprec = ifelse(time == 0, time, lag(time))) %>%
      dplyr::mutate(time_ID = row_number()) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(profile, time) %>%
      dplyr::mutate(i_row = row_number()) %>%
      dplyr::mutate(i_prec = ifelse(time_ID == 1, i_row, dplyr::lag(i_row)))



  } else{

    gm_survData_interpolate = gm_survData_interpolate(gm_survData,
                                                      extend.time = extend_time) %>%
      dplyr::arrange(profile, time)

    gm_survData = gm_survData_interpolate %>%
      dplyr::filter(!is.na(Nsurv)) %>%
      dplyr::rename(time_ID_red = time_ID_long) %>%
      # Group by profile to profile an indice of profile:
      dplyr::mutate(profile_ID = group_indices_(., .dots="profile")) %>%
      dplyr::group_by(profile) %>%
      dplyr::mutate(time_ID = row_number()) %>%
      dplyr::ungroup()%>%
      dplyr::mutate(i_row = row_number()) %>%
      dplyr::mutate(i_prec = ifelse(time_ID == 1, i_row, dplyr::lag(i_row))) %>%
      dplyr::arrange(profile, time)
  }

  ##
  ## ============================= Construction of modelData
  ##

  modelData = priorsData$priorsList

  ##
  ##  observations
  ##

  modelData$Nsurv = gm_survData$Nsurv
  modelData$profile = gm_survData$profile

  ##
  ##  parameters
  ##
  modelData$Nprec = gm_survData$Nprec
  modelData$time = gm_survData$time

  modelData$profile_ID = gm_survData$profile_ID
  modelData$time_ID = gm_survData$time_ID

  ##
  ## other parameters specific to model SD vs. IT and cst vs. var
  ##

  if(model_type == "IT"){
    if(cst_conc){ # return TRUE if conc are constant by replicat, FALSE otherwise

      modelData$n_data = nrow(gm_survData)

      modelData$conc = gm_survData$conc
      modelData$i_prec = gm_survData$i_prec

    } else{

      modelData$n_dataRed = nrow(gm_survData)
      modelData$n_dataLong = nrow(gm_survData_interpolate)

      ### Integration
      modelData$profile_ID_long  = gm_survData_interpolate$profile_ID_long
      modelData$time_ID_long = gm_survData_interpolate$time_ID_long
      #modelData$tprec_ID_long = gm_survData_interpolate$tprec_ID_long
      modelData$conc_long  = gm_survData_interpolate$conc
      modelData$time_long = gm_survData_interpolate$time

      modelData$tprec_long = gm_survData_interpolate$tprec_long
      modelData$concprec_long = gm_survData_interpolate$concprec_long

      ### Interpolation
      modelData$time_ID_red = gm_survData$time_ID_red
      # modelData$i_row = gm_survData$i_row
      modelData$i_prec = gm_survData$i_prec

    }
  } else if (model_type == "SD"){
    if(cst_conc){

      modelData$n_data = nrow(gm_survData)

      modelData$conc = gm_survData$conc

      modelData$i_prec = gm_survData$i_prec
      modelData$profile_ID = NULL
      modelData$time_ID = NULL

    } else{

      modelData$n_dataRed = nrow(gm_survData)
      modelData$n_dataLong = nrow(gm_survData_interpolate)

      ### Integration
      modelData$profile_ID_long  = gm_survData_interpolate$profile_ID_long
      modelData$time_ID_long = gm_survData_interpolate$time_ID
      modelData$tprec_ID_long = gm_survData_interpolate$tprec_ID_long
      modelData$conc_long  = gm_survData_interpolate$conc
      modelData$time_long = gm_survData_interpolate$time

      modelData$tprec_long = gm_survData_interpolate$tprec_long
      modelData$concprec_long = gm_survData_interpolate$concprec_long

      ### Interpolation
      modelData$time_ID_red = gm_survData$time_ID_red
      # modelData$i_row = gm_survData$i_row
      modelData$i_prec = gm_survData$i_prec

    }
  } else stop("Please provide 'model_type': 'SD' or 'IT'")

  ##
  ## =========== Object to return
  ##

  ##
  ## Model data Null (without observations of Nsurv)
  ##
  modelData_NULL = modelData
  modelData_NULL$Nsurv = NULL # remove Nsurv from mdelData

  ### OUT ================
  OUT_modelDATA = list(modelData = modelData,
                       modelData_Null = modelData_NULL,
                       priorsList = priorsData$priorsList,
                       priorsMinMax = priorsData$priorsMinMax,
                       cst_conc = cst_conc)


  class(OUT_modelDATA) <- c("gm_modelData", "list")

  return(OUT_modelDATA)

}
