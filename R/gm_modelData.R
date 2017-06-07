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

      modelData$tprec = gm_survData$tprec
      modelData$conc = gm_survData$conc

      modelData$bigtime = max(gm_survData$time)+10

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
