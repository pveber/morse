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
    dplyr::do(data.frame(profile = .$profile, time = seq(.$min_time, .$max_time, length = extend.time)))

  gm_survData.Interpolate = dplyr::full_join(df_MinMax,
                      gm_survData,
                      by = c("profile", "time")) %>%
    dplyr::group_by(profile) %>%
    dplyr::arrange(time) %>%
    dplyr::mutate(conc.origin = conc) %>%
    dplyr::mutate(conc = na.approx(conc,time, na.rm=FALSE)) %>%
    dplyr::mutate(id_conc_interp = ifelse(is.na(Nsurv), NA, row_number()))%>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    dplyr::mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc)) %>%# FOR THE LAST VALUE
    dplyr::arrange(profile, time) %>%
    dplyr::ungroup() %>%
    # Group by profile to profile an indice of profile:
    dplyr::mutate(profile_ID_long = group_indices_(., .dots="profile")) %>%
    dplyr::group_by(profile) %>%
    dplyr::mutate(time_ID_long = row_number()) %>%
    dplyr::ungroup()

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
                        model_type = NULL){

  test_cst_conc = test_cst_conc(gm_survData = gm_survData)

  if(test_cst_conc){

    # return priors for model.
    priorsData = gm_priors(gm_survData = gm_survData, model_type = model_type)

    gm_survData = gm_survData %>%
      # Group by profile to profile an indice of profile:
      dplyr::mutate(profile_ID = group_indices_(., .dots="profile")) %>%
      dplyr::group_by(profile) %>%
      dplyr::mutate(time_ID = row_number()) %>%
      dplyr::ungroup()


  } else{

    gm_survData_interpolate = df.long_interpolate(gm_survData,
                                   extend.time=extend_time)

    gm_survDataRed = gm_survData_interpolate %>%
      dplyr::filter(!is.na(Nsurv)) %>%
      dplyr::rename(time_ID_red = time_ID_long) %>%
      # Group by profile to profile an indice of profile:
      dplyr::mutate(profile_ID = group_indices_(., .dots="profile")) %>%
      dplyr::group_by(profile) %>%
      dplyr::mutate(time_ID = row_number()) %>%
      dplyr::ungroup()
  }

  ##
  ## ============================= Construction of modelData
  ##

  modelData = priorsData

  ##
  ##  observations
  ##
  modelData$conc = gm_survData$conc
  modelData$Nsurv = gm_survData$Nsurv

  ##
  ##  parameters
  ##
  modelData$Nprec = gm_survData$Nprec
  modelData$time = gm_survData$time

  modelData$profile_ID = gm_survData$profile_ID
  modelData$time_ID = gm_survData$time_ID

  modelData$n_data = nrow(gm_survData)

  ##
  ## other parameters specific to model SD vs. IT and cst vs. var
  ##

  if(model_type == "IT"){
    if(test_cst_conc){ # return TRUE if conc are constant by replicat, FALSE otherwise

    } else{

      # Nécessaire lorsque la concentration est variable en fonction du temps
      modelData$intC = MYintCList$intC
      modelData$intC.time = MYintCList$intC.time
      modelData$id.intCtime = MYintCList$id.intCtime
      modelData$N.intC = MYintCList$N.intC

    }
  } else if (model_type == "SD"){
    if(test_cst_conc){

      modelData$bigtime = max(gm_survData$time)+10

    } else{

      # Nécessaire lorsque la concentration est variable en fonction du temps
      modelData$intC = MYintCList$intC
      modelData$intC.time = MYintCList$intC.time
      modelData$id.intCtime = MYintCList$id.intCtime
      modelData$N.intC = MYintCList$N.intC

    }
  } else stop("Please provide 'model_type': 'SD' or 'IT'")

  ##
  ## =========== Object to return
  ##

  ##
  ## Model data Null (without observations of Nsurv)
  ##
  modelData_NULL = modelData
  modelData_NULL$Nsurv = NULL

  return(list(modelData = modelData,
              modelData_Null = modelData_NULL,
              priorsData = priorsData))

}
