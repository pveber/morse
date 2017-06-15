#' Create a list of scalars giving Data to use in Bayesian modelling (JAGS or Stan)
#'
#' @param x An object of class 'survData'
#'
#' @return A list of scalar for parameterization of priors for Bayesian modeling
#'
#' @export
#'

modelData <- function(x, ...){
  UseMethod("modelData")
}


#' survDataCstC
#' 
#' 
#' 
modelData.survDataCstC <- function(x_survData, model_type = NULL){
  
  ##
  ## 0. Creation of additional variable
  ## - tprec: previous time
  ## - Nprec: previous number of survivors
  ## - time_ID: identification of row number inside a group
  ## - i_row: identification of row number (for every group)
  ## - i_row: identification of previous row number (for every group) exeptc when time_ID (in group) is 1
  
  x_survData = x_survData %>%
    # Add an indice of replicate:
    dplyr::mutate(replicate_ID = group_indices_(., .dots="replicate")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate(tprec = ifelse(time == 0, time, dply::lag(time)),
                  Nprec = ifelse( time == 0, Nsurv, dplyr::lag(Nsurv) ),
                  time_ID = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate(i_row = row_number(),
                  i_prec = ifelse(time_ID == 1, i_row, dplyr::lag(i_row)))
  
  
  ##
  ## ====== Construction of modelData
  ##
  
  # return priors for model
  priorsData = priors(x_survData, model_type)
  
  modelData = priorsData$priorsList
  
  ##
  ##  observations
  ##
  
  modelData$Nsurv = x_survData$Nsurv
  modelData$replicate = x_survData$replicate
  modelData$conc = x_survData$conc
  
  modelData$n_data = nrow(x_survData)
  
  ##
  ##  parameters
  ##
  modelData$Nprec = x_survData$Nprec
  modelData$time = x_survData$time
  
  modelData$i_prec = x_survData$i_prec
  
  ##
  ## other parameters specific to model SD vs. IT and cst vs. var
  ##
  
  if(model_type == "IT"){
    modelData$replicate_ID = x_survData$replicate_ID
    modelData$time_ID = x_survData$time_ID
  }
  
  ##
  ## =========== Object to return
  ##
  
  ##
  ## Model data Null (without observations of Nsurv)
  ##
  modelData_NULL = modelData
  modelData_NULL$Nsurv = NULL # remove Nsurv from modelData
  
  ### OUT ================
  OUT_modelDATA = list(modelData = modelData,
                       modelData_Null = modelData_NULL,
                       priorsList = priorsData$priorsList,
                       priorsMinMax = priorsData$priorsMinMax)
  
  
  class(OUT_modelDATA) <- c("modelDataCstC", "list")
  
  return(OUT_modelDATA)
  
}


#' survDataVarC
#' 
#' 
#' 
#' 
modelData.survDataVarC <- function(x_survData,
                                   model_type = NULL,
                                   extend_time = 100){
  
  ##
  ## 0. Creation of additional variable
  ## - tprec: previous time
  ## - Nprec: previous number of survivors
  ## - time_ID: identification of row number inside a group
  ## - i_row: identification of row number (for every group)
  ## - i_row: identification of previous row number (for every group) exeptc when time_ID (in group) is 1
  
  x_survData_interpolate = x_survData_interpolate(x_survData,
                                                  extend.time = extend_time) %>%
    dplyr::arrange(replicate, time)
  
  x_survData = x_survData_interpolate %>%
    dplyr::filter(!is.na(Nsurv)) %>%
    dplyr::rename(time_ID_red = time_ID_long) %>%
    # Group by replicate to replicate an indice of replicate:
    dplyr::mutate(replicate_ID = group_indices_(., .dots="replicate")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate( tprec = ifelse( time == 0, time, dplyr::lag(time) ) ) %>%
    dplyr::mutate( Nprec = ifelse( time == 0, Nsurv, dplyr::lag(Nsurv) ) ) %>%
    dplyr::mutate(time_ID = row_number()) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(i_row = row_number()) %>%
    dplyr::mutate(i_prec = ifelse(time_ID == 1, i_row, dplyr::lag(i_row))) %>%
    dplyr::arrange(replicate, time)
  
  
  # 'lag' function copy values lagged by 1 (see 'dplyr' package)
  
  
  ##
  ## ============================= Construction of modelData
  ##
  
  # return priors for model
  priorsData = gm_priors(x_survData = x_survData, model_type = model_type)
  
  modelData = priorsData$priorsList
  
  ##
  ##  observations
  ##
  
  modelData$Nsurv = x_survData$Nsurv
  modelData$replicate = x_survData$replicate
  
  ##
  ##  parameters
  ##
  modelData$Nprec = x_survData$Nprec
  
  modelData$replicate_ID = x_survData$replicate_ID
  modelData$time_ID = x_survData$time_ID
  
  
  modelData$n_dataRed = nrow(x_survData)
  modelData$n_dataLong = nrow(x_survData_interpolate)
  
  ### Integration
  modelData$replicate_ID_long  = x_survData_interpolate$replicate_ID_long
  modelData$time_ID_long = x_survData_interpolate$time_ID
  
  modelData$conc_long  = x_survData_interpolate$conc
  modelData$time_long = x_survData_interpolate$time
  
  modelData$tprec_long = x_survData_interpolate$tprec_long
  modelData$concprec_long = x_survData_interpolate$concprec_long
  
  ### Interpolation
  modelData$time_ID_red = x_survData$time_ID_red
  # modelData$i_row = x_survData$i_row
  modelData$i_prec = x_survData$i_prec
  
  ##
  ## other parameters specific to model SD vs. IT
  ##
  
  if(model_type == "SD"){
    
    modelData$tprec_ID_long = x_survData_interpolate$tprec_ID_long
    
  }
  if (model_type == "IT"){
    
    modelData$time = x_survData$time
    
  }
  
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
                       priorsMinMax = priorsData$priorsMinMax)
  
  class(OUT_modelDATA) <- c("modelDataVarC", "list")
  
  return(OUT_modelDATA)
  
}


#' Create a dataset for survival analysis when the replicate of concentration is variable
#'
#' @param x_survData An object of class \code{x_survData}
#'
#' @return A dataframe of class \code{x_survData}.
#'
#' @export
#'
#'

x_survData_interpolate = function(x_survData,
                                  extend.time = 100){
  
  ## data.frame with time
  
  df_MinMax = x_survData %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(min_time = min(time, na.rm = TRUE),
                     max_time = max(time, na.rm = TRUE)) %>%
    dplyr::group_by(replicate) %>%
    # dplyr::do(data.frame(replicate = .$replicate, time = seq(.$min_time, .$max_time, length = extend.time)))
    dplyr::do(data_frame(replicate = as.character(.$replicate), time = seq(.$min_time, .$max_time, length = extend.time)))
  
  x_survData.Interpolate = dplyr::full_join(df_MinMax,
                                            x_survData,
                                            by = c("replicate", "time")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>% # organize in replicate and time
    dplyr::mutate(conc = zoo::na.approx(conc, time, na.rm=FALSE)) %>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    dplyr::mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc)) %>%# FOR THE LAST VALUE
    # identification of time point index for Nsurv
    dplyr::mutate(id_conc_interp = ifelse(is.na(Nsurv), NA, row_number()))%>%
    # 'lag' function copy values lagged by 1 (see 'dplyr' package)
    dplyr::mutate( tprec_long = ifelse( time == 0, time, dplyr::lag(time) )
                   concprec_long = ifelse( time == 0, conc, dplyr::lag(conc) ) ) %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(time_ID_long = row_number()
                  tprec_ID_long = ifelse(time_ID_long==1, time_ID_long,  dplyr::lag(time_ID_long))) %>%
    dplyr::ungroup() %>%
    # Group by replicate to replicate an indice of replicate:
    dplyr::mutate(replicate_ID_long = group_indices_(., .dots="replicate"))
  
  
  class(x_survData.Interpolate) = c("x_survData", "data.frame")
  return(x_survData.Interpolate)
}

