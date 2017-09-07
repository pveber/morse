#' Create a dataset for survival analysis for survDataVarExp object
#'
#' @param x An object of class \code{survData}
#' @param model_type TKTD model type ('SD' or 'IT')
#'
#' @return A list
#'
#' @export
#' 
modelData.survDataVarExp <- function(x,
                                     model_type = NULL,
                                     extend_time = 100, ...){
  
  
  ## 0. Creation of additional variable
  ## - tprec: previous time
  ## - Nprec: previous number of survivors
  ## - time_ID_red: identification of row number inside a group
  ## - i_row: identification of row number (for every group)
  ## - i_prec: identification of previous row number (for every group) exept when time_ID_red (in group) is 1
  
  x_interpolate <- survData_interpolate(x,  extend_time = extend_time) %>%
    dplyr::arrange(replicate, time)
  
  x_reduce <- x_interpolate %>%
    dplyr::filter(!is.na(Nsurv)) %>%
    # Group by replicate to replicate an indice of replicate:
    dplyr::mutate(replicate_ID = group_indices_(., .dots="replicate")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate( tprec = ifelse( time == 0, time, dplyr::lag(time) ) ) %>%
    dplyr::mutate( Nprec = ifelse( time == 0, Nsurv, dplyr::lag(Nsurv) ) ) %>%
    dplyr::mutate(time_ID_red = row_number()) %>%
    dplyr::ungroup()%>%
    dplyr::mutate(i_row = row_number()) %>%
    dplyr::mutate(i_prec = ifelse(time_ID_red == 1, i_row, dplyr::lag(i_row))) %>%
    dplyr::arrange(replicate, time)
  
  ##
  ## ============================= Construction of modelData
  ##
  
  ### return priors for model
  priorsData <- priors_survData(x = x, model_type = model_type)
  
  modelData <- priorsData$priorsList
  
  ### reduce data set: To remove NA in Nsurv column
  modelData$time <-  x_reduce$time
  modelData$conc <-  x_reduce$conc
  modelData$replicate <-  x_reduce$replicate
  modelData$Nsurv <-  x_reduce$Nsurv
  modelData$Nprec <-  x_reduce$Nprec
  
  ###  parameters
  
  modelData$replicate_ID <- x_reduce$replicate_ID
  modelData$time_ID_red <- x_reduce$time_ID_red
  
  modelData$n_data_red <- nrow(x_reduce)
  
  ## Interpolation
  modelData$time_ID_long_red <- x_reduce$time_ID_long
  modelData$i_prec <- x_reduce$i_prec
  
  ## Interpolation
  
  modelData$n_data_long <- nrow(x_interpolate)
  
  ### Integration
  modelData$replicate_ID_long  <- x_interpolate$replicate_ID_long
  modelData$time_ID_long <- x_interpolate$time_ID_long
  
  modelData$conc_long <- x_interpolate$conc
  modelData$time_long <- x_interpolate$time
  modelData$replicate_long <- x_interpolate$replicate
  
  modelData$tprec_long <- x_interpolate$tprec_long
  modelData$concprec_long <- x_interpolate$concprec_long
  
  
  ##
  ## other parameters specific to model SD vs. IT
  ##
  
  if(model_type == "SD"){
    
    modelData$tprec_ID_long <- x_interpolate$tprec_ID_long
    
  }
  if (model_type == "IT"){
    
    modelData$time <- x_reduce$time
    
  }
  
  ##
  ## =========== Object to return
  ##
  

  ### OUT ================
  OUT_modelDATA <- list(modelData = modelData,
                        priorsList = priorsData$priorsList,
                        priorsMinMax = priorsData$priorsMinMax)
  
  return(OUT_modelDATA)
  
}

# Create a dataset for survival analysis when the replicate of concentration is variable
#
# @param x An object of class \code{survData}
#
# @return A dataframe
#
survData_interpolate <- function(x, extend_time = 100){
  
  ## data.frame with time
  
  df_MinMax <- x %>%
    dplyr::group_by(replicate) %>%
    dplyr::summarise(min_time = min(time, na.rm = TRUE),
                     max_time = max(time, na.rm = TRUE)) %>%
    dplyr::group_by(replicate) %>%
    # dplyr::do(data.frame(replicate = .$replicate, time = seq(.$min_time, .$max_time, length = extend_time)))
    dplyr::do(data_frame(replicate = .$replicate, time = seq(.$min_time, .$max_time, length = extend_time)))
  
  x_interpolate <- dplyr::full_join(df_MinMax, x,
                                    by = c("replicate", "time")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>% # organize in replicate and time
    dplyr::mutate(conc = zoo::na.approx(conc, time, na.rm = FALSE)) %>%
    # from package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    dplyr::mutate(conc = ifelse(is.na(conc),zoo::na.locf(conc),conc),
                  # identification of time point index for Nsurv
                  id_conc_interp = ifelse(is.na(Nsurv), NA, row_number()),
                  # 'lag' function copy values lagged by 1 (see 'dplyr' package)
                  tprec_long = ifelse( time == 0, time, dplyr::lag(time) ),
                  concprec_long = ifelse( time == 0, conc, dplyr::lag(conc) ) ) %>%
    dplyr::group_by(replicate) %>%
    dplyr::mutate(time_ID_long = row_number(),
                  tprec_ID_long = ifelse(time_ID_long==1, time_ID_long,  dplyr::lag(time_ID_long))) %>%
    dplyr::ungroup() %>%
    # Group by replicate to replicate an indice of replicate:
    dplyr::mutate(replicate_ID_long = group_indices_(., .dots="replicate"))
  
  return(x_interpolate)
}

