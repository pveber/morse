#' Create a dataset for survival analysis for survDataCstExp object
#'
#' @param x An object of class \code{survData}
#' @param model_type TKTD model type ('SD' or 'IT')
#'
#' @return A list
#'
#' @export
#' 
modelData.survDataCstExp <- function(x, model_type = NULL){
  
  ## 1. Gather replicate when there is the same constante concentration
  x_gather <- gather_survDataCstExp(x)
  
  ## 2. Creation of additional variable require in the Bayesian model
  x_dev <- addVariable_survDataCstExp(x_gather)
  
  ##
  ## ====== Construction of modelData
  ##
  
  ### return priors for model
  
  priorsData <- priors_survData(x = x, model_type = model_type)

  ###  observations
  
  dataList <- list(
    replicate = x_dev$replicate,
    time =  x_dev$time,
    conc = x_dev$conc,
    Nsurv = x_dev$Nsurv,
    Nprec = x_dev$Nprec
  ) 

  ###  parameters
  
  dataList$n_data <- nrow(x_dev)
  
  ### other parameters specific to model IT
  
  if(model_type == "IT"){
    dataList$i_prec <- x_dev$i_prec
    dataList$replicate_ID <- x_dev$replicate_ID
    dataList$time_ID <- x_dev$time_ID
  }
  if(model_type == "SD"){
    dataList$tprec <- x_dev$tprec
  }
  
  ##
  ## =========== Object to return
  ##
  
  OUT_modelDATA <- list(dataList = dataList,
                        priorsList = priorsData$priorsList,
                        priorsMinMax = priorsData$priorsMinMax)
  
  return(OUT_modelDATA)
}


#' Gather replicates with the same concentration
#'
#' @param x An object of class \code{survData}
#'
#' @return A dataframe
#'

gather_survDataCstExp <- function(x){
  
  bool_checkTimeReplicate <- checkTimeReplicate(x)
  
  if( bool_checkTimeReplicate ){
    ### Sum Nsurv data for each (conc, time) couple
    x_dev <- x %>%
      dplyr::group_by(conc, time) %>%
      dplyr::summarise(Nsurv = sum(Nsurv)) %>%
      # concate replicate in the same replicate using factor (conc)
      dplyr::mutate(replicate = as.character(conc)) %>%
      as.data.frame()
  } else{
    x_dev <- x
  }
  return(x_dev)
}


#' Check the same number of (time, replicate) for a common concentration
#' 
#' @param x 

checkTimeReplicate <- function(x){
  df_checkTimeReplicate <- x %>%
    dplyr::group_by(conc, time) %>%
    dplyr::summarise(nbrReplicate_ConcTime = n_distinct(replicate)) %>%
    dplyr::group_by(conc) %>%
    dplyr::summarise(nbrReplicate_uniqueConc = length(unique(nbrReplicate_ConcTime)))
  
  return(all(df_checkTimeReplicate$nbrReplicate_uniqueConc == 1))
  
}
 

#' Add variables for Bayesian fitting
#'
#' @param x An object of class \code{data.frame}
#'
#' @return A dataframe
#'
addVariable_survDataCstExp <- function(x){
  ## Creation of additional variable:
  ## - tprec: previous time
  ## - Nprec: previous number of survivors
  ## - time_ID: identification of row number inside a group
  ## - i_row: identification of row number (for every group)
  ## - i_prec: identification of previous row number (for every group) except when time_ID (in group) is 1
  
  x_dev <- x %>%
    # Add an indice of replicate:
    dplyr::mutate(replicate_ID = group_indices_(., .dots = "replicate")) %>%
    dplyr::group_by(replicate) %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate(tprec = ifelse(time == 0, time, dplyr::lag(time)),
                  Nprec = ifelse( time == 0, Nsurv, dplyr::lag(Nsurv) )) %>%
    # remove time = 0 for the analysis
    # dplyr::filter(time != 0) %>%
    dplyr::mutate(time_ID = row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(replicate, time) %>%
    dplyr::mutate(i_row = row_number(),
                  i_prec = ifelse(time_ID == 1, i_row, dplyr::lag(i_row)))
  
  return(x_dev)
  
}

