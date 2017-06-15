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
