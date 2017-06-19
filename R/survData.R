#' Creates a dataset for survival analysis
#'
#' This function creates a \code{survData} object from experimental data
#' provided as a \code{data.frame}. The resulting object
#' can then be used for plotting and model fitting.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame}
#' corresponds to one experimental measurement, that is a number of alive
#' individuals for a given concentration of pollutant at a certain time
#' during the assay in a certain replicate. The function fails if
#' \code{data} does not meet the
#' expected requirements. Please run \code{\link{survDataCheckCstC}} to ensure
#' \code{data} is well-formed.
#'
#' @param data a \code{data.frame} containing the following four columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or factor for replicate
#' identification
#' \item \code{conc}: a vector of class \code{numeric} with tested concentrations
#' (positive values)
#' \item \code{time}: a vector of class \code{integer} with time points, min value must be 0
#' \item \code{Nsurv}: a vector of class \code{integer} providing the number of
#' alive individuals at each time point for each concentration and each replicate
#' }
#'
#' @return A dataframe of class \code{survDataCstC} or class \code{survDataVarC}.
#'
#' @keywords transformation
#'
#' @examples
#'
#' # (1) Load the survival dataset
#' data(zinc)
#'
#' # (2) Create an objet of class 'survData'
#' dat <- survData(zinc)
#' class(dat)
#' 
#' @importFrom dplyr left_join rename
#' 
#' @export
#' 
survData <- function(data){
  ### INPUT
  # [data]: a [data.frame] with above mentionned requirements
  #
  ### OUTPUT
  # a [data.frame] with additional class [survData]
  # containing columns:
  # - replicate, conc, time, Nsurv as in [data]
  
  # test the integrity of the data with survDataCheck
  if (dim(survDataCheck(data, diagnosis.plot = FALSE))[1] > 0)
    stop("The [data] argument is not well-formed, please use [survDataCheck] for details.")
 
  tst_cst <- test_cst_conc(data)
  
  if(tst_cst){
    class(data) <- c("survDataCstC", "survData", "data.frame")
  } else{
    class(data) <- c("survDataVarC", "survData", "data.frame")
  }

  return(data)
}

#' Test if concentration of each replicate in a 'survData' object is constant
#'
#' @param x A data.frame
#' @return a boolean \code{TRUE} if concentration in replicate is constant,
#'  or \code{FALSE} if the concentration in at least one of the replicate is variable
#'
#' @export
#' 
test_cst_conc = function(x){
  
  # Test if concentration is constant in a same replicate
  
  df.test <- x %>%
    group_by(replicate) %>%
    summarise(count = n_distinct(conc))
  
  return(all(df.test$count==1))
  
}

#' Join two datasets and creates a survData object for survival analysis
#'
#' This function join two datasets and creates a \code{survData} object from experimental data
#' provided by two \code{data.frame}. The resulting object
#' can then be used for plotting and model fitting.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame}
#' corresponds to one experimental measurement, that is a number of alive
#' individuals for a given concentration of pollutant at a certain time
#' during the assay in a certain replicate. The function fails if
#' \code{data} does not meet the expected requirements. 
#' Please run \code{\link{survDataCheck}} to ensure
#' \code{data} is well-formed.
#'
#' @param data a \code{data.frame} containing the following four columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or factor for replicate
#' identification
#' \item \code{conc}: a vector of class \code{numeric} with tested concentrations
#' (positive values)
#' \item \code{time}: a vector of class \code{integer} with time points, min value must be 0
#' \item \code{Nsurv}: a vector of class \code{integer} providing the number of
#' alive individuals at each time point for each concentration and each replicate
#' }
#'
#' @return A dataframe of class \code{survDataCstC} or class \code{survDataVarC}.
#'
#' @keywords transformation
#'
#' @examples
#'
#' # (1) Load the two survival dataset
#' data(propiconazole_varC_conc) ; data(propiconazole_varC_Nsurv)
#' 
#' # (2) Create an objet of class 'survData'
#' dat_join <- survData_join(propiconazole_varC_conc,propiconazole_varC_Nsurv)
#' class(dat_join)
#' 
#' @importFrom dplyr left_join rename
#' 
#' @export
#' 
survData_join(data_x, data_y){
  
  ##
  ## 1 check column number
  ##
  
  if ( ncol(data_x) != 3  | ncol(data_y) != 3  ) {
    stop("Object must have 3 columns: 'replicate','time' and 'Nsurv' and 'replicate','time' and 'conc'")
  }
  
  ##
  ## 2 assert column names are correct and assign join data_x or data_y to corresponding data.frame
  ##
  ref.names_surv <- c("replicate","time","Nsurv")
  ref.names_conc <- c("replicate","time","conc")
  
  ### data_surv
  if( (all(colnames(data_x) %in% ref.names_surv) & all(colnames(data_y) %in% ref.names_conc))
      ||
      (all(colnames(data_x) %in% ref.names_conc) & all(colnames(data_y) %in% ref.names_surv))){
    
    
    join_df <- dplyr::full_join(data_x, data_y,
                                by=c("replicate","time")) %>%
      dplyr::arrange(replicate, time)
    
    
  } else {
    stop("'data_x' and 'data_y' must both have 3 columns: one with 'replicate','time' and 'Nsurv'
         and the other with 'replicate','time' and 'conc'.")
    
  }
  
  OUT_df <- survData(join_df) 
  
  return(OUT_df)
  
}
