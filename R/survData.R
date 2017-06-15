#' Creates a dataset for survival analysis
#'
#' Dataset for survival analysis using Bayesian inference for \code{survFitTT},
#' \code{survFitCstC} and \code{survFitVarC} objects.
#'
#' @param \dots \code{data.frame} used to select a class

#' @export
survData <- function(...){
  
  ls.data = list(...)
  
  if(length(ls.data) == 1){
    
    survDataCstC(...)
    
  } else if(length(ls.data) == 2){
    
    survDataVarC(...)
    
  } else stop("Only 1 or 2 'data.frame' can be pass as argument.\n
              Upload 1 'data.frame' if exposure concentration are constant per
              replicate (column name should be 'replicate', 'time', 'conc', and 'Nsurv'),
              and upload 2 'data.frame's if exposure concentration 
              are variable in at least one of the replicate (with a common 'replicate' column name:\n
              - 'replicate', 'time' and 'conc' for one 'data.frame', and\n
              - 'replicate', 'time' and 'Nsurv' for the other 'data.frame'.")
}


#' Creates a dataset for survival analysis
#'
#' This function creates a \code{survDataCstC} object from experimental data
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
#' @return A dataframe of class \code{survDataCstC}.
#'
#' @seealso \code{\link{survDataCheckCstC}}
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
survDataCstC <- function(data) {
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

      class(data) <- c("survDataCstC", "survData", "data.frame")
      return(data)
}
      
#' Creates a dataset for survival analysis
#'
#' This function creates a \code{survDataVarC} object from experimental data
#' provided as two \code{data.frame}s. The resulting object can then be used for
#' plotting and model fitting.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame} corresponds to one experimental
#' measurement, that is a number of alive individuals for a given concentration
#' of pollutant at a certain time during the assay in a certain replicate. The
#' function fails if \code{data} does not meet the expected requirements. Please
#' run \code{\link{survDataCheck}} to ensure \code{data} is well-formed.
#'
#' @param data_surv a \code{data.frame} containing the following four columns:
#'   \itemize{ #' \item \code{profile}: a vector of class \code{integer},
#'   \code{character} or \code{factor} for profile identification \item
#'   \code{replicate}: a vector of class \code{integer} or \code{factor} for
#'   replicate identification \item \code{time}: a vector of class
#'   \code{integer} with time points, min value must be 0 \item \code{Nsurv}: a
#'   vector of class \code{integer} providing the number of alive individuals at
#'   each time point for each profile and each replicate }
#' @param data_conc a \code{data.frame} containing the following four columns:
#'   \itemize{ #' \item \code{profile}: a vector of class \code{integer},
#'   \code{character} or \code{factor} for profile identification \item
#'   \code{replicate}: a vector of class \code{integer} or \code{factor} for
#'   replicate identification \item \code{conc}: a vector of class
#'   \code{numeric} with tested concentrations (positive values) at each time
#'   point for each profile and each replicate \item \code{time}: a vector of
#'   class \code{integer} with time points, min value must be 0 }
#'
#'
#' @return A dataframe of class \code{survDataVarC}.
#'
#' @seealso \code{\link{survDataCheckVarC}}
#'
#' @keywords transformation
#'
#' @examples
#'
#'
#' @export

survDataVarC = function(data_x, data_y){
  
  # test the integrity of the data with survDataCheck
  if (dim(survDataCheck(data_x, data_y, diagnosis.plot = FALSE))[1] > 0)
    stop("The [data] argument is not well-formed, please use [survDataCheck] for details.")
  
  ##
  ## Join data_x and data_y
  ##
  
  df.out <- dplyr::full_join(data_x, data_y,
                             by=c("replicate","time"))
  
  class(df.out) <- c("survDataVarC", "survData", "data.frame")
  
  return(df.out)
}

