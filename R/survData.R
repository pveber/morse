#' Transformed dataset for the \code{survFitTt} function from survival dataset
#' 
#' The \code{survData} function creates a \code{survData} object needed
#' to run the \code{\link{survFitTt}} function. A new dataframe called
#' \code{transformed.data} with six columns is created (see the section
#' \bold{Value}). The generic methods are \code{print} and \code{summary}.
#' 
#' The \code{survData} function builds a new dataframe called \code{transformed.data}
#' with six columns: \code{ID}, \code{replicate}, \code{conc}, \code{time},
#' \code{Ninit} and \code{Nsurv}.
#' The new dataframe is automatically reordered by \code{time}, \code{concentration}
#' and \code{replicate}. The \code{Ninit} value is the number of survivors
#' (\code{Nsurv}) at the beginning of the bioassay for the control concentration
#' for each replicate.
#' 
#' @aliases survData print.survData summary.survData
#' 
#' @param data The raw dataframe with four columns passed to the function in
#' argument \code{data}:
#' \describe{
#' \item{replicate}{A vector of class \code{integer} or factor for replicate
#' identification.}
#' \item{conc}{A vector of class \code{numeric} with tested concentrations
#' (positive values).}
#' \item{time}{A vector of class \code{integer} with time points (positive values).
#' The first time must be 0.}
#' \item{Nsurv}{A vector of class \code{integer} with positive values
#' of the number of alive individuals (positive values) at each time point for
#' each concentration and each replicate.}
#' }
#' @param x An object of class survData.
#' @param object An object of class survData.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return A dataframe of class \code{survData}. A list of two objects:
#' \item{raw.data}{The raw dataframe with four columns corresponding to the
#' argument passed in the function.}
#' \item{transformed.data}{A dataframe with six columns:
#' \describe{
#' \item{ID}{A vector of class \code{factor}. It is a character string who
#' identify the triplet \code{replicat_conc_time}.}
#' \item{replicate}{A vector of class \code{integer} or \code{factor} for
#' replicate identification.}
#' \item{conc}{A vector of class \code{numeric} with tested concentrations
#' (positive values).}
#' \item{time}{A vector of class \code{integer} with time points (positive values).}
#' \item{Ninit}{A vector of class \code{integer} with the number of individuals
#' at the beginning of the bioassay (positive values).}
#' \item{Nsurv}{A vector of class \code{integer} with positive values of the number
#' of alive individuals (positive values) for each concentration and each replicate.}
#' }}
#' 
#' Generic functions:
#' \describe{
#' \item{\code{summary}}{The summary provides information about the structure
#' of the dataset and the experimental design:
#' the number of datapoints per replicate, concentration and time both for the
#' raw dataset and the transformed dataset.}
#' \item{\code{print}}{Print of a \code{survData} object with the transformed
#' dataframe.}}
#' 
#' @author Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link{survDataCheck}}, \code{\link{survFitTt}}
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
#' # (3) Print and summarize object dat
#' print(dat)
#' summary(dat)
#' 
#' @export
#' 
survData <- function(data) {
  
  # check the data class 
  if (!is.data.frame(data))
    stop("data.frame expected!")
  
  # test the integrity of the data with survDataCheck
  if (!is.null(survDataCheck(data)$id))
    stop("There is one or more error in the data! Please use the survDataCheck function before the survData function!")

  # raw data
  raw.data <- data
  
  # reorder dataset by replicate concentration and time
  data <- data[order(data$replicate, data$conc, data$time), ]
  
  # create an ID column of triplet replicate_conc_time
  data[, "ID"] <- idCreate(data, notime = FALSE)
  
  # calcul of Ninit create transformed.data
  transformed.data <- survTransformData(data)
  
  # output of class survData
  out <- list(raw.data = raw.data,
              transformed.data = transformed.data)
  
  class(out) <- "survData"
  
  return(out)
}
