#' Summary for reproData objects
#' 
#' The generic \code{summary} S3 method for the \code{reproData} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param x an object of class \code{reproData}
#' @param quiet if \code{TRUE}, does no prints
#' 
#' @return The function returns a list with the following fields:
#' \item{NbrepTimeConc}{nb of replicates for all concentrations and time points}
#' \item{NbsurvTimeConc}{nb of surviving ind. for all concentrations and time points}
#' \item{NboffTimeConc}{nb of offspring for all concentrations and time points}
#' \item{NbdataConc}{nb of datapoints per concentration}
#' \item{NbdataTime}{nb of datapoints per time}
#' 
#' @seealso \code{\link{reproData}}
#' 
#' @examples
#' # (1) Load the data
#' data(cadmium1)
#' 
#' # (2) Create a reproData object
#' cadmium1 <- reproData(cadmium1)
#' 
#' # (3) Summarize the dataset
#' summary(cadmium1)
#' 
#' @keywords summary
#' 
#' @export
#' 
summary.reproData <- function(x, quiet = FALSE) {
  res <- summary.survData(x, quiet = quiet)

  # matrix of number of offspring (sum of all replicate) by time / conc
  ans3 <- tapply(x$Nrepro, list(as.factor(x$conc), as.factor(x$time)), sum)
  
  if (! quiet) {
    cat("\nNumber of offspring (sum of replicate) per time and concentration: \n")
    print(ans3)
  }
  
  invisible(c(res, list(NboffTimeConc = ans3)))
}

#' Creates a dataset for reproduction toxicity analysis
#'
#' This function creates a \code{reproData} object from experimental
#' data provided as a \code{data.frame}. The resulting object can then be used
#' for plotting and model fitting. The \code{reproData} class is a sub-class
#' of \code{survData}, meaning that all functions and method available for
#' survival analysis can be used with \code{reproData} objects.
#'
#' The \code{x} argument contains the experimental data, and should have
#' the same structure that the argument of \code{survData}, plus a single
#' additional colum providing the total number of offspring observed since the
#' last time point. The function fails if \code{x} does not meet the
#' expected requirements. Please run \code{\link{reproDataCheck}} to ensure
#' \code{x} is well-formed.
#'
#' @aliases reproData
#'
#' @param x a dataframe as expected by \code{survData} containing one
#' additional \code{Nrepro} column of class \code{integer} with positive
#' values only. This column should
#' provide the number of offspring produced since the last observation.
#'
#' @return An object of class \code{reproData}.
#'
#' @author Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#'
#' @keywords transformation
#'
#' @examples
#'
#' # (1) Load reproduction dataset
#' data(cadmium1)
#'
#' # (2) Create an object of class "reproData"
#' dat <- reproData(cadmium1)
#' dat
#' class(dat)
#'
#' @export
#'
reproData <- function(x) {

  # test the integrity of the data with reproDataCheck
  if (dim(reproDataCheck(x))[1] > 0)
    stop("The [x] argument is not well-formed, please use [reproDataCheck] for details.")

  x <- survData(x)

  T <- sort(unique(x$time)) # observation times
  Nreprocumul <- x$Nrepro
  for (i in 2:length(T)) {
    now <- x$time == T[i]
    before <- x$time == T[i - 1]
    Nreprocumul[now] <- Nreprocumul[before] + x$Nrepro[now]
  }

  x <- cbind(x,Nreprocumul)
  class(x) <- c("reproData", "survData","data.frame")
  return(x)
}
