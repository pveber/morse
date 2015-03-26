#' Creates a dataset for reproduction toxicity analysis
#'
#' This function creates a \code{reproData} object from experimental
#' data provided as a \code{data.frame}. The resulting object can then be used
#' for plotting and model fitting. The \code{reproData} class is a sub-class
#' of \code{survData}, meaning that all functions and method available for
#' survival analysis can be used with \code{reproData} objects.
#'
#' The \code{data} argument contains the experimental data, and should have
#' the same structure that the argument of \code{survData}, plus a single
#' additional colum providing the total number of offspring observed since the
#' last time point. The function fails if code{data} does not meet the
#' expected requirements. Please run \code{\link{reproDataCheck}} to ensure
#' \code{data} is well-formed.
#'
#' @aliases reproData print.reproData summary.reproData
#'
#' @param data a dataframe as expected by \code{survData} containing one
#' additional \code{Nrepro} column of class \code{integer} with positive
#' values only. This column should
#' provide the number of offspring produced since the last observation.
#'
#' @param x An object of class \code{reproData}.
#' @param object An object of class \code{reproData}.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return An object of class \code{reproData}. #'
#' Generic functions:
#' \describe{
#' \item{\code{summary}}{The summary provides information about the structure
#' of the dataset and the experimental design:
#' the number of datapoints per replicate, concentration and time both for the
#' raw dataset and the transformed dataset.}
#' \item{\code{print}}{Print of a \code{reproData} object with the transformed
#' dataframe.}}
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
#' # (3) Print and summarize object dat
#' print(dat)
#' summary(dat)
#'
#' @export
#'
reproData <- function(data) {

  # test the integrity of the data with reproDataCheck
  if (dim(reproDataCheck(data))[1] > 0)
    stop("The [data] argument is not well-formed, please use [reproDataCheck] for details.")

  data <- survData(data)

  T <- sort(unique(data$time)) # observation times
  Nreprocumul <- data$Nrepro
  for (i in 2:length(T)) {
    now <- data$time == T[i]
    before <- data$time == T[i - 1]
    Nreprocumul[now] <- Nreprocumul[before] + data$Nrepro[now]
  }

  data <- cbind(data,Nreprocumul)
  class(data) <- c("reproData", "survData","data.frame")
  return(data)
}
