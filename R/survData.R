#' Creates a dataset for survival analysis
#'
#' This function creates a \code{survData} object from experimental data
#' provided as a \code{data.frame}. The resulting object
#' can then be used for plotting and model fitting. It can also be used
#' to generate \emph{individual-time} estimates.
#'
#' Survival datasets can be in either constant or variable exposure. The
#' resulting object, in addition to its \code{survData} class, inherits the
#' class \code{survDataCstExp} or \code{survDataVarExp} respectively.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame}
#' corresponds to one experimental measurement, that is a number of alive
#' individuals for a given concentration of pollutant at a certain time
#' during the assay in a certain replicate. Note that either the concentration
#' or the number of alive individuals may be missing. The dataset is inferred
#' to be in constant exposure if the concentration is constant for each
#' replicate and systematically available. The function \code{survData} fails if
#' \code{data} does not meet the
#' expected requirements. Please run \code{\link{survDataCheck}} to ensure
#' \code{data} is well-formed.
#'
#' @param data a \code{data.frame} containing the following four columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or factor for replicate
#' identification. A given replicate value should identify the same group of
#' individuals followed in time
#' \item \code{conc}: a vector of class \code{numeric} with tested concentrations
#' (positive values, may contain NAs)
#' \item \code{time}: a vector of class \code{integer} with time points, min value must be 0
#' \item \code{Nsurv}: a vector of class \code{integer} providing the number of
#' alive individuals at each time point for each concentration and each replicate
#' (may contain NAs)
#' }
#'
#' @return A dataframe of class \code{survData}.
#'
#' @seealso \code{\link{survDataCheck}}
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
survData <- function(x) {

  # test the integrity of the data with survDataCheck
  if (dim(survDataCheck(x, diagnosis.plot = FALSE))[1] > 0)
    stop("The [x] argument is not well-formed, please use [survDataCheck] for details.")

  child_class <-
    if (is_exposure_constant(x)) "survDataCstExp"
    else "survDataVarExp"

  class(x) <- c(child_class, "survData", "data.frame")

  return(x)
}


#' Tests in a well-formed argument to function 'survData' if the concentration
#' is constant and different from NA for each replicate
#'
#' @param x a data.frame
#' @return a boolean \code{TRUE} if concentration in replicate is constant,
#'  or \code{FALSE} if the concentration in at least one of the replicates is variable
#'
#' @export
#'
is_exposure_constant <- function(x) {

  # Test if concentration is constant in a same replicate
  df.test <- x %>%
    group_by(replicate) %>%
    summarise(count = n_distinct(conc))

  has_NA <- sum(is.na(x$conc)) > 0
  return(!has_NA && all(df.test$count==1))
}


#' Computes the effective period of observation in individual days for a
#' survival dataset
#'
#' @param x an object of class \code{survData}
#' @return a numeric vector
#'
#' @export
#'
Nindtime <- function(x) {
  x <- x[!is.na(x$Nsurv),]
  T <- sort(unique(x$time)) # observation times
  Nindtime <- rep(0,dim(x)[1])
  for (i in 2:length(T)) {
    now <- x$time == T[i]
    before <- x$time == T[i - 1]
    Nindtime[now] <-
      Nindtime[before] +
      (x$Nsurv[before] - x$Nsurv[now]) * ((T[i] - T[i - 1]) / 2) +
      x$Nsurv[now] * (T[i] - T[i - 1])
  }
  return(Nindtime)
}
