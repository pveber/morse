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
#' @keywords transformation
#'
#' @examples
#'
#' # (1) Load reproduction dataset
#' data(cadmium1)
#'
#' # (2) Create an object of class "reproData"
#' dat <- reproData(cadmium1)
#' class(dat)
#'
#' @export
#' 
reproData <- function(data) {

  # test the integrity of the data with reproDataCheck
  if (dim(reproDataCheck(data, diagnosis.plot = FALSE))[1] > 0)
    stop("The [x] argument is not well-formed, please use [reproDataCheck] for details.")

  data <- survData(data)
  data <- data[order(data$replicate, data$conc, data$time), ]
  
  data.t0 <- data[data$time == 0, c("replicate", "conc", "Nsurv")]
  data.t0 <- rename(data.t0, Ninit = Nsurv)
  out <- left_join(data, data.t0, by = c("replicate", "conc"))
  
  Tab <- sort(unique(data$time)) # observation times
  Nindtime <- rep(0,dim(out)[1])
  for (i in 2:length(Tab)) {
    now <- out$time == Tab[i]
    before <- out$time == Tab[i - 1]
    Nindtime[now] <-
      Nindtime[before] +
      (out$Nsurv[before] - out$Nsurv[now]) * ((Tab[i] - Tab[i - 1]) / 2) +
      out$Nsurv[now] * (Tab[i] - Tab[i - 1])
  }
  
  x <- cbind(out, Nindtime)
  # force concentration as type double
  x$conc <- as.double(out$conc)
  Tab_2 <- sort(unique(x$time)) # observaTab_2ion Tab_2imes
  Nreprocumul <- x$Nrepro
  for (i in 2:length(Tab_2)) {
    now <- x$time == Tab_2[i]
    before <- x$time == Tab_2[i - 1]
    Nreprocumul[now] <- Nreprocumul[before] + x$Nrepro[now]
  }

  x <- cbind(x,Nreprocumul)
  # force concentration as type double
  x$conc <- as.double(x$conc)
  class(x) <- c("reproData", "survData", "data.frame")
  return(x)
}
