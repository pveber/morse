#' Summarize of reproData object
#' 
#' The generic \code{summary} S3 method for the \code{reproData} class provides
#' information about the structure of the dataset and the experimental design.
#' 
#' @param x An object of class \code{reproData}
#' @param quiet If \code{TRUE}, make silent all prints
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
  # matrix of number of replicate by time / conc
  ans1 <- table(x[, c("conc", "time")])
  
  # matrix of number of survival (sum of all replicate) by time / conc
  ans2 <- tapply(x$Nsurv, list(as.factor(x$conc), as.factor(x$time)), sum)
  
  # matrix of number of offspring (sum of all replicate) by time / conc
  ans3 <- tapply(x$Nrepro, list(as.factor(x$conc), as.factor(x$time)), sum)
  
  # table of datapoints per concentration
  ans4 <- table(x$conc)
  
  # table of datapoints per time
  ans5 <- table(x$time)
  
  if (! quiet) {
    cat("Number of replicate per time and concentration: \n")
    print(ans1)
    cat("\nNumber of survival (sum of replicate) per time and concentration: \n")
    print(ans2)
    cat("\nNumber of offspring (sum of replicate) per time and concentration: \n")
    print(ans3)
    cat("\nNumber of datapoints per concentration: \n")
    print(ans4)
    cat("\nNumber of datapoints per time: \n")
    print(ans5)
  }
  
  invisible(list(NbrepTimeConc = ans1,
                 NbsurvTimeConc = ans2,
                 NboffTimeConc = ans3,
                 NbdataConc = ans4,
                 NbdataTime = ans5))
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
