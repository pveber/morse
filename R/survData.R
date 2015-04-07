#' Creates a dataset for survival analysis
#'
#' This function creates a \code{survData} object from experimental data
#' provided as a \code{data.frame}. The resulting object
#' can then be used for plotting and model fitting. It can also be used
#' to generate \emph{individual-time} estimates.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame}
#' corresponds to one experimental measurement, that is a number of alive
#' individuals for a given concentration of pollutant at a certain time
#' during the assay in a certain replicate. The function fails if
#' \code{data} does not meet the
#' expected requirements. Please run \code{\link{survDataCheck}} to ensure
#' \code{data} is well-formed.
#'
#' @aliases survData print.survData summary.survData
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
#' @param firstTimeClutch the first day of reproduction activity
# FIXME @param x An object of class survData.
# FIXME @param object An object of class survData.
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return A dataframe of class \code{survData}.
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
#' @seealso \code{\link{survDataCheck}}
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
#FIXME
# print(dat)
# summary(dat)
#
#' @export
#' @importFrom dplyr left_join rename
# @importFrom plyr rename
survData <- function(data, firstTimeClutch = NULL) {
  ### INPUT
  # [data]: a [data.frame] with above mentionned requirements
  #
  ### OUTPUT
  # a [data.frame] with additional class [survData]
  # containing columns:
  # - replicate, conc, time, Nsurv as in [data]
  # - ID: a string identifier built from replicate, conc and time columns
  # - Ninit: number of initial individuals for the corresponding time series

  # test the integrity of the data with survDataCheck
  if (dim(survDataCheck(data))[1] > 0)
    stop("The [data] argument is not well-formed, please use [survDataCheck] for details.")

  data <- data[order(data$replicate, data$conc, data$time), ]

  # create an ID column of triplet replicate_conc_time
  data[, "ID"] <- idCreate(data, notime = FALSE)

  # select the start time to calculat NID
  if (is.null(firstTimeClutch)) { # default firstTimeClutch
    firstTimeClutch <- 0
  }
  
  # estimate first clutch observed time
  beforeObservedCluthTime <- min(unique(data[!unlist(lapply(split(data, data$time),
                                 function(x) {all(x$Nrepro == 0)})),
                                             "time"])) - 1
  
  if (beforeObservedCluthTime < 0) {
    beforeObservedCluthTime <- 0
  }
  
  # save original time for graph
  data$timeOrigine <- data$time
  
  if (firstTimeClutch > beforeObservedCluthTime)
    stop("The selected time is equal or higger than the first time where clutches are
         oberved !")
  
  if (firstTimeClutch != 0) { # selected firstime clutch
    if (firstTimeClutch %in% unique(data$time)) { # firstTimeClutch is one of the
      # oberserved time
      data <- data[data$time >= firstTimeClutch,]
      
    } else { # firstimeclutch is not one of the observed time
      lstPreviousTime <- max(unique(data[data$time < firstTimeClutch, "time"]))
      data <- data[c(which(data$time == lstPreviousTime),
                     which(data$time >= firstTimeClutch)),]
    }
  
    # change first time to 0
    data[data$time == min(data$time), "time"] <- 0
    # change other time
    data[data$time != 0, "time"] <- data[data$time != 0, "time"] - firstTimeClutch
  }
  
  data.t0 <- data[data$time == 0, c("replicate", "conc", "Nsurv")]
  data.t0 <- rename(data.t0, Ninit = Nsurv)
  out <- left_join(data, data.t0, by = c("replicate", "conc"))

  T <- sort(unique(out$time)) # observation times
  Nindtime <- rep(0, dim(out)[1])
  NindtimeNoCumul <- rep(0, dim(out)[1])
  for (i in 2:length(T)) {
    now <- out$time == T[i]
    before <- out$time == T[i - 1]
    Nindtime[now] <- Nindtime[before] +
      (out$Nsurv[before] - out$Nsurv[now]) * ((T[i] - T[i - 1]) / 2) +
      out$Nsurv[now] * (T[i] - T[i - 1])
    NindtimeNoCumul[now] <- Nindtime[now] - Nindtime[before]
  }

  out <- cbind(out, Nindtime, NindtimeNoCumul)

  class(out) <- c("survData", "data.frame")
  return(out)
}
