#' Transformed dataset for the \code{reproFitTt} function from reprodution dataset
#' 
#' The \code{reproData} function creates a \code{reproData} object
#' needed to run the \code{\link{reproFitTt}} function. A new dataframe
#' called \code{transformed.data} with nine columns is created
#' (see the section \bold{Value}). The generic methods are \code{print} and
#' \code{summary}.
#' 
#' The \code{reproData} function builds a new dataframe called
#' \code{transformed.data}, with nine columns: \code{ID}, \code{replicate},
#' \code{conc}, \code{time}, \code{Ninit}, \code{Nsurv}, \code{Nrepro},
#' \code{Nreprocumul} and \code{Nindtime}. The new dataframe is automatically
#' reordered by \code{time}, \code{replicate} and \code{concentration}.
#' The number of individual-days (\code{Nindtime}) is needed to account for the
#' time-contribution of each individual to the cumulative reproduction
#' (Delignette-Muller et al., 2014).
#' The \code{Ninit} value is the number of survivors (\code{Nsurv}) at the beginning
#' of the bioassay for the control concentration for each replicate. The
#' \code{Nreprocumul} value is the cumulative number of offsrping calculated by:
#' 
#' \deqn{Nreprocumul_{1} = Nrepro_{1}}{Nreprocumul_{1} = Nrepro_{1}}
#' 
#' \deqn{Nreprocumul_{t} = Nreprocumul_{t-1} + Nrepro_{t}}{Nreprocumul_{t} =
#' Nreprocumul_{t-1} + Nrepro_{t}} with \eqn{t \in [2, target.time]}{t in [2,
#' target.time]}.
#' 
#' The \code{Nindtime} value is the number of individual-days at the time point.
#' Using the survival data, it is thus possible to calculate the period during
#' which each parent animal has stay alive, or the period during which it may
#' have reproduced. For such a calculation, an animal discovered dead at time point
#' \eqn{t_{i+1}}{t_{i+1}} was assumed to be alive from the beginning of the
#' experiment to time \eqn{\frac{t_{i+1} + t_{i}}{2}}{(t_{i+1} + t_{i}) / 2}.
#' For each replicate, the sum of periods of observation of each animal before
#' its death is calculated (Delignette-Muller et al., 2014).
#' 
#' @aliases reproData print.reproData summary.reproData
#' 
#' @param data A dataframe with five columns:
#' \describe{
#' \item{replicate}{A vector of class \code{integer} or \code{factor} for
#' replicate identification.}
#' \item{conc}{A vector of class \code{numeric} with tested concentrations
#' (positive values).}
#' \item{time}{A vector of class \code{integer} with time points (positive values).
#' The first time must be 0.}
#' \item{Nsurv}{A vector of class \code{integer} with positive values of
#' the number of alive individuals (positive values) at each time point for
#' each concentration and each replicate.}
#' \item{Nrepro}{A vector of class \code{integer} with the number of offspring
#' (positive values) at each time point for each concentration and each replicate.}
#' }
#' @param x An object of class \code{reproData}.
#' @param object An object of class \code{reproData}.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return A dataframe of class \code{reproData}. A list of two objects:
#' \item{raw.data}{The raw dataframe with five columns corresponding to the
#' argument passed in the function.}
#' \item{transformed.data}{A dataframe with nine columns:
#' \describe{
#' \item{ID}{A vector of class \code{factor}. It is a character string who
#' identify the triplet \code{replicat_conc_time}.}
#' \item{replicate}{A vector of class \code{integer} or \code{factor} for replicate
#' identification.}
#' \item{conc}{A vector of class \code{numeric} with tested concentrations
#' (positive values).}
#' \item{time}{A vector of class \code{integer} with time points (positive values).}
#' \item{Ninit}{A vector of class \code{integer} with the number of individuals
#' at the beginning of the bioassay (positive values).}
#' \item{Nsurv}{A vector of class \code{integer} with positive values of the number
#' of alive individuals (positive values) for each concentration and each replicate.}
#' \item{Nrepro}{A vector of class \code{integer} with the number of offspring
#' (positive values) at each time point for each concentration and each replicate.}
#' \item{Nreprocumul}{A vector of class \code{integer} with the cumulative number
#' of offsrping calculated for each time point, concentration and replicat.}
#' \item{Nindtime}{A vector of class \code{numeric} with the number of
#' individual-days for each time point, concentration and replicat.}
#' }} 
#' 
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
#' # (1) Laod reproduction dataset
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
  # check the data class
  if (!is.data.frame(data))
    stop("data.frame expected !")

  # create a column ID as a string: "replicate_conc_time"
  # check if the three columns exist
  ref.names <- c("replicate", "conc", "time")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]

  if (length(missing.names) != 0)
    stop(paste("The column ", missing.names,
               " is missing or have a wrong name.", sep = ""))
  
  # test the integrity of the data with survDataCheck
  if (!is.null(reproDataCheck(data)$id))
    stop("There is one or more error in the data ! Please use the reproDataCheck function before the reproData function !")

  # raw data
  raw.data <- data
  
  # reorder dataset by replicate concentration and time
  data <- data[order(data$replicate, data$conc, data$time), ]
  
  # create an ID column of triplet replicate_conc_time
  data[, "ID"] <- idCreate(data, notime = FALSE)

  # calcul of Ninit and Nreprocumul and Nindtime create transformed.data
  transformed.data <- reproTransformData(data)

  # output of class survData
  out <- list(raw.data = raw.data,
              transformed.data = transformed.data)

  class(out) <- "reproData"

  return(out)
}
