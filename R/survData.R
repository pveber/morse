#' Creates a data set for survival analysis
#'
#' This function creates a \code{survData} object from experimental data
#' provided as a \code{data.frame}. The resulting object
#' can then be used for plotting and model fitting. It can also be used
#' to generate \emph{individual-time} estimates.
#'
#' Survival data sets can be under either constant or time-variable exposure profile. The
#' resulting object, in addition to its \code{survData} class, inherits the
#' class \code{survDataCstExp} or \code{survDataVarExp} respectively.
#'
#' The \code{x} argument describes experimental results from a survival
#' toxicity test. Each line of the \code{data.frame}
#' corresponds to one experimental measurement, that is a number of alive
#' individuals at a given concentration at a given time point and in a given replicate.
#'  Note that either the concentration
#' or the number of alive individuals may be missing. The data set is inferred
#' to be under constant exposure if the concentration is constant for each
#' replicate and systematically available. The function \code{survData} fails if
#' \code{x} does not meet the
#' expected requirements. Please run \code{\link{survDataCheck}} to ensure
#' \code{x} is well-formed.
#'
#' @param x a \code{data.frame} containing the following four columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or \code{factor} for replicate
#' identification. A given replicate value should identify the same group of
#' individuals followed in time
#' \item \code{conc}: a vector of class \code{numeric} with tested concentrations
#' (positive values, may contain NAs)
#' \item \code{time}: a vector of class \code{integer} with time points, minimal value must be 0
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
#' @importFrom tibble as_tibble
#' 
#' @examples
#'
#' # (1) Load the survival data set
#' data(zinc)
#'
#' # (2) Create an objet of class 'survData'
#' dat <- survData(zinc)
#' class(dat)
#'
#' @export
#' 
survData <- function(x) {

  x <- as_tibble(x)
  # test the integrity of the data with survDataCheck
  if (dim(survDataCheck(x, diagnosis.plot = FALSE))[1] > 0)
    stop("The [x] argument is not well-formed, please use [survDataCheck] for details.")

  if (!("Ninit" %in% colnames(x))){
    x$Ninit <- Ninit(x)
  }
  
  child_class <-
    if (is_exposure_constant(x)) "survDataCstExp"
    else "survDataVarExp"

  class(x) <- c(child_class, "survData", class(x))

  return(x)
}


#' Test in a well-formed argument to function 'survData' if the concentration
#' is constant and different from \code{NA} for each replicate (each time-serie)
#'
#' @param x an object of class \code{data.frame}
#' @return a boolean \code{TRUE} if concentration in \code{replicate} is constant,
#'  or \code{FALSE} if the concentration in at least one of the replicates is time-variable,
#'  and/or if \code{NA} occures. 
#'
#' @examples
#'
#' # (1) Load the survival data set and test if concentration in replicates is constant
#' data("propiconazole")
#' is_exposure_constant(propiconazole)
#' is_exposure_constant(survData(propiconazole))
#'
#'  # (1) Load the survival data set and test if concentration in replicates is constant
#' data("propiconazole_pulse_exposure") 
#' is_exposure_constant(propiconazole_pulse_exposure)
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


# Computes the effective period of observation in individual days for a
# survival data set
#
# @param x an object of class \code{survData}
# @return a numeric vector
#
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

# Computes a vector associating to each measurement in a \code{survData}
# object the initial number of individuals in the corresponding replicate
#
# @param x an object of class \code{survData}
# @return an integer vector
#
# @importFrom dplyr left_join
#
Ninit <- function(x) {
  x.t0 <- x[x$time == 0 & !is.na(x$Nsurv), c("replicate", "Nsurv")]
  x.t0 <- data.frame(replicate = x.t0$replicate,
                     Ninit = x.t0$Nsurv)
  out <- left_join(x, x.t0, by = c("replicate"))
  return(out$Ninit)
}


#' Joins a concentration with a survival data set into an argument for 'survData'
#' when the concentration varies over time
#'
#' This function joins two data sets, one for exposure measurements, the other
#' for survival measurements, into a single dataframe that can be used
#' with the \code{survData} function.
#'
#' @param x a \code{data.frame} containing the following three columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or \code{factor} for replicate
#' identification
#' \item \code{time}: a vector of class \code{integer} with time points, min value must be 0
#' \item \code{Nsurv}: a vector of class \code{integer} providing the number of
#' alive individuals at some or all time points for each replicate
#' }
#' @param y a \code{data.frame} containing the following three columns:
#' \itemize{
#' \item \code{replicate}: a vector of class \code{integer} or factor for replicate
#' identification
#' \item \code{time}: a vector of class \code{integer} with time points, min value must be 0
#' \item \code{conc}: a vector of class \code{numeric} providing the concentration
#'  at some or all time points for each replicate
#' }
#' 
#' @return a dataframe suitable for `survData`
#'
#' @examples
#'
#' # (1) Load the two survival data sets
#' data(propiconazole_pulse_exposure)
#' exposure <- propiconazole_pulse_exposure[,c("replicate", "time", "conc")]
#' survival <- propiconazole_pulse_exposure[,c("replicate", "time", "Nsurv")]
#'
#' # (2) Create an objet of class 'survData'
#' dat_join <- survData(survData_join(exposure, survival))
#' class(dat_join)
#'
#' @export
#'
survData_join <- function(x, y) {

  ##
  ## 1 assert column names are correct and assign join x or y to corresponding data.frame
  ##
  ref.names_surv <- c("replicate","time","Nsurv")
  ref.names_conc <- c("replicate","time","conc")

  ### data_surv
  if( (all(colnames(x) %in% ref.names_surv) & all(colnames(y) %in% ref.names_conc))
      ||
      (all(colnames(x) %in% ref.names_conc) & all(colnames(y) %in% ref.names_surv))){


    join_df <- dplyr::full_join(x, y,
                                by=c("replicate","time")) %>%
      dplyr::arrange(replicate, time)


  } else {
    stop("'x' and 'y' must both have 3 columns: one with 'replicate','time' and 'Nsurv'
         and the other with 'replicate','time' and 'conc'.")

  }

  return(join_df)
}
