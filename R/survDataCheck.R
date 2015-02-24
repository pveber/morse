#' Check the consistency of a survival dataset
#' 
#' The \code{survDataCheck} function performs several tests on the integrity of
#' the dataset (column headings, type of data\dots) and returns an object of
#' class \code{survDataCheck}, which is basically a dataframe of error
#' messages. This dataframe is non-empty if the dataset is not in the correct
#' format. The aim of this function is to check the consistency of the dataframe
#' before using function \code{\link{survData}}. This function highlights
#' possible errors in the data structure that would disturb or prevent the
#' execution of the function \code{\link{survFitTt}}.
#' 
#' For a given dataframe, the function checks if:
#' \describe{
#' \item{1)}{column headings are correct: \code{replicate} for the column of
#' replicates, \code{conc} for the column of concentrations, \code{time}
#' for the column of time points and \code{Nsurv} for the column of the number
#' of alive individuals,}
#' \item{2)}{the first time point of the dataset is 0,}
#' \item{3)}{the class of column \code{conc} is \code{numeric},}
#' \item{4)}{the classes of columns \code{Nsurv} is \code{integer},}
#' \item{5)}{values of the dataframe are all positive,}
#' \item{6)}{the number of survivor is not 0 at \eqn{t = 0},}
#' \item{7)}{there is only one triplet \code{replicate} - \code{conc} - \code{time},}
#' \item{8)}{each replicate appears only once per concentration and per time point,}
#' \item{9)}{the number of replicates is the same at any concentration and any
#' time point,}
#' \item{10)}{the number of alive individuals never increases with time,}
#' }
#' 
#' @aliases survDataCheck print.survDataCheck
#' 
#' @param data Raw dataframe with four columns. See \code{\link{survData}}
#' function for details on the required data format.
#' @param diagnosis.plot If \code{TRUE}, calls the default \code{\link{survFullPlot}}
#' function if the number of survivors increases at some time points.
#' @param x An object of class survDataCheck.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return The function returns an object of class \code{survDataCheck}. A
#' dataframe with two columns of character string, \code{id} and \code{msg}.
#' The \code{id} is invisible when displaying the function. Print only shows
#' error messages \code{msg}.
#' \item{id}{The identifier of the test, equals to:
#' \describe{
#' \item{\code{missingColumn}}{if one or more columns are missing or if the column
#' headings are not \code{replicate}, \code{conc},\code{time} and \code{Nsurv}.}
#' \item{\code{firstTime0}}{if the first time point is not 0 at each concentration
#' and each replicate.}
#' \item{\code{concNumeric}}{if column \code{conc} does not contain values of
#' class \code{numeric} only.}
#' \item{\code{NsurvInteger}}{if column \code{Nsurv} does not contain values of
#' class \code{integer} only.}
#' \item{\code{tablePositive}}{if there are negative values within the data.}
#' \item{\code{Nsurv0T0}}{if \code{Nsurv} is 0 at time 0 for one or more
#' concentration and replicate.}
#' \item{\code{duplicateID}}{if there are two or more triplet \code{replicate} -
#' \code{conc} - \code{time}}
#' \item{\code{uniqueReplicateNumberPerCondition}}{if a replicate is duplicated
#' on different lines for the same time points and the same concentration.}
#' \item{\code{missingReplicate}}{if a replicate is missing for at least one time
#' points at one concentration.}
#' \item{\code{NsurvMonotone}}{if \code{Nsurv} increases at some time points
#' compared to the previous one.}
#' }}
#' \item{msg}{One or more user friendly error messages are generated.}
#' 
#' @note If an error of type \code{missingColumn} is detected, the function
#' \code{suvDataCheck} is stopped. When no error is detected the \code{survDataCheck}
#' function returns an empty dataframe.
#' 
#' @author Marie Laure Delignette-Muller <marielaure.delignettemuller@@vetagro-sup.fr>,
#' Philippe Veber <philippe.veber@@univ-lyon1.fr>,
#' Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link{survFullPlot}}, \code{\link{survData}}
#' 
#' @keywords check
#' 
#' @examples
#' # Run the check data function
#' data(zinc)
#' survDataCheck(zinc)
#' 
#' # Example with an error in the dataframe
#' # (1) Load the data
#' data(zinc)
#' 
#' # (2) Insert an error (increase the number of survivors at a certain time
#' # point compared to its value at the previous time point within the same
#' # replicate)
#' zinc[25, "Nsurv"] <- 20
#' zinc$Nsurv <- as.integer(zinc$Nsurv)
#' check <- survDataCheck(zinc, diagnosis.plot = TRUE)
#' 
#' # (3) Check for potential errors in the dataframe
#' check
#' 
#' @export
#' @importFrom stringr str_c
#' 
survDataCheck <- function(data, diagnosis.plot = TRUE) {
  # make a singleton error dataframe (taking care of string/factor conversion
  # issues)
  error <- function(id, msg) {
    data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  }
  # return value: errors will be stored there once found
  errors <- data.frame(stringsAsFactors = FALSE)
  # 1 test if the column names are correct
  ref.names <- c("replicate","conc","time","Nsurv")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]
  if (length(missing.names) != 0) {
    errors <- error("missingColumn",
                    paste("The column ", missing.names,
                          " is missing or have a wrong name.", sep = ""))
    class(errors) <- c("survDataCheck", "data.frame")
    return(errors)
  }
  # 2 test if the first time point is zero
  subdata <- split(data, list(data$replicate, data$conc), drop = TRUE)
  if (any(unlist(lapply(subdata, function(x) x$time[1] != 0)))) {
    err <- error("firstTime0",
                 "Data are required at time 0 for each concentration and each replicate.")
    errors <- rbind(errors, err)
  }
  # 3 test if concentrations are numeric
  if (!is.numeric(data$conc)) {
    err <- error("concNumeric",
                 "Column 'conc' must contain only numerical values.")
    errors <- rbind(errors, err)
  }
  # 4 test if Nsurv are integer
  if (!is.integer(data$Nsurv)) {
    err <- error("NsurvInteger",
                 "Column 'Nsurv' must contain only integer values.")
    errors <- rbind(errors, err)
  }
  # 5 positivity test table
  table <- subset(data, select = -c(replicate)) # remove replicate column
  if (any(table < 0.0)) {
    err <- error("tablePositive",
                 "Data must contain only positive values.")
    errors <- rbind(errors, err)
  }
  # 6 test Nsurv != 0 at time 0
  datatime0 <- data[data$time == 0, ]  # select data for initial time points
  if (any(datatime0$Nsurv == 0)) { # test if Nsurv != 0 at time 0
    err <- error("Nsurv0T0",
                 "Nsurv should be different to 0 at time 0 for each concentration and each replicate.")
    errors <- rbind(errors, err)
  }
  # 7 test unique triplet Replicat - conc - time only if survDataCheck is called
  # by survData or reproData
  # FIXME: replace ID creation by function call
  ID <- str_c(data[, "replicate"], data[, "conc"], data[, "time"],
              sep = "_")
  if (any(duplicated(ID))) {
    err <- error("duplicatedID",
                 paste("The triplet Replicate - conc - time: ",
                       ID[duplicated(ID)], " is duplicated.",
                       sep = ""))
    errors <- rbind(errors, err)
  }
  consistency <- function(subdata) {
    # Function to be used on a subdataset corresponding to one replicate at one
    # concentration.
    # This function checks:
    #   - if each replicate appears once and only once at each time
    #   - if Nsurv is never increasing with time
    
    # errors consistency dataframe
    consistency.errors <- data.frame(stringsAsFactors = FALSE) 
    # 8 test if each replicate appears once and only once at each conc and time
    if (nrow(subdata) != length(unique(subdata$time))) {
      err2 <- error("uniqueReplicateNumberPerCondition",
                    paste("Replicate ",
                          unique(subdata$replicate),
                          " appears on different lines for the same time point and the same concentration ", 
                          unique(subdata$conc), ".", sep = ""))
      consistency.errors <- rbind(consistency.errors, err2)
    }
    # 9 test if there is no lack of replicate at each conc and time
    if (length(subdata$replicate) != length(unique(data$time))) {
      err2 <- error("missingReplicate",
                    paste("Replicate ", unique(subdata$replicate),
                          " is missing for at least one time points at concentration ", 
                          unique(subdata$conc), ".", sep = ""))
      consistency.errors <- rbind(consistency.errors, err2)
    }
    # 10 test Nsurv monotony
    nonmonotonous <- subdata$Nsurv[-length(subdata$Nsurv)] < subdata$Nsurv[-1]
    if (any(nonmonotonous)) {
      err2 <- error("NsurvMonotone",
                    paste("For replicate ", unique(subdata$replicate),
                          " and concentration ", unique(subdata$conc),
                          ", Nsurv increases at some time points compared to the previous one.",
                          sep = ""))
      consistency.errors <- rbind(consistency.errors, err2)
    }
    return(consistency.errors)
  }
  res <- by(data, list(data$replicate, data$conc), consistency)
  err <- do.call("rbind", res)
  if (length(err) != 0) {
    errors <- rbind(errors, err)
  }
  # call function survFullPlot
  if (length(err) != 0 && diagnosis.plot && "NsurvMonotone" %in% err) {
    survFullPlot(data)
  }
  class(errors) <- c("survDataCheck", "data.frame")
  return(errors)
}
