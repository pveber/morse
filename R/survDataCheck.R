#' Checks if an object can be used to perform survival analysis
#'
#' The \code{survDataCheck} function can be used to check if an object
#' containing survival data is formatted according to the expectations of the
#' \code{survData} function.
#'
#'
#' @aliases survDataCheck print.survDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, the function may produce diagnosis plots
# FIXME @param x An object of class survDataCheck.
# FIXME @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns an object of class \code{survDataCheck}, which
#' basically is a dataframe with two columns \code{id} and \code{msg} of
#' character strings. When no error is detected the \code{survDataCheck}
#' function returns an empty dataframe.
#' Here is the list of possible error \code{id}s and their signification:
#' \tabular{rl}{
#' \code{missingColumn} \tab at least one expected column heading is missing \cr
#' \code{firstTime0} \tab the first time point for some (concentration, replicate) is not 0 \cr
#' \code{concNumeric} \tab column \code{conc} contains a value of class other than \code{numeric} \cr
#' \code{NsurvInteger} \tab column \code{Nsurv} contains a value of class other than \code{integer} \cr
#' \code{tablePositive} \tab some data are negative \cr
#' \code{Nsurv0T0} \tab \code{Nsurv} is 0 at time 0 for some (concentration, replicate) \cr
#' \code{duplicateID} \tab there are two identical (\code{replicate}, \code{conc}, \code{time}) triplets \cr
#' \code{missingReplicate} \tab a replicate is missing at some (time point, concentration) \cr
#' \code{NsurvIncrease} \tab \code{Nsurv} increases at some time point of some (concentration, replicate) \cr
#' \code{ReplicateLabel} \tab replicate labels differ between two time points at some concentration \cr
#' }
#'
#' @note If an error of type \code{missingColumn} is detected, the function
#' \code{suvDataCheck} is stopped.
#'
#'
# FIXME @seealso \code{\link{survFullPlot}}, \code{\link{survData}}
#'
#' @keywords check
#'
# FIXME @examples
# # Run the check data function
# data(zinc)
# survDataCheck(zinc)
#
# # Example with an error in the dataframe
# # (1) Load the data
# data(zinc)
#
# # (2) Insert an error (increase the number of survivors at a certain time
# # point compared to its value at the previous time point within the same
# # replicate)
# zinc[25, "Nsurv"] <- 20
# zinc$Nsurv <- as.integer(zinc$Nsurv)
# check <- survDataCheck(zinc, diagnosis.plot = TRUE)
#
# # (3) Check for potential errors in the dataframe
# check
#'
#' @importFrom stringr str_c
#' @export
#'
survDataCheck <- function(data, diagnosis.plot = TRUE) {
  # make a singleton error dataframe (taking care of string/factor conversion
  # issues)
  error <- function(id, msg) {
    data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  }
  # return value: errors will be stored there once found
  errors <- data.frame(stringsAsFactors = FALSE)

  ##
  ## 1. assert column names are correct
  ##
  ref.names <- c("replicate","conc","time","Nsurv")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]
  if (length(missing.names) != 0) {
    errors <- error("missingColumn",
                    paste("The column ", missing.names,
                          " is missing or have a wrong name.", sep = ""))
    class(errors) <- c("survDataCheck", "data.frame")
    return(errors)
  }

  ##
  ## 2. assert the first time point is zero for each (replicate, concentration)
  ##
  subdata <- split(data, list(data$replicate, data$conc), drop = TRUE)
  if (any(unlist(lapply(subdata, function(x) x$time[1] != 0)))) {
    err <- error("firstTime0",
                 "Data are required at time 0 for each concentration and each replicate.")
    errors <- rbind(errors, err)
  }

  ##
  ## 3. assert concentrations are numeric
  ##
  if (!is.numeric(data$conc)) {
    err <- error("concNumeric",
                 "Column 'conc' must contain only numerical values.")
    errors <- rbind(errors, err)
  }

  ##
  ## 4. assert Nsurv contains integer
  ##
  if (!is.integer(data$Nsurv)) {
    err <- error("NsurvInteger",
                 "Column 'Nsurv' must contain only integer values.")
    errors <- rbind(errors, err)
  }

  ##
  ## 5. assert all data are positive
  ##
  table <- subset(data, select = -c(replicate)) # remove replicate column
  if (any(table < 0.0)) {
    err <- error("tablePositive",
                 "Data must contain only positive values.")
    errors <- rbind(errors, err)
  }

  ##
  ## 6. assert Nsurv != 0 at time 0
  ##
  datatime0 <- data[data$time == 0, ]  # select data for initial time points
  if (any(datatime0$Nsurv == 0)) { # test if Nsurv != 0 at time 0
    err <- error("Nsurv0T0",
                 "Nsurv should be different to 0 at time 0 for each concentration and each replicate.")
    errors <- rbind(errors, err)
  }

  ##
  ## 7 assert each (replicate, concentration, time) triplet is unique
  ##
  ID <- idCreate(data, notime = FALSE) # ID vector
  if (any(duplicated(ID))) {
    err <- error("duplicatedID",
                 paste("The triplet Replicate - conc - time: ",
                       ID[duplicated(ID)],
                       " is duplicated.", sep = ""))
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

    ##
    ## 8. assert there is the same number of replicates for each conc and time
    ##
    if (length(subdata$replicate) != length(unique(data$time))) {
      err2 <- error("missingReplicate",
                    paste("Replicate ", unique(subdata$replicate),
                          " is missing for at least one time points at concentration ",
                          unique(subdata$conc), ".", sep = ""))
      consistency.errors <- rbind(consistency.errors, err2)
    }

    ##
    ## 9. assert Nsurv never increases with time
    ##
    nsurv.increase <- subdata$Nsurv[-length(subdata$Nsurv)] < subdata$Nsurv[-1]
    if (any(nsurv.increase)) {
      err2 <- error("NsurvIncrease",
                    paste("For replicate ", unique(subdata$replicate),
                          " and concentration ", unique(subdata$conc),
                          ", Nsurv increases at some time points.",
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

  ##
  ## 10. assert the label of replicate for each time and concentration
  ##
  reslab <- by(data,
               list(data$conc, data$time),
               function(x) {str_c(sort(x$replicate), collapse = "")})

  if (any(reslab != reslab[[1]])) {
    err <- error("ReplicateLabel",
                 "For at least one time and one concentration a replicate label is different from the control.")
    errors <- rbind(errors, err)
  }

  # call function survFullPlot
# FIXME when fn available
#   if (length(err) != 0 && diagnosis.plot && "NsurvIncrease" %in% err) {
#     survFullPlot(data)
#   }
  class(errors) <- c("survDataCheck", "data.frame")
  return(errors)
}
