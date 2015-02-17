#' Check the consistency of a reproduction dataset
#' 
#' The \code{reproDataCheck} function performs several tests on the integrity
#' of the dataset (column headings, type of data\dots) and returns an object
#' of class \code{reproDataCheck}, which is basically a dataframe of error
#' messages. This dataframe is non-empty if the dataset is not in the correct
#' format. The aim of this function is to check the consistency of the
#' dataframe before using function \code{\link{reproData}}. This function
#' highlights possible errors in the data structure that would disturb or
#' prevent the execution of the function \code{\link{reproFitTt}}.
#' 
#' For a given dataframe, the function checks if: \describe{
#' \item{1)}{column headings are correct: \code{replicate} for the column of
#' replicates, \code{conc} for the column of concentrations, \code{time}
#' for the column of time points, \code{Nsurv} for the column of the number of
#' alive individuals and \code{Nrepro} for the column of the number of collected
#' offspring at each time point,}
#' \item{2)}{the first time point of the dataset is 0,}
#' \item{3)}{the class of column \code{conc} is \code{numeric},}
#' \item{4)}{the classes of columns \code{Nsurv} and \code{Nrepro} are
#' \code{integer},}
#' \item{5)}{values of the dataframe are all positive,}
#' \item{6)}{the number of collected offspring is 0 at \eqn{t = 0},}
#' \item{7)}{the number of survivor is not 0 at \eqn{t = 0},}
#' \item{8)}{there is only one triplet \code{replicate} - \code{conc} - \code{time},}
#' \item{9)}{each replicate appears only once per concentration and per time
#' point,}
#' \item{10)}{the number of replicates is the same at any concentration and any
#' time point,}
#' \item{11)}{the number of alive individuals never increases with time,}
#' \item{12)}{at each time \eqn{T}, if the number of alive individuals is null,
#' the number of collected offspring is also null at time \eqn{T+1}.} }
#' 
#' @aliases reproDataCheck print.reproDataCheck
#' 
#' @param data Raw dataframe with five columns. See \code{\link{reproData}}
#' function for details on the required data format.
#' @param diagnosis.plot If \code{TRUE}, calls the default \code{\link{survFullPlot}}
#' function if the number of survivors increases at some time points.
#' @param x An object of class repro.check.data.
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @return The function returns an object of class \code{reproDatacheck}. A
#' dataframe with two columns of character string, \code{id} and \code{msg}.
#' The \code{id} is invisible when displaying the function. Print only shows
#' error messages \code{msg}.
#' \item{id}{The identifier of the test, equals to:
#' \describe{
#' \item{\code{missingColumn}}{if one or more columns are missing or if the column
#' headings are not \code{replicate}, \code{conc}, \code{time}, \code{Nsurv}
#' and \code{Nrepro}.}
#' \item{\code{firstTime0}}{if the first time point is not 0 at each concentration
#' and each replicate.}
#' \item{\code{concNumeric}}{if column \code{conc} does not contain values of
#' class \code{numeric} only.}
#' \item{\code{NsurvInteger}}{if column \code{Nsurv} does not contain values of
#' class \code{integer} only.}
#' \item{\code{NreproInteger}}{if column \code{Nrepro} does not contain values
#' of class \code{integer} only.}
#' \item{\code{tablePositive}}{if there are negative values within the data.}
#' \item{\code{Nrepro0T0}}{if \code{Nrepro} is not 0 at time 0 for each concentration
#' and each replicate.}
#' \item{\code{Nsurv0T0}}{if \code{Nsurv} is 0 at time 0 for one or more
#' concentration and replicate.}
#' \item{\code{duplicateID}}{if there are two or more triplet \code{replicate} -
#' \code{conc} - \code{time}}
#' \item{\code{onlyReplicate}}{if a replicate is duplicated on different lines
#' for the same time points and the same concentration.}
#' \item{\code{missingReplicate}}{if a replicate is missing for at least one time
#' points at one concentration.}
#' \item{\code{NsurvMonotone}}{if \code{Nsurv} increases at some time points
#' compared to the previous one.}
#' \item{\code{Nsurvt0Nreprotp1P}}{if at a giving time \eqn{T}, the number of
#' alive individuals is null and the number of collected offspring is not null
#' for the same replicate and the same concentration at time \eqn{T+1}.}
#' }}
#' \item{msg}{One or more user friendly error messages are generated.}
#' 
#' @note If an error of type \code{missingColumn} is detected, the function
#' \code{reproDataCheck} is stopped. When no error is detected the \code{reproDataCheck}
#' function returns an empty dataframe.
#' 
#' @author Marie Laure Delignette-Muller <marielaure.delignettemuller@@vetagro-sup.fr>,
#' Philippe Veber <philippe.veber@@univ-lyon1.fr>,
#' Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link{survFullPlot}}, \code{\link{reproData}}
#' 
#' @keywords check
#' 
#' @examples
#' 
#' # Run the check data function
#' data(zinc)
#' reproDataCheck(zinc)
#' 
#' # Example with an error in the dataframe
#' 
#' # (1) Load the data
#' data(zinc)
#' 
#' # (2) Insert an error (increase the number of survivors at a certain time
#' # point compared to its value at the previous time point within the same
#' # replicate)
#' zinc[25, "Nsurv"] <- 20
#' zinc$Nsurv <- as.integer(zinc$Nsurv)
#' check <- reproDataCheck(zinc, diagnosis.plot = TRUE)
#' 
#' # (3) Check for potential errors in the dataframe
#' check
#' 
#' @export
#' 
reproDataCheck <-  function(data, diagnosis.plot = TRUE) {
  
  # make a singleton error dataframe (taking care of string/factor conversion
  # issues)
  error <- function(id, msg) {
    data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  }
  
  # return value: errors will be stored there once found
  errors <- data.frame(stringsAsFactors = FALSE)
  
  # 1 run the tests of the survDataCheck
  errors <- survDataCheck(data, diagnosis.plot = FALSE)
  if ("missingColumn" %in% errors$id){
    return(errors)
  }
  
  # 1' test if the column names "Nrepro" exist  
  if (!"Nrepro" %in% colnames(data)) {
    errors <- error("missingColumn",
                    "The column Nrepro is missing or have a wrong name.")
    class(errors) <- c("reproDataCheck", "data.frame")
    return(errors)
  }
  
  # 2' test if Nrepro are integer
  if (!is.integer(data$Nrepro)) {
    err <- error("NreproInteger",
                 "Column 'Nrepro' must contain only integer values.")
    errors <- rbind(errors, err)
  }
  
  # 3'test Nrepro = 0 at time 0
  datatime0 <- data[data$time == 0, ] # select data for initial time points
  if (any(datatime0$Nrepro > 0)) { # test if Nrepro > 0 at time 0
    
    err <- error("Nrepro0T0",
                 "Nrepro should be 0 at time 0 for each concentration and each replicate.")
    errors <- rbind(errors, err)
  }
  
  subdata <- split(data, list(data$replicate, data$conc), drop = TRUE)
  
  .consistency <- function(subdata) {
    # Function to be used on a subdataset corresponding to one replicate at one
    # concentration.
    # This function checks:
    #   - if at each time T for which Nsurv = 0, Nrepro = 0 at time T+1

    # errors consitency dataframe
    errors2 <- data.frame(stringsAsFactors = FALSE)
    
    # 4' test Nsurv = 0 at time t => Nrepro > 0 at time t-1
    NsurvT <- subdata$Nsurv[-length(subdata$Nsurv)]
    NreproTplus1 <- subdata$Nrepro[-1]
    if (any(NreproTplus1[NsurvT == 0] > 0)) {
      err2 <- error("Nsurvt0Nreprotp1P",
                    paste("For replicate ", 
                          unique(subdata$replicate),
                          " and concentration ", unique(subdata$conc),
                          ", there are some Nsurv = 0 followed by Nrepro > 0 at the next time point.",
                          sep = ""))
      errors2 <- rbind(errors2, err2)
    }
    return(errors2)
  }
  
  res <- by(data, list(data$replicate, data$conc), .consistency)
  err <- do.call("rbind", res)
  
  if (length(err)!= 0) {
    errors <- rbind(errors, err)
  }
  
  # call function survFullPlot
  if (length(err)!= 0 && diagnosis.plot && "NsurvMonotone" %in% err) {
      survFullPlot(data)
  }
  
  class(errors) <- c("reproDataCheck", "data.frame")
  return(errors)
}
