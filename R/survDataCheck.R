#' Checks if an object can be used to perform survival analysis
#'
#' The \code{survDataCheck} function can be used to check if an object
#' containing survival data is formatted according to the expectations of the
#' \code{survData} function.
#'
#'
#' @aliases survDataCheck
#'
#' @param data any object
#' @param diagnosis.plot if \code{TRUE}, the function may produce diagnosis plots
#'
#' @return The function returns a dataframe of class \code{msgTable} and \code{data.frame} with two columns: \code{id} and \code{msg} of
#' character strings. When no error is detected the object is empty.
#' Here is the list of possible error \code{id}s with their meaning:
#' \tabular{rl}{
#' \code{dataframeExpected} \tab an object of class \code{data.frame} is expected \cr
#' \code{missingColumn} \tab at least one expected column heading is missing \cr
#' \code{firstTime0} \tab the first time point for some (concentration, replicate) couples is not 0 \cr
#' \code{concNumeric} \tab column \code{conc} contains a value of class other than \code{numeric} \cr
#' \code{timeNumeric} \tab column \code{time} contains a value of class other than \code{numeric} \cr
#' \code{NsurvInteger} \tab column \code{Nsurv} contains a value of class other than \code{integer} \cr
#' \code{tablePositive} \tab some data are negative \cr
#' \code{Nsurv0T0} \tab \code{Nsurv} is 0 at time 0 for some (concentration, replicate) \cr
#' \code{duplicateID} \tab there are two identical (\code{replicate}, \code{conc}, \code{time}) triplets \cr
#' \code{NsurvIncrease} \tab \code{Nsurv} increases at some time point of some (concentration, replicate) \cr
#' \code{maxTimeDiffer} \tab maximum time for concentration is lower than maximum time for survival \cr
#' }
#'
#' @note If an error of type \code{dataframeExpected} or \code{missingColumn} is
#' detected, the function \code{survDataCheck} is stopped before looking for
#' other errors.
#'
#' @seealso \code{\link{survData}}
#'
#' @examples
#' # Run the check data function
#' data(zinc)
#' survDataCheck(zinc)
#'
#' # Now we insert an error in the dataset, by artificially increasing the
#' # number of survivors at a given time point, in such a way that the number
#' # of indivuals increases in the corresponding replicate
#' zinc[25, "Nsurv"] <- as.integer(20)
#' survDataCheck(zinc, diagnosis.plot = TRUE)
#'
#' @importFrom stringr str_c
#' @importFrom magrittr '%>%'
#' @importFrom dplyr arrange
#' @importFrom dplyr mutate
#'
#' @export
survDataCheck <- function(data, diagnosis.plot = FALSE) {


  ##
  ## 0. check we have a data.frame
  ##
  if (!("data.frame" %in% class(data))) {
    return(msgTableSingleton("dataframeExpected",
                                "A dataframe is expected."))
  }

  ##
  ## 1. assert column names are correct
  ##
  ref.names <- c("replicate","conc","time","Nsurv")
  missing.names <- ref.names[which(is.na(match(ref.names, names(data))))]
  if (length(missing.names) != 0) {
    msg <- paste("The column ", missing.names,
                 " is missing.", sep = "")
    return(msgTableSingleton("missingColumn",msg))
  }

  # Next errors do not prevent from checking others
  errors <- msgTableCreate()

  ##
  ## 2. assert the first time point is zero for each (replicate, concentration)
  ##
  subdata <- split(data, list(data$replicate), drop = TRUE)
  if (any(unlist(lapply(subdata, function(x) x$time[1] != 0)))) {
    msg <- "Data are required at time 0 for each replicate."
    errors <- msgTableAdd(errors, "firstTime0", msg)
  }

  ##
  ## 3. assert concentrations are numeric
  ##
  if (!is.double(data$conc) && !is.integer(data$conc)) {
    msg <- "Column 'conc' must contain only numerical values."
    errors <- msgTableAdd(errors, "concNumeric", msg)
  }

  ##
  ## 4. assert time is numeric
  ##
  if (!is.numeric(data$time)) {
    msg <- "Column 'time' must contain only numerical values."
    errors <- msgTableAdd(errors, "timeNumeric", msg)
  }

  ##
  ## 5. assert Nsurv contains integer
  ##
  if (!is.integer(data$Nsurv)) {
    msg <- "Column 'Nsurv' must contain only integer values."
    errors <- msgTableAdd(errors, "NsurvInteger", msg)
  }

  ##
  ## 6. assert all data are positive
  ##
  table <- subset(data, select = -c(replicate)) # remove replicate column
  if (any(table < 0.0, na.rm = TRUE)) {
    msg <- "Data must contain only positive values."
    errors <- msgTableAdd(errors, "tablePositive", msg)
  }

  ##
  ## 7. assert Nsurv != 0 at time 0
  ##
  datatime0 <- data[data$time == 0, ]  # select data for initial time points
  if (any(datatime0$Nsurv == 0)) { # test if Nsurv != 0 at time 0
    msg <- "Nsurv should be different to 0 at time 0 for each concentration and each replicate."
    errors <- msgTableAdd(errors, "Nsurv0T0", msg)
  }

  ##
  ## 8. assert each (replicate, time) pair is unique
  ##
  ID <- idCreate(data) # ID vector
  if (any(duplicated(ID))) {
    msg <- paste("The (replicate, time) pair ",
                 ID[duplicated(ID)],
                 " is duplicated.", sep = "")
    errors <- msgTableAdd(errors, "duplicatedID", msg)
  }

  ##
  ## 9. assert all replicates are available for each time point
  ##
  # df_repl <- data %>%
  #   filter(!is.na(Nsurv)) %>%
  #   group_by(time) %>%
  #   summarise(set = paste(sort(replicate), collapse = '-'))
  # 
  # if (length(unique(df_repl$set)) > 1) {
  #   sets <- df_repl %>%
  #     group_by(set) %>%
  #     summarise(cardinal = length(set))
  #   reference_set <- sets$set[which.max(sets$cardinal)]
  #   diffs <- which(df_repl$set != reference_set)
  #   msg <- paste(
  #       "Changing set of replicates for time point(s) ",
  #       paste(diffs, collapse = ", "),
  #       ".",
  #       sep = "")
  #   errors <- msgTableAdd(errors, "missingReplicate", msg)
  # }

  ##
  ## 10. assert Nsurv never increases with time
  ##
  df_variation <- data %>%
    filter(!is.na(Nsurv)) %>%
    group_by(replicate) %>%
    arrange(time) %>%
    mutate(Nprec = ifelse(time == min(time), Nsurv, lag(Nsurv))) %>%
    mutate(decrease = Nsurv <= Nprec) %>%
    summarise(decreasing = all(decrease))

  if (! all(df_variation$decreasing)) {
    replicates <- df_variation$replicate[! df_variation$decreasing]
    msg <- paste(
        "'Nsurv' increases at some time points in replicate(s) ",
        paste(replicates, collapse=", "),
        ".",
        sep = "")
    errors <- msgTableAdd(errors, "NsurvIncrease", msg)
  }

  ##
  ## 11. Assert max(time in data_conc) >= max(time in data_surv)
  ##
  df_checkMaxTimeSurv <- data %>%
    filter(!is.na(Nsurv)) %>%
    group_by(replicate) %>%
    filter(time == max(time))
  
  df_checkMaxTimeConc <- data %>%
    filter(!is.na(conc)) %>%
    group_by(replicate) %>%
    filter(time == max(time))
  
  if(!all(df_checkMaxTimeConc$time >= df_checkMaxTimeSurv$time) ){
    msg <- "In each 'replicate', maximum time for concentration record should
    be greater or equal to maximum time in survival data observation."
    errors <- msgTableAdd(errors, "maxTimeDiffer", msg)
  }

  if (diagnosis.plot && "NsurvIncrease" %in% errors$id) {
    survDataPlotFull(data, ylab = "Number of surviving individuals")
  }
  return(errors)
}
