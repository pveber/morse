#' @importFrom stringr str_c
idCreate <- function(data, notime = FALSE) {
  # INPUT
  # data: a raw data frame
  # notime: if TRUE there is no column time in the dataset
  # OUTPUT
  # vector of triplet replicate_conc_time or couple replicate_conc
  if (!notime) {
    return(str_c(data[, "replicate"], data[, "conc"], data[, "time"],
                 sep = "_"))
  } else {
    return(str_c(data[, "replicate"], data[, "conc"], sep = "_"))
  }
}

#' @importFrom dplyr right_join %>% rename
survTransformData <- function(data) {
  # create Ninit 
  # INPUTS
  # - data: a dataframe of class survData with:
  #   - replicate: replicate indentification
  #   - conc: tested concentration
  #   - time: observation time
  #   - Nsurv: number of alive individuals at concentration "conc" and at time
  #     "time"
  #  !!!!!! data are supposed to be sorted by replicate, conc and time !!!!!!
  # OUTPUT: a dataframe with 6 columns
  #   - ID: ID replicate - conc - time
  #   - replicate: replicate indentification
  #   - conc: tested concentration
  #   - time: observation time
  #   - Ninit: number of individuals at the beginning of the bioassay
  #   - Nsurv: number of alive individuals at the time specified in the input of
  #     the function

  # split dataset by time to calculate Ninit
  temp <- split(data, data$time)

  # add column Ninit
  tabletime0 <- data[data$time == 0, ] # control dataset
  tabletime0[, "Ninit"] <- tabletime0[, "Nsurv"]

  survCalculNinit <- function(t1, t2) {
    # calcul the correct number Ninit
    # for each replicate con and time
    # INPUTS
    # - t1: list of splited dataframe
    # - t2: tabletime0
    # OUTPUTS
    # - list of splited dataframe with the news column Ninit

    . = NULL
    ID.x = NULL
    time.x = NULL
    Nsurv.x = NULL

    right_join(t1, t2,
               by = c("replicate","conc"))[, c("ID.x", "replicate", "conc", "time.x",
               "Nsurv.x", "Ninit")] %>% rename(., time = time.x) %>% rename(.,
               Nsurv = Nsurv.x) %>% rename(., ID = ID.x)
  }
  
  res <- lapply(temp, function(x) survCalculNinit(x, tabletime0)) # Ninit
  res2 <- do.call(rbind, res) # return a dataframe
  rownames(res2) <- 1:dim(res2)[1]
  res2$time <- as.numeric(res2$time) # change type of time

  return(res2)
}
