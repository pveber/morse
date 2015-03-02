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
  data <- res2[order(res2$replicate, res2$conc, res2$time), ]  # reorder dataset

  T <- unique(data$time) # times of obs without repetitions
  finalnbr <- match(max(data$time), T) # index of the time at which we want the
  #  results in vector T
  if (finalnbr == 1)
    stop("!!!! It isn't possible to use the first observation time as the last observation time !!!!")
  
  tableTi = list()
  
  for(i in 2:finalnbr) {
    # original dataset at T[i-1]
    dataTim1 <- subset(data, data$time==T[i-1])
    # original dataset at T[i]
    dataTi <- subset(data, data$time==T[i])
    
    # check if data have been properly captured if replicate exists
    if (any(dataTim1$replicate!=dataTi$replicate) || any(dataTim1$conc!=dataTi$conc))
      warning("!!!! BE CAREFUL concentrations and/or replicates are not identical at each time !!!!")
    
    tableTi[[i]] <- data.frame(ID = dataTi$ID,
                               replicate = dataTi$replicate,
                               conc = dataTi$conc,
                               time = dataTi$time,
                               Ninit = dataTi$Ninit,
                               Nsurv = dataTi$Nsurv)
  }
  
  tablefinale <- do.call("rbind", tableTi)
  return(tablefinale)
}