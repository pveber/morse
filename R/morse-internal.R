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