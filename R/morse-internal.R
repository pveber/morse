# Ugly hack to get rid of spurious notes in package check, caused by uses
# of dplyr::{rename, filter}. R is such a sad language.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("Nsurv","conc","Ninit"))


# Generates a character string vector from a data.frame using its replicate,
# conc and time columns. The result can be used as identifiers for the rows
# of the data.set.
#
#' @importFrom stringr str_c
#'
idCreate <- function(data) {
  str_c(data[, "replicate"],
        data[, "conc"],
        data[, "time"],
        sep = "_")
}
