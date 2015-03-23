errorTableCreate <- function() {
  t <- data.frame(stringsAsFactors = FALSE)
  class(t) <- c("errorTable", "data.frame")
  t
}

errorTableAppend <- function(s, t) {
  u <- rbind(s,t)
  class(u) <- c("errorTable", "data.frame")
  u
}

errorTableAdd <- function(t, id, msg) {
  newlines <- data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  errorTableAppend(t, newlines)
}

errorTableSingleton <- function(id, msg) {
  errorTableAdd(errorTableCreate(), id, msg)
}

