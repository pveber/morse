errorTableCreate <- function() {
  t <- data.frame(stringsAsFactors = FALSE)
  class(t) <- c("errorTable", "data.frame")
  t
}

errorTableAppend <- function(...) {
  u <- rbind(...)
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

#' @export
print.errorTable <- function(x, ...) {
  if (is.null(x$id)) {
    cat("No error detected.\n")
  } else {
    cat("Error(s):\n")
    for (m in x$msg) {
      cat(paste("\t",m,"\n",sep=""))
    }
  }
}
