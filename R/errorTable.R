errorTableCreate <- function() {
  t <- data.frame(id = character(0), msg = character(0), stringsAsFactors = FALSE)
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

warningTableCreate <- function() {
  t <- errorTableCreate()
  class(t) <- c("warningTable", class(t))
  t
}

warningTableAppend <- function(...) {
  u <- errorTableAppend(...)
  class(u) <- c("warningTable", class(u))
  u
}

warningTableAdd <- function(t, id, msg) {
  newlines <- data.frame(id = id, msg = msg, stringsAsFactors = FALSE)
  warningTableAppend(t, newlines)
}

errorTableSingleton <- function(id, msg) {
  errorTableAdd(errorTableCreate(), id, msg)
}

errorTableIsEmpty <- function(x)
  dim(x)[1] == 0
  
#' @export
print.errorTable <- function(x, ...) {
  if (errorTableIsEmpty(x)) {
    cat("No error detected.\n")
  }
  else {
    cat("Error(s):\n")
    for (m in x$msg) {
      cat(paste("\t",m,"\n",sep=""))
    }
  }
}


#' @export
print.warningTable <- function(x, ...) {
  if (errorTableIsEmpty(x)) {
    cat("No warning detected.\n")
  }
  else {
    cat("Warning(s):\n")
    for (m in x$msg) {
      cat(paste("\t",m,"\n",sep=""))
    }
  }
}