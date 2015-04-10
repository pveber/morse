#' Plot of survival data
#' 
#' This function plots the survival rate in function for a fixed concentration
#' and a fixed timepoint.
#' 
#' @param data an object of class \code{survData}.
#' @param target.time a numeric value corresponding to some observed time in
#' \code{data}.
#' @param concentration a numeric value corresponding to some observed concentration
#' in \code{data}.
#' @param xlab X-axis label.
#' @param ylab Y-axis label.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' @param addlegend If \code{TRUE}, a default legend is added to the plot.
#' 
#' @export
#' @import ggplot2
#' @importFrom dplyr %>% filter
survDataPlotReplicates <- function(data,
                                   target.time = NULL,
                                   concentration = NULL,
                                   xlab = NULL,
                                   ylab = NULL,
                                   style = "generic",
                                   addlegend = TRUE) {
  
  # response variable
  data$response <- data$Nsurv / data$Ninit
  
  # default argument
  if (is.null(xlab)) {
    xlab <- "Replicates"
  }
  if (is.null(ylab)) {
    ylab <- "Survival rate"
  }
  if (is.null(concentration)) {
    concentration <- max(data$conc)
  }
  if (is.null(target.time)) {
    target.time <- max(data$time)
  }
  
  # check [target.time] and [concentration]
  if (!target.time %in% data$time)
    stop("[target.time] is not one of the possible time !")
  
  if (!concentration %in% data$conc)
    stop("The choosen value for [concentration] is not possible !")
  
  # select the concnetration and the target.time
  data <- filter(data, conc == concentration & time == target.time)
  
  if (style == "generic") {
    plot(factor(data$replicate), data$response,
         type = "n",
         xlab = xlab,
         ylab = ylab)
  }
  
  if (style == "ggplot") {
    df <- ggplot(data, aes(x = replicate, y = response))
    df + geom_point() + labs(x = xlab, y = ylab)
  }
}
