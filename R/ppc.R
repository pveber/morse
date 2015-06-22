#' Predicted vs obseved plot function
#' 
#' The \code{ppc} function plots the predicted versus observed response values.
#' 
#' @param x An object of class \code{survFitTT} or \code{reproFitTT}
#' 
#' @export
#' @import ggplot2
#' 
ppc <- function(x) {
  # test class object
  if (! is(x, "reproFitTT") && ! is(x, "survFitTT"))
    stop("The object passed in argument [out] is not of class 'reproFitTT' or 'survFitTT' !\n")
  
  dataTT <- x$dataTT
  
  if (is(x, "reproFitTT")) {
    dataTT$obs <- dataTT$Nreprocumul / dataTT$Nindtime
    dataTT$pred <- reproEvalFit(x, dataTT$conc)
    CI <- reproLlmCI(x, dataTT$conc)
    
    
    plot(dataTT$obs, dataTT$obs,
         type = "n",
         bty = "n",
         xaxt = "n",
         yaxt = "n",
         xlim = c(0, max(pretty(c(0, max(CI$qsup95))))),
         ylim = c(0, max(pretty(c(0, max(CI$qsup95))))),
         ylab = "Predicted survival rate",
         xlab = "Observed survival rate")
    abline(a = 0, b = 1)
    
    # axis
    axis(side = 2, at = pretty(c(0, max(CI$qsup95))))
    axis(side = 1, at = pretty(c(0, max(CI$qsup95))))
    
    for (i in 1:length(dataTT$conc)) {
      points(dataTT$obs[i], dataTT$pred[i], pch = 16)
      segments(dataTT$obs[i], CI$qinf95[i], dataTT$obs[i], CI$qsup95[i],
               col = ifelse(CI$qinf95[i] > dataTT$obs[i] && CI$qsup95[i] < dataTT$obs[i], "red", "green"))
    }
  }
  if (is(x, "survFitTT")) {
    dataTT$obs <- dataTT$Nsurv / dataTT$Ninit
    
  }
}