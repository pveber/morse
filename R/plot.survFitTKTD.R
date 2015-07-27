Surv <- function (Cw, time , ks, ke, NEC, m0)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- exp(-m0*time) # survie de base avec mortalite naturelle seule
  if(Cw > NEC) {
    tNEC <- -(1/ke)*log(1 - NEC/Cw)
    if (time > tNEC) {
      # ajoute de la mortalite due au toxique
      S <- S * exp( ks/ke*Cw*(exp(-ke*tNEC) -exp(-ke*time)) - ks*(Cw-NEC)*(time - tNEC) )
    }
  }
  return(S)
}

survFitPlotDataTKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  dtheo <- data.frame(conc = numeric(), t = numeric(), psurv = numeric())
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # parameters
  ks <- x$estim.par["ks", "median"]
  ke <- x$estim.par["ke", "median"]
  nec <- x$estim.par["nec", "median"]
  m0 <- x$estim.par["m0", "median"]
  
  for (i in 1:length(concobs)) {
    for (j in 1:npoints) {
      psurv <- Surv(Cw = concobs[i], time = tfin[j],
                    ks = ks, ke = ke,
                    NEC = nec,
                    m0 = m0)
      dtheo <- rbind(dtheo, data.frame(conc = concobs[i],
                                       t = tfin[j],
                                       psurv = psurv))
    }
  }
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     t = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init)
  
  return(list(dtheo = dtheo,
              dobs = dobs))
}

#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTKTD <- function(x,
                             xlab,
                             ylab,
                             main,
                             style = "generic") {
  
  # default axis parameters
  if (missing(xlab)) {
    xlab <- "Time"
  }
  if (missing(ylab)) {
    ylab <- "Response"
  }
  if (missing(main)) {
    main <- NULL
  }
  
  # create observed and theoretical values
  data <- survFitPlotDataTKTD(x)
  
  if (style == "generic") {
    # vector color
    data[["dobs"]]$color <- as.numeric(as.factor(data[["dobs"]][["conc"]]))
    data[["dtheo"]]$color <- as.numeric(as.factor(data[["dtheo"]][["conc"]]))
    plot(data[["dobs"]][["t"]],
         data[["dobs"]][["psurv"]],
         xlab = xlab,
         ylab = ylab,
         pch = 16,
         col = data[["dobs"]]$color,
         main = main)
    
    # one line by replicate
    by(data[["dtheo"]], list(data[["dtheo"]]$conc),
       function(x) {
         lines(x$t, x$psurv, # lines
               col = x$color)
       })
  }
  
  if (style == "ggplot") {
    
    plt1 <- ggplot(data$dobs, aes(x = t, y = psurv, colour = conc,
                                  group = as.factor(conc))) +
      labs(x = xlab, y = ylab) + ggtitle(main) +
      geom_point() + geom_line(data = data$dtheo) + theme_minimal()
    
    plt1
  }
}
