#' Convergence check of the MCMC chains
#' 
#' The \code{convergence} function checks the convergence of the MCMC chains
#' from the JAGS estimate with the Gelman and Rubin convergence diagnostic
#' (Gelman and Rubin, 1992). It summarizes the \code{mcmc} or \code{mcmc.list}
#' object with a trace of the sampled output, a density estimate and an
#' autocorrelation plot for each parameter in the chain.
#' 
#' 
#' @param out An object of class \code{reproFitTt} or \code{survFitTt}.
#' @param trace If \code{TRUE}, the function traces the sampled output estimate
#' for each parameter in the chain.
#' @param density If \code{TRUE}, the function plots the density estimate for
#' each parameter in the chain.
#' @param autocorr If \code{TRUE}, the function plots the autocorrelation for
#' each parameter in each chain.
#' @param  ppc If \code{TRUE} the function plots the Posterior predictive check.
#' @param style Graphical method: \code{generic} or \code{ggplot}.
#' 
#' @return The function returns an object of class list with the point estimate
#' of the multivariate potential scale reduction factor and the point estimate
#' of the potential scale reduction factor (Rhat) for each parameter of the
#' Gelman and Rubin test (Gelman and Rubin, 1992). A value close to 1 is
#' expected when convergence is reached. See the
#' \code{\link[coda]{gelman.diag}} help for more details.
#' 
#' @note When \code{style = "ggplot"}, the function calls packages \code{ggmcmc}
#' and \code{gridExtra} and returns a graphical object of class \code{ggplot}.
#' 
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>
#' 
#' @seealso \code{\link{reproFitTt}} and \code{\link[coda]{gelman.diag}},
#' \code{\link[coda]{plot.mcmc}}, \code{\link[coda]{autocorr.plot}} from the
#' \code{rjags} package and \code{\link[ggmcmc]{ggs_traceplot}},
#' \code{\link[ggmcmc]{ggs_density}} and
#' \code{\link[ggmcmc]{ggs_autocorrelation}} from the \code{ggmcmc} package
#' (\url{http://xavier-fim.net/packages/ggmcmc})
#' 
#' @references Gelman, A. and Rubin, D.B. (1992) \emph{Inference from iterative
#' simulation using multiple sequences}, Statistical Science, 7, 457-511.
#' 
#' @keywords mcmc-analysis
#' 
#' @examples
#' 
#' # (1) Load the data
#' data(zinc)
#' 
#' # (2) Create an object of class "reproData"
#' dat <- reproData(zinc)
#' 
#' \dontrun{
#' # (3) Run the reproFitTt function
#' out <- reproFitTt(dat)
#' 
#' # (4) Check the convergence
#' convergence(out, trace = TRUE, density = FALSE, 
#' autocorr = TRUE)
#' 
#' # (5) Check the convergence using the "ggmcmc" package
#' convergence(out, trace = TRUE, density = TRUE, 
#' autocorr = TRUE, style = "ggplot")
#' 
#' }
#' 
#' @export
#' 
#' @import ggplot2
#' @importFrom coda autocorr.plot gelman.diag
#' @importFrom ggmcmc ggs ggs_traceplot ggs_density ggs_autocorrelation
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom reshape2 melt
#' 
convergence <- function(out,
                        trace = TRUE,
                        density = TRUE,
                        autocorr = TRUE,
                        ppc = TRUE,
                        style = "generic") {
  # test class object
  if (class(out) != "reproFitTT" && class(out) != "survFitTT")
    stop("The object passed in argument [out] is not of class 'reproFitTT' or 'survFitTT' !\n")
  
  # keep MCMC list
  mcmc <- out$mcmc
  #  Gelamn and Rubin diagnostic
  GelRubmulti <- gelman.diag(mcmc)$mpsrf
  GelRubesti <- gelman.diag(mcmc)$psrf[, "Point est."]
  
  # return the psrf and mprsf value
  cat("Gelman and Rubin:\n")
  cat("Potential scale reduction factor for each parameter:\n")
  print(GelRubesti)
  cat("\nMultivariate potential scale reduction factor:\n", GelRubmulti,"\n")
  
  # PPC
  if (ppc) {
    # create null model
    nodata <- out$jags.data
    nodata["Nsurv"] <- NULL
    # select model text
    if (out$det.part == "loglogisticbinom_2") {
      model.text <- llbinom2.model.text
    }
    if (out$det.part == "loglogisticbinom_3") {
      model.text <- llbinom3.model.text
    }
    # load model text in a temporary file
    model.file <- tempfile() # temporary file address
    fileC <- file(model.file) # open connection
    writeLines(model.text, fileC) # write text in temporary file
    close(fileC) # close connection to temporary file
    # null model
    M0 <- jags.model(model.file, data = nodata, n.chains = out$n.chains,
                     n.adapt = 5000)
    # Sampling
    M0.mcmc <- coda.samples(M0, out$parameters,
                            n.iter = 5000,
                            progress.bar = "none")
    
    M0.mcmc <- do.call("rbind", M0.mcmc)
    mcmcTot <- do.call("rbind", out$mcmc)
    M0.mcmcTot <- data.frame(cbind(M0.mcmc[1:length(mcmcTot[,1]),],
                                   mcmcTot))
    
    # PPC
    # Define data
    concentrations <- out$dataTT$conc
    observation <- out$dataTT$Nsurv / out$dataTT$Ninit
    parameters <- out$parameters
    CI.ppc <- survLlbinomCi(out)
  }
  
  # generic plot
  if (style == "generic") {
    # trace and density
    if (trace || density) {
      plot(mcmc, trace = trace, density = density)
    }
    # autocorrelation
    if (autocorr ) {
      if (trace || density) {
        if (Sys.getenv("RSTUDIO") == "") dev.new() # create a new page plot
        # when not use RStudio
      }
      autocorr.plot(mcmc, ask = TRUE)
    }
    # posterior predictive check
    if (ppc) {
      if (trace || density || autocorr) {
        if (Sys.getenv("RSTUDIO") == "") dev.new() # create a new page plot
        # when not use RStudio
      }
      
      par(mfrow = c(2, 2))
      
      # CPPS
      for (i in 1:length(parameters)) {
        plot(density(M0.mcmcTot[, length(parameters) + i]),
             xlim = range(density(M0.mcmcTot[, i])$x),
             main = "",
             xlab = colnames(M0.mcmcTot)[i],
             lty = 1) # posterior
        lines(density(M0.mcmcTot[, i]),
              lty = 2) # prior
        legend("topleft", legend = c("prior",
                                     "posterior"),
               lty = c(2, 1), bty = "n")
      }
    }
  }
  
  # ggplot
  if (style == "ggplot") {
    # creat ggs objects
    D <- ggs(mcmc)
    if (trace) trp <- ggs_traceplot(D) # trace plot
    if (density) dns <- ggs_density(D) # density plot
    if (autocorr) atc <- ggs_autocorrelation(D) #autocorr plot
    if(! trace || ! density || ! autocorr) {
      blank <- grid.rect(gp = gpar(col = "white"))
    }
    if(trace || density || autocorr) {
      do.call(grid.arrange, list(if (trace) trp else blank,
                                 if (density) dns else blank,
                                 if (autocorr) atc else blank,
                                 ncol = 2))
    }
    
    if (ppc) {
      if (Sys.getenv("RSTUDIO") == "") dev.new() # create a new page plot
      # when not use RStudio
      
      # PPC
      
      plt = list()
      
      mcmc.ppc <- melt(list(prior = M0.mcmc,
                            posterior = do.call("rbind", out$mcmc)))
      mcmc.ppc <- split(mcmc.ppc, mcmc.ppc$Var2)
      
      # CPPS
      for (i in parameters) {
        plt[[i]] <- ggplot(mcmc.ppc[[i]], aes(x = value, color = L1)) +
          geom_density() + labs(x = i) + theme_minimal() +
          theme(legend.title = element_blank())
      }
      
      do.call(grid.arrange, plt)
      
    }
    
  }
  
  return(invisible(list(mpsrf = GelRubmulti,
                        psrf = GelRubesti)))
}
