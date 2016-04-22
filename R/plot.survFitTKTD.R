#' Plotting method for \code{survFitTKTD} objects
#' 
#' This is the generic \code{plot} S3 method for the \code{survFitTKTD}.
#' It plots time-exposure-response fits for each concentration of
#' survival analysis.
#' 
#' The fitted curves represent the \strong{estimated survival rate} as a function
#' of time for each concentration;
#' When \code{adddata = TRUE} the black dots depict the \strong{observed survival
#' rate} at each time point. Note that since our model does not take
#' inter-replicate variability into consideration, replicates are systematically
#' pooled in this plot.
#' The function plots both 95 \% credible intervals for the estimated survival
#' rate (by default the red area around the fitted curve) and 95 \% confidence
#' intervals for the observed survival rate (as black error bars if
#' \code{adddata = TRUE}).
#' Both types of intervals are taken at the same level. Typically
#' a good fit is expected to display a large overlap between the two intervals.
#' If \code{spaghetti = TRUE}, the credible intervals are represented by two
#' dotted ines limiting the credible band, and a spaghetti plot is added to this
#' band.
#' It consists of the representation of simulated curves using parameter values
#' sampled in the posterior distribution (2 \% of the MCMC chains are randomly
#' taken for this sample).
#' 
#' @param x An object of class \code{survFitTKTD}.
#' @param xlab A label for the \eqn{X}-axis, by default \code{Time}.
#' @param ylab A label for the \eqn{Y}-axis, by default \code{Survival rate}.
#' @param main A main title for the plot.
#' @param spaghetti if \code{TRUE}, the credible interval is represented by 
#' multiple curves
#' @param one.plot if \code{TRUE}, draws all the estimeted curves in one plot.
#' @param adddata if \code{TRUE}, adds the observed data with confidence interval
#' to the plot
#' @param addlegend if \code{TRUE}, adds a default legend to the plot.
#' @param style graphical backend, can be \code{'generic'} or \code{'ggplot'}
#' @param \dots Further arguments to be passed to generic methods.
#' 
#' @keywords plot 
#' 
#' @examples
#' 
#' # (1) Load the survival data
#' data(propiconazole)
#' 
#' # (2) Create an object of class "survData"
#' dat <- survData(propiconazole)
#' 
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#'
#' # (4) Plot the fitted curves in one plot 
#' plot(out)
#'
#' # (5) Plot one fitted curve by concentration with credible limits as
#' # spaghetti, data and confidence intervals
#' # and with a ggplot style
#' plot(out, spaghetti = TRUE , adddata = TRUE, one.plot = FALSE,
#'      style = "ggplot")
#' }
#' 
#' @export
#' 
#' @import ggplot2
#' @import grDevices
#' @importFrom reshape2 melt
#' @importFrom gridExtra grid.arrange arrangeGrob
#' @importFrom grid grid.rect gpar
#' @importFrom graphics plot
#' 
plot.survFitTKTD <- function(x,
                             xlab = "Time",
                             ylab = "Survival rate",
                             main = NULL,
                             spaghetti = FALSE,
                             one.plot = TRUE,
                             adddata = FALSE,
                             addlegend = FALSE,
                             style = "generic", ...) {
  
  if (one.plot && !is.null(concentration))
    one.plot <- FALSE
  
  if ((addlegend && is.null(concentration)) ||
      (addlegend && !one.plot))
    warning("The legend is available only if [one.plot] is TRUE or if [concentration] is not NULL !", call. = FALSE)

  if (!is.null(concentration) && !any(x$transformed.data$conc == concentration))
    stop("The [concentration] argument is not one of the possible concentration !")
  
  if (one.plot)
    warning("The credible limits and confidence intervals are not drawn when 'one.plot' = TRUE.",
            call. = FALSE)

  conf.int <- survTKTDConfInt(x)
  
  data.credInt <- survFitPlotCITKTD(x)
  
  data.credInt[["dobs"]] <- data.frame(data.credInt[["dobs"]],
                                       qinf95 = conf.int["qinf95", ],
                                       qsup95 = conf.int["qsup95",],
                                       Conf.Int = "Confidence interval",
                                       Points = "Observed values")
  
  data.credInt[["dtheoQ"]]$Cred.Lim <- "Credible limits"
  
  # vector color
  data.credInt[["dobs"]]$color <-
    as.numeric(as.factor(data.credInt[["dobs"]][["conc"]]))
  data.credInt[["dtheoQ"]]$color <-
    as.numeric(as.factor(data.credInt[["dtheoQ"]][["conc"]]))
  
  dataCIm <- melt(data.credInt[["dtheoSp"]],
                  id.vars = c("conc", "time"))
  
  if (style == "generic") {
    survFitPlotTKTDGeneric(data.credInt, xlab, ylab, main, one.plot, spaghetti,
                           dataCIm, adddata, addlegend)
  }
  else if (style == "ggplot") {
    survFitPlotTKTDGG(data.credInt, xlab, ylab, main, one.plot, spaghetti,
                      dataCIm, adddata, addlegend)
  }
  else stop("Unknown style")
}

Surv <- function(Cw, time, ks, kd, NEC, m0)
  # Fonction S ecrite en R pour la validation en simu ensuite
  # Cw est la concentration dans le milieu
{
  S <- exp(-m0*time) # survie de base avec mortalite naturelle seule
  if (Cw > NEC) {
    tNEC <- -(1/kd)*log(1 - NEC/Cw)
    if (time > tNEC) {
      # ajoute de la mortalite due au toxique
      S <- S * exp( ks/kd*Cw*(exp(-kd*tNEC) -exp(-kd*time))
                    - ks*(Cw-NEC)*(time - tNEC) )
    }
  }
  return(S)
}

#' @importFrom stats aggregate binom.test
survTKTDConfInt <- function(x) {
  # create confidente interval on observed data
  # binomial model by a binomial test
  # INPUT:
  # - x : object of class survFitTT
  # OUTPUT:
  # - ci : confidente interval
  
  ci <- apply(x$transformed.data, 1, function(x) {
    binom.test(x["N_alive"], x["N_init"])$conf.int
  })
  rownames(ci) <- c("qinf95", "qsup95")
  colnames(ci) <- x$transformed.data$conc
  
  return(ci)
}

survFitPlotCITKTD <- function(x) {
  # INPUT
  # x : An object of class survFitTKTD
  # OUTPUT
  # A list of - dobs : observed values
  #           - dtheo : estimated values
  npoints <- 100
  
  concobs <- unique(x$transformed.data$conc)
  tfin <- seq(0, max(x$jags.data$t), length.out = npoints)
  
  # prameters
  mctot <- do.call("rbind", x$mcmc)
  sel <- sample(nrow(mctot))[1:ceiling(nrow(mctot) / 50)]
  ks <- 10^mctot[, "log10ks"][sel]
  kd <- 10^mctot[, "log10kd"][sel]
  m0 <- 10^mctot[, "log10m0"][sel]
  nec <- 10^mctot[, "log10NEC"][sel]
  
  # all theorical
  dtheo = list()
  for (k in 1:length(concobs)) {
    dtheo[[k]] <- array(data = NA, dim = c(npoints, length(nec)))
    for (i in 1:length(nec)) {
      for (j in 1:npoints) {
        dtheo[[k]][j, i] <- Surv(Cw = concobs[k], time = tfin[j],
                                 ks = ks[i], kd = kd[i],
                                 NEC = nec[i],
                                 m0 = m0[i])
      }
    }
  }
  
  dtheoSp <- do.call("rbind", dtheo)
  dtheoSp <- as.data.frame(cbind(rep(concobs, rep(npoints, length(concobs))),
                                rep(tfin, length(concobs)),
                                dtheoSp))
  names(dtheoSp) <- c("conc", "time", paste0("X", 1:length(sel)))
  
  # quantile
  qinf95 = NULL
  qsup95 = NULL
  q50 = NULL
  
  for (i in 1:dim(dtheoSp)[1]) {
    qinf95[i] <- quantile(dtheoSp[i, 3:length(dtheoSp)],
                          probs = 0.025, na.rm = TRUE)
    qsup95[i] <- quantile(dtheoSp[i, 3:length(dtheoSp)],
                          probs = 0.975, na.rm = TRUE)
    q50[i] <- quantile(dtheoSp[i, 3:length(dtheoSp)],
                       probs = 0.5, na.rm = TRUE)
  }
  
  dtheoQ <- data.frame(conc = dtheoSp[, "conc"], time = dtheoSp[, "time"],
                       qinf95 = qinf95, qsup95 = qsup95, q50 = q50)
  
  dobs <- data.frame(conc = x$transformed.data$conc,
                     time = x$transformed.data$time, 
                     psurv = x$transformed.data$N_alive / x$transformed.data$N_init)
  
  return(list(dtheoQ = dtheoQ,
              dtheoSp = dtheoSp,
              dobs = dobs))
}

survFitPlotTKTDGeneric <- function(data, xlab, ylab, main, one.plot,
                                   spaghetti, dataCIm, adddata,
                                   addlegend) {

    if (one.plot) {
        survFitPlotTKTDGenericOnePlot(data, xlab, ylab, main, adddata, addlegend)
    } else if (!one.plot && is.null(concentration)) {
        par(mfrow = plotMatrixGeometry(length(unique(data[["dobs"]][["conc"]]))))

        survFitPlotTKTDGenericNoOnePlot(data, xlab, ylab, spaghetti,
                                        dataCIm, adddata, concentration)
        
        par(mfrow = c(1, 1))
    } else {
        survFitPlotTKTDGenericNoOnePlot(data, xlab, ylab, spaghetti,
                                        dataCIm, adddata, concentration, addlegend)
    }
}

survFitPlotTKTDGenericOnePlot <- function(data, xlab, ylab, main, adddata,
                                          addlegend) {
  plot(data[["dobs"]][["time"]],
       data[["dobs"]][["psurv"]],
       xlab = xlab,
       ylab = ylab,
       type = "n",
       main = main)
  
  # one line by replicate
  by(data[["dtheoQ"]], list(data[["dtheoQ"]]$conc),
     function(x) {
       lines(x$time, x$q50, # lines
             col = x$color)
     })
  
  # points
  if (adddata) {
    points(data[["dobs"]][["time"]],
           data[["dobs"]][["psurv"]],
           pch = 20,
           col = "black")
  }
  if (addlegend) {
    legend("bottomleft",
           legend = unique(data[["dobs"]]$conc),
           pch = ifelse(adddata, 20, NA),
           lty = 1,
           bty = "n",
           cex = 1,
           ncol = 2,
           col = "black",
           title = "Concentrations")
  }
}

survFitPlotTKTDGenericNoOnePlot <- function(data, xlab, ylab, spaghetti,
                                            dataCIm, adddata, concentration,
                                            addlegend = FALSE) {
  
  if (is.null(concentration)) {
    dobs <- split(data[["dobs"]], data[["dobs"]]$conc)
    dtheoQ <- split(data[["dtheoQ"]], data[["dtheoQ"]]$conc)
    dataCIm <- split(dataCIm, dataCIm$conc)
  } else {
    dobs <- list(filter(data[["dobs"]], conc == concentration))
    dtheoQ <- list(filter(data[["dtheoQ"]], conc == concentration))
    dataCIm <- list(filter(dataCIm, conc == concentration))
  }
  
  delta <- 0.01 * (max(data[["dobs"]]$time) - min(data[["dobs"]]$time))
  
  mapply(function(x, y, z) {
    plot(x[, "time"],
         x[, "q50"],
         xlab = xlab,
         ylab = ylab,
         type = "n",
         ylim = c(0, 1),
         main = paste0("Concentration = ", unique(x[, "conc"])),
         col = x[, "color"])
    
    if (spaghetti) {
      color <- "gray"
      color_transparent <- adjustcolor(color, alpha.f = 0.05)
      by(z, z$variable, function(x) {
        lines(x[, "time"], x[, "value"], col = color_transparent)
      })
    } else {
      polygon(c(x[, "time"], rev(x[, "time"])), c(x[, "qinf95"],
                                                  rev(x[, "qsup95"])),
              col = "pink", border = NA)
    }
    
    lines(x[, "time"], x[, "q50"], # lines
          col = x[, "color"])
    lines(x[, "time"], x[, "qinf95"],
          col = x[, "color"], lty = 2)
    lines(x[, "time"], x[, "qsup95"], 
          col = x[, "color"], lty = 2)
    
    if (adddata) {
      points(y[, "time"],
             y[, "psurv"],
             pch = 20,
             col = "black") # points
      segments(y[, "time"], y[, "qinf95"],
               y[, "time"], y[, "qsup95"],
               col = "gray")
      segments(y[, "time"] - delta, y[, "qinf95"],
               y[, "time"] + delta, y[, "qinf95"],
               col = "gray")
      segments(y[, "time"] - delta, y[, "qsup95"],
               y[, "time"] + delta, y[, "qsup95"],
               col = "gray")
    }
    
    if (addlegend) {
      legend("bottomleft", pch = c(ifelse(adddata, 20, NA), NA, NA, NA),
             lty = c(NA, ifelse(adddata, 1, NA), 1, 2),
             col = c(ifelse(adddata, "black", NA),
                     ifelse(adddata, "gray", NA),
                     x[, "color"], x[, "color"]),
             legend = c(ifelse(adddata, "Observed values", NA),
                        ifelse(adddata, "Confidence interval", NA),
                        "Credible limits", "Mean curve"),
             bty = "n")
    }
  }, x = dtheoQ, y = dobs, z = dataCIm)
}

survFitPlotTKTDGG <- function(data, xlab, ylab, main, one.plot, spaghetti,
                              dataCIm, adddata, addlegend) {
  
  if (one.plot) {
    survFitPlotTKTDGGOnePlot(data, xlab, ylab, main, adddata, addlegend)
  } else {
    survFitPlotTKTDGGNoOnePlot(data, xlab, ylab, main, spaghetti,
                               dataCIm, adddata, concentration, addlegend)
  }
}

survFitPlotTKTDGGOnePlot <- function(data, xlab, ylab, main, adddata, addlegend) {
  gf <- ggplot(data[["dobs"]]) +
    geom_line(aes(x = time, y = q50, colour = factor(conc)),
              data = data[["dtheoQ"]]) +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
  if (adddata) {
    gf <- gf + geom_point(aes(x = time, y = psurv, colour = factor(conc)),
                          data = data[["dobs"]])
  }
  if (addlegend) {
    gf + scale_color_discrete("Concentrations")
  } else {
    gf + scale_color_discrete(guide = "none")
  }
}

survFitPlotTKTDGGNoOnePlot <- function(data, xlab, ylab, main, spaghetti,
                                       dataCIm, adddata, concentration,
                                       addlegend = FALSE) {
  if (!is.null(concentration)) {
    data[["dobs"]] <- filter(data[["dobs"]], conc == concentration)
    data[["dtheoQ"]] <- filter(data[["dtheoQ"]], conc == concentration)
    dataCIm <- filter(dataCIm, conc == concentration)
    curv_resp <- data.frame(time = data[["dtheoQ"]]$time,
                            resp = data[["dtheoQ"]]$q50,
                            Line = "Mean curve")
  }
  
  if (spaghetti) {
    gf <- ggplot(data[["dobs"]]) +
      geom_line(data = dataCIm, aes(x = time, y = value, group = variable),
                alpha = 0.05)
  } else {
    gf <- ggplot(data[["dobs"]]) +
      geom_ribbon(data = data[["dtheoQ"]], aes(x = time, ymin = qinf95,
                                               ymax = qsup95),
                  fill = "pink", col = "pink", alpha = 0.4)
  }
  
  if (!is.null(concentration)) {
    gf <- gf + geom_line(data = curv_resp, aes(x = time, y = resp,
                                               color = Line),
                         linetype = 'dashed')  +
      scale_color_hue("")
  } else {
    gf <- gf + geom_line(data = data[["dtheoQ"]], aes(x = time, y = q50),
                                                      color = data[["dtheoQ"]]$color,
                         linetype = 'dashed')
  }
  gf <- gf + geom_line(data = data[["dtheoQ"]], aes(x = time, y = qinf95,
                                                    color = Cred.Lim),
                       linetype = 'dashed') +
    geom_line(data = data[["dtheoQ"]], aes(x = time, y = qsup95,
                                           color = Cred.Lim),
              linetype = 'dashed') +
    facet_wrap(~conc) +
    scale_linetype(name = "") +
    labs(x = xlab, y = ylab) + ggtitle(main) +
    ylim(c(0, 1)) +
    theme_minimal()
  
  if (adddata) {
    # dataframes points (data) and curve (curv)
    gf <- gf +
      geom_point(aes(x = time, y = psurv, fill = Points),
                 data = data[["dobs"]], col = "black") +
      geom_segment(aes(x = time, xend = time, y = qinf95, yend = qsup95,
                       linetype = "Conf.Int"),
                   data[["dobs"]], col = "gray",
                   size = 0.5) +
      scale_fill_hue("")
  } else {
     gf <- gf
  }
  
  if (addlegend) {
    gf
  } else {
    gf + theme(legend.position = "none")
  }
  
}
