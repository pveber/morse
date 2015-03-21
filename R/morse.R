#' MOdelling tools for Reproduction and Survival data in Ecotoxicology
#'
#' The package offers tools for ecotoxicologists and regulators based on
#' advanced and innovative methods for a valuable quantitative environmental
#' risk assessment. The package allows the analysis of bioassay reproduction
#' data accounting for mortality all along the bioassay.  Such data are
#' commonly used to estimate Effective Concentration (\eqn{EC_{x}}) values from
#' chronic toxicity tests. The aim is to fit an exposure-response curve to
#' reproduction data by Bayesian inference while taking into account mortality
#' among parents without loosing valuable data (Delignette-Muller et al.,
#' 2014). Models are characterized by a deterministic log-logistic part
#' associated with a stochastic part. Two different stochastic parts can be
#' chosen: Poisson or Gamma-Poisson. The package also allows the analysis of
#' bioassay survival data at the final time. Such data are commonly used to
#' estimate Lethal Concentration (\eqn{LC_{x}}) values from chronic toxicity
#' tests. The aim is to fit an exposure-response curve to survival data by
#' Bayesian inference. Models are characterized by a deterministic log-logistic
#' part with 2 or 3 parameters associated with a binomial stochastic part. The
#' number of parameters depending on the mortality in the control at the final
#' time. The package uses the \code{rjags} package (Plummer, 2013), an
#' interface from R to the JAGS library for Bayesian data analysis. Note that
#' the \code{rjags} package does not include a copy of the JAGS library: you
#' must install it separately. For instructions on downloading JAGS, see the
#' home page at \url{http://mcmc-jags.sourceforge.net}.
#'
#' \tabular{ll}{ Package: \tab morse\cr Type: \tab Package\cr Version: \tab
#' 2.0.0\cr Date: \tab 2015-01-14\cr License: \tab GPL (>=2)\cr }
#'
#' @name morse-package
#' @aliases morse-package morse
#' @docType package
#' @author Marie Laure Delignette-Muller
#' <marielaure.delignettemuller@@vetagro-sup.fr>, Philippe Ruiz
#' <philippe.ruiz@@univ-lyon1.fr>, Sandrine Charles
#' <sandrine.charles@@univ-lyon1.fr>, Wandrille Duchemin
#' <wandrille.duchemin@@insa-lyon.fr>, Christelle Lopes
#' <christelle.lopes@@univ-lyon1.fr>, Guillaume Kon-Kam-king
#' <guillaume.kon-kam-king@@univ-lyon1.fr>, Philippe Veber
#' <philippe.veber@@univ-lyon1.fr>
#'
#' Maintainer: Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>
#' @seealso \code{\link[rjags]{rjags}}, \code{\link[lattice]{lattice}},
#' \code{\link[ggplot2]{ggplot}}
#' @references Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S.
#' (2014) \emph{Statistical handling of reproduction data for exposure-response
#' modelling}.
#' \url{http://pubs.acs.org/doi/abs/10.1021/es502009r?journalCode=esthag}.
#'
#' Plummer, M. (2013) \emph{JAGS Version 3.4.0 user manual}.
#' \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/3.x/jags_user_manual.pdf/download}.
#' @keywords package
#FIXME @examples
#
# # (1) Load the data
# data(cadmium1)
#
# # Reproduction analysis at final time with log-logistic model
# # on reproduction data
#
# # (2) Check data
# reproDataCheck(cadmium1)
#
# # (3) Create an object of class "reproData"
# dat1 <- reproData(cadmium1)
#
# # (4) Plot raw data
# survPlotTt(dat1, log.scale = TRUE, type = "generic")
# survFullPlot(cadmium1, type = "lattice")
# reproCumulPlotTt(dat1, type = "ggplot")
#
# \dontrun{
# # (5) Fit the exposure-response log-logistic model on reproduction data at
# # final time
# out1 <- reproFitTt(dat1, n.chains = 3)
#
# # (6) Check the mcmc convergence
# convergence(out1, type = "generic")
#
# # (7) Summarize the results
# plot(out1, ci = TRUE, log.scale = TRUE)
# summary(out1)
# print(out1)
# }
#
# # Survival analysis at final time with log-logistic binomial
# # model on survival data
#
# # (8) Create an object of class "survData"
# dat2 <- survData(cadmium1)
#
# \dontrun{
# # (9) Fit the exposure-response log-logistic model on survival data
# out2 <- survFitTt(dat2, det.part = "loglogisticbinom_3")
#
# # (12) Check the mcmc convergence
# convergence(out2, type = "generic")
#
# # (13) Summarize the results
# plot(out2, ci = TRUE, log.scale = TRUE,
# pool.replicate = FALSE)
# summary(out2)
# print(out2)
# }
#
NULL


#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' cadmium during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of one metal contaminant (cadmium) during 21 days. Five concentrations were
#' tested, with four replicates per concentration. Each replicate contained 10
#' organisms. Reproduction and survival were monitored at 10 time points.
#'
#'
#' @name cadmium1
#' @docType data
#' @usage data(cadmium1)
#' @format A data frame with 200 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{integer} with
#' the replicate code (\code{1} to \code{4}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E., Delhaye, H., Forfait, C., Clement, B.,
#' Triffault-Bouchet, G., Charles, S. and Delignette-Muller, M.L. (2012)
#' Comparison of bioassays with different exposure time patterns: The added
#' value of dynamic modelling in predictive ecotoxicology, \emph{Ecotoxicology
#' and Environmental Safety}, 75, 80-86.
#' @keywords datasets
#FIXME @examples
#
# # (1) Load the data
# data(cadmium1)
#
# # (2) Create an object of class "reproData"
# dat <- reproData(cadmium1)
#
# # (3) Plot the number of survivors during time for each concentration
# survFullPlot(cadmium1)
#
# # (4) Plot the survival data depending on the concentration
# survPlotTt(dat, log.scale = TRUE)
#
# # (5) Plot the cumulated number of offspring depending on the concentration
# reproCumulPlotTt(dat, log.scale = TRUE)
#'
#'
NULL





#' Reproduction and survival datasets for snails exposed to cadmium during 56
#' days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' snails exposed to six concentrations of one metal contaminant (cadmium)
#' during 56 days. Six concentrations were tested, with six replicates per
#' concentration. Each replicate contained five organisms. Reproduction and
#' survival were monitored at 17 time points.
#'
#'
#' @name cadmium2
#' @docType data
#' @usage data(cadmium2)
#' @format A data frame with 612 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{F}).} \item{\code{conc}}{A vector of
#' class \code{integer} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @keywords datasets
# @examples
#
# # (1) Load the data
# data(cadmium2)
#
# # (2) Create an object of class "reproData"
# dat <- reproData(cadmium2)
#
# # (3) Plot the number of survivors during time for each concentration
# survFullPlot(cadmium2)
#
# # (4) Plot the survival data depending on the concentration
# survPlotTt(dat, log.scale = TRUE)
#
# # (5) Plot the cumulated number of offspring depending on the concentration
# reproCumulPlotTt(dat, log.scale = TRUE)
#
#
NULL





#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' chlordan during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to six concentrations
#' of one organochlorine insecticide during 21 days. Six concentrations were
#' tested, with 10 replicates per concentration. Each replicate contained one
#' organism. Reproduction and survival were monitored at 22 time points.
#'
#'
#' @name chlordan
#' @docType data
#' @usage data(chlordan)
#' @format A data frame with 1320 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{integer} with
#' the replicate code (\code{1} to \code{10}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the chlordan concentrations in \eqn{\mu
#' g.L^{-1}}.} \item{\code{time}}{A vector of class \code{integer} with the
#' time points (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Manar, R., Bessi, H. and Vasseur, P. (2009) Reproductive effects
#' and bioaccumulation of chlordan in Daphnia magna, \emph{Environmental
#' Toxicology and Chemistry}, 28, 2150-2159.
#' @keywords datasets
# @examples
#
# # (1) Load the data
# data(chlordan)
#
# # (2) Create an object of class "reproData"
# dat <- reproData(chlordan)
#
# # (3) Plot the number of survivors during time for each concentration
# survFullPlot(chlordan)
#
# # (4) Plot the survival data depending on the concentration
# survPlotTt(dat, log.scale = TRUE)
#
# # (5) Plot the cumulated number of offspring depending on the concentration
# reproCumulPlotTt(dat, log.scale = TRUE)
#
#
NULL





#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to
#' copper during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of one metal contaminant (copper) during 21 days. Five concentrations were
#' tested, with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 16 time points.
#'
#'
#' @name copper
#' @docType data
#' @usage data(copper)
#' @format A data frame with 240 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{C}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with the copper concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E., Delignette-Muller, M.L., Pery, A.R.R. and
#' Charles, S. (2008) A Bayesian Approach to Analyzing Ecotoxicological Data,
#' \emph{Environmental Science & Technology}, 42 (23), 8978-8984.
#' @keywords datasets
# @examples
#
# # (1) Load the data
# data(copper)
#
# # (2) Create an object of class "reproData"
# dat <- reproData(copper)
#
# # (3) Plot the number of survivors during time for each concentration
# survFullPlot(copper)
#
# # (4) Plot the survival data depending on the concentration
# survPlotTt(dat, log.scale = TRUE)
#
# # (5) Plot the cumulated number of offspring depending on the concentration
# reproCumulPlotTt(dat, log.scale = TRUE)
#
#
NULL









#' Reproduction and survival datasets for \emph{Daphnia magna} exposed to zinc
#' during 21 days
#'
#' Reproduction and survival datasets of chronic laboratory bioassays with
#' \emph{Daphnia magna} freshwater invertebrate exposed to four concentrations
#' of one metal contaminant (zinc) during 21 days. Four concentrations were
#' tested with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 15 time points.
#'
#'
#' @name zinc
#' @docType data
#' @usage data(zinc)
#' @format A data frame with 180 observations on the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{A} to \code{C}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with zinc concentrations in \eqn{mg.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of offspring at each time point for each concentration and each
#' replicate.} }
#' @references Billoir, E.,Delignette-Muller, M.L., Pery, A.R.R. and
#' Charles S. (2008) A Bayesian Approach to Analyzing Ecotoxicological Data,
#' \emph{Environmental Science & Technology}, 42 (23), 8978-8984.
#' @keywords datasets
# @examples
#
# # (1) Load the data
# data(zinc)
#
# # (2) Create an object of class "reproData"
# dat <- reproData(zinc)
#
# # (3) Plot the number of survivors during time for each concentration
# survFullPlot(zinc)
#
# # (4) Plot the survival data depending on the concentration
# survPlotTt(dat, log.scale = TRUE)
#
# # (5) Plot the cumulated number of offspring depending on the concentration
# reproCumulPlotTt(dat, log.scale = TRUE)
#
#
#
NULL



