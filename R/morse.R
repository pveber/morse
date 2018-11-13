#' MOdelling tools for Reproduction and Survival data in Ecotoxicology
#'
#' Provides tools for the analysis of survival/reproduction
#' toxicity test data in quantitative environmental risk assessment. It can be
#' used to explore/visualize experimental data, and to get estimates
#' of \eqn{LC_{x}} (\eqn{X}\% Lethal Concentration) or,
#' \eqn{EC_{x}} (\eqn{X}\% Effective Concentration) by fitting exposure-response
#' curves. The \eqn{LC_{x}}, \eqn{EC_{x}} and parameters of the curve are
#' provided along with an indication of the uncertainty of the estimation.
#' \code{morse} can also be used to get an estimation of the \eqn{NEC} (No Effect Concentration)
#' by fitting a Toxico-Kinetic Toxico-Dynamic (TKTD) model (GUTS: General Unified Threshold
#' model of Survival). Within the TKTD-GUTS approach, \eqn{LC(x,t)}, \eqn{EC(x,t)} and \eqn{MF(x,t)}
#' (\eqn{x}\% Multiplication Factors aka Lethal Profiles) can be explored in proportion \eqn{x} and 
#' time \eqn{t}.
#'
#' Estimation procedures in \code{morse} can be used without a deep knowledge of
#' their underlying probabilistic model or inference methods. Rather, they
#' were designed to behave as well as possible without requiring a user to
#' provide values for some obscure parameters. That said, \code{morse} models can also
#' be used as a first step to tailor new models for more specific situations.
#'
#' The package currently handles survival and reproduction data. Functions
#' dedicated to survival (resp. reproduction) analysis start with a
#' \code{surv} (resp. \code{repro}) prefix. \code{morse} provides a similar
#' workflow in both cases:
#' \enumerate{
#' \item create and validate a data set
#' \item explore a data set
#' \item plot a data set
#' \item fit a model on a data set and output the expected estimates
#' \item check goodness of fit with posterior preditive check plot (ppc)
#' }
#' 
#' More specifically, for survival data handles with TKTD `GUTS` model, \code{morse}
#' provides:
#' \enumerate{
#' \item plot \eqn{LC(x,t)} and \eqn{MF(x,t)}.
#' \item compute goodness-of-fit measures (PPC percent, NRMSE and SPPE)
#' }
#' 
#' Those steps are presented in more details in the "Tutorial" vignette, while
#' a more formal description of the estimation procedures are provided in the
#' vignette called "Models in \code{morse} package". Please refer to these documents
#' for further introduction to the use of \code{morse}.
#'
#' This reference manual is a detailed description of the functions exposed in
#' the package.
#'
#' \strong{Getting started} The package uses the \code{rjags} package
#' (Plummer, 2013), an R interface to the JAGS library for Bayesian model
#' estimation. Note that the \code{rjags} package does not include a copy
#' of the JAGS library: you need to install it separately. For instructions
#' on downloading JAGS, see the home page at
#' \url{http://mcmc-jags.sourceforge.net}. Once done, simply follow the steps
#' described in the tutorial vignette.
#'
#' \tabular{ll}{ Package: \tab morse\cr Type: \tab Package\cr Version: \tab
#' 3.2.0\cr Date: \tab 2018-11-15\cr License: \tab GPL (>=2)\cr }
#'
#' @name morse-package
#' @aliases morse-package morse
#' @docType package
#' @author
#' Virgile Baudrot  <virgile.baudrot@@posteo.net>,
#' Sandrine Charles <sandrine.charles@@univ-lyon1.fr>,
#' Marie Laure Delignette-Muller <marielaure.delignettemuller@@vetagro-sup.fr>,
#' Wandrille Duchemin <wandrille.duchemin@@insa-lyon.fr>,
#' Benoit Goussen <Benoit.Goussen@@ibacon.com>,
#' Guillaume Kon-Kam-king <guillaume.kon-kam-king@@univ-lyon1.fr>,
#' Christelle Lopes <christelle.lopes@@univ-lyon1.fr>,
#' Philippe Ruiz <philippe.ruiz@@univ-lyon1.fr>,
#' Philippe Veber <philippe.veber@@univ-lyon1.fr>
#' 
#' Maintainer: Philippe Veber <philippe.veber@@univ-lyon1.fr>
#' @seealso \code{\link[rjags]{rjags}},
#' \code{\link[ggplot2]{ggplot2}}
#' 
#' @references 
#' Delignette-Muller, M.L., Lopes, C., Veber, P. and Charles, S. (2014)
#' \emph{Statistical handling of reproduction data for exposure-response modelling}.
#' \url{http://pubs.acs.org/doi/abs/10.1021/es502009r?journalCode=esthag}.
#' 
#' Forfait-Dubuc, C., Charles, S., Billoir, E. and Delignette-Muller, M.L. (2012)
#' \emph{Survival data analyses in ecotoxicology: critical effect concentrations, methods and models. What should we use?}
#' \url{https://doi.org/10.1007/s10646-012-0860-0}.
#'
#' Plummer, M. (2013) \emph{JAGS Version 4.0.0 user manual}.
#' \url{http://sourceforge.net/projects/mcmc-jags/files/Manuals/4.x/jags_user_manual.pdf/download}
#' 
#' Delignette-Muller, M. L., Ruiz, P. and Veber, P. (2017)
#' \emph{Robust Fit of Toxicokinetic--Toxicodynamic Models Using Prior Knowledge Contained in the Design of Survival Toxicity Tests}
#' \url{https://pubs.acs.org/doi/abs/10.1021/acs.est.6b05326}
#'
#' Baudrot, V., Preux, S., Ducrot, V., Pavé, A. and Charles, S. (2018)
#' \emph{New insights to compare and choose TKTD models for survival based on an inter-laboratory study for \emph{Lymnaea stagnalis} exposed to Cd}.
#' \url{https://pubs.acs.org/doi/abs/10.1021/acs.est.7b05464}.
#' 
#' EFSA PPR Scientific Opinion (2018)
#' \emph{Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms}
#' \url{https://www.efsa.europa.eu/en/efsajournal/pub/5377}.
#' 
NULL


#' Reproduction and survival data sets for \emph{Daphnia magna} exposed to
#' cadmium during 21 days
#'
#' Reproduction and survival data sets of chronic laboratory toxicity tests with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of cadmium during 21 days. Five concentrations were
#' tested, with four replicates per concentration. Each replicate contained 10
#' organisms. Reproduction and survival were monitored at 10 time points.
#'
#'
#' @name cadmium1
#' @docType data
#' @usage data(cadmium1)
#' @format A data frame with 200 observations of the following five variables:
#' \describe{
#' \item{\code{replicate}}{A vector of class \code{numeric} with
#' the replicate code (\code{1} to \code{20}).}
#' \item{\code{conc}}{A vector of
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
#' Comparison of toxicity tests with different exposure time patterns: The added
#' value of dynamic modelling in predictive ecotoxicology, \emph{Ecotoxicology
#' and Environmental Safety}, 75, 80-86.
#' @keywords data set
NULL




#' Reproduction and survival data sets for \emph{Lymnaea stagnalis} exposed to cadmium during 28
#' days
#'
#' Reproduction and survival data sets of chronic laboratory toxicity tests with
#' snails (\emph{Lymnaea stagnalis}) exposed to six concentrations of cadmium
#' during 28 days. Six concentrations were tested, with six replicates per
#' concentration. Each replicate contained five organisms. Reproduction and
#' survival were monitored at 17 time points.
#'
#'
#' @name cadmium2
#' @docType data
#' @usage data(cadmium2)
#' @format A data frame with 612 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{numeric} with the
#' replicate code (\code{1} to \code{36}).} \item{\code{conc}}{A vector of
#' class \code{integer} with the cadmium concentrations in \eqn{\mu g.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.} \item{\code{Nrepro}}{A vector of class \code{integer} with the
#' number of clutches at each time point for each concentration and each
#' replicate.} }
#' @references Ducrot, V., Askem, C., Azam, D., Brettschneider, D., Brown,
#' R., Charles, S., Coke, M., Collinet, M., Delignette-Muller, M.L.,
#' Forfait-Dubuc, C., Holbech, H., Hutchinson, T., Jach, A., Kinnberg, K.L.,
#' Lacoste, C., Le Page, G., Matthiessen, P., Oehlmann, J., Rice, L.,
#' Roberts, E., Ruppert, K., Davis, J.E., Veauvy, C., Weltje, L., Wortham, R.
#' and Lagadic, L. (2014)
#' Development and validation of an OECD reproductive toxicity test guideline with
#' the pond snail Lymnaea stagnalis (Mollusca, Gastropoda),
#' \emph{Regulatory Toxicology and Pharmacology}, 70(3), 605-14.
#' 
#' Charles, S., Ducrot, V., Azam, D., Benstead, R., Brettschneider, D., De Schamphelaere, K.,
#' Filipe Goncalves, S., Green, J.W., Holbech, H., Hutchinson, T.H., Faber, D., Laranjeiro, F.,
#' Matthiessen, P., Norrgren, L., Oehlmann, J., Reategui-Zirena, E., Seeland-Fremer, A., Teigeler, M.,
#' Thome, J.P., Tobor Kaplon, M., Weltje, L., Lagadic, L. (2016) 
#' Optimizing the design of a reproduction toxicity test with the pond snail Lymnaea stagnalis,
#' \emph{Regulatory Toxicology and Pharmacology}, vol. 81 pp.47-56.
#' 
#' @keywords data set
NULL




#' Reproduction and survival data sets for \emph{Daphnia magna} exposed to
#' chlordan during 21 days
#'
#' Reproduction and survival data sets of chronic laboratory toxicity tests with
#' \emph{Daphnia magna} freshwater invertebrate exposed to six concentrations
#' of one organochlorine insecticide (chlordan) during 21 days. Six concentrations were
#' tested, with 10 replicates per concentration. Each replicate contained one
#' organism. Reproduction and survival were monitored at 22 time points.
#'
#'
#' @name chlordan
#' @docType data
#' @usage data(chlordan)
#' @format A data frame with 1320 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{numeric} with
#' the replicate code (\code{1} to \code{60}).} \item{\code{conc}}{A vector of
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
#' @keywords data set
NULL




#' Reproduction and survival data sets for \emph{Daphnia magna} exposed to
#' copper during 21 days
#'
#' Reproduction and survival data sets of chronic laboratory toxicity tests with
#' \emph{Daphnia magna} freshwater invertebrate exposed to five concentrations
#' of copper during 21 days. Five concentrations were
#' tested, with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 16 time points.
#'
#'
#' @name copper
#' @docType data
#' @usage data(copper)
#' @format A data frame with 240 observations of the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{numeric} with the
#' replicate code (\code{1} to \code{15}).} \item{\code{conc}}{A vector of
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
#' @keywords data set
NULL




#' Survival data set for \emph{Daphnia magna} exposed to dichromate
#' during 21 days
#'
#' Survival data set of chronic laboratory toxicity tests with
#' \emph{Daphnia magna} freshwater invertebrate exposed to six concentrations
#' of one oxidizing agent (potassium dichromate) during 21 days. Six
#' concentrations were tested with one replicate of 50 organisms per concentration.
#' Survival is monitored at 10 time points.
#'
#'
#' @name dichromate
#' @docType data
#' @usage data(dichromate)
#' @format A data frame with 60 observations on the following four variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{numeric} with the
#' replicate code (\code{1}).} \item{\code{conc}}{A vector of
#' class \code{numeric} with dichromate concentrations in \eqn{mg.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.}}
#' @references Bedaux, J., Kooijman, SALM (1994) Statistical analysis of toxicity tests,
#' based on hazard modeling, \emph{Environmental and Ecological Statistics}, 1,
#' 303-314.
#' @keywords data set
NULL



#' Survival data set for \emph{Gammarus pulex} exposed to propiconazole
#' during four days
#'
#' Survival data set of chronic laboratory toxicity tests with
#' \emph{Gammarus pulex} freshwater invertebrate exposed to eight concentrations
#' of one fungicide (propiconazole) during four days. Eight
#' concentrations were tested with two replicates of 10 organisms per concentration.
#' Survival is monitored at five time points.
#'
#'
#' @name propiconazole
#' @docType data
#' @usage data(propiconazole)
#' @format A dataframe with 75 observations on the following four variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{SC} for the control and \code{A1} to \code{G2} for other profiles).}
#' \item{\code{conc}}{A vector of class \code{numeric} with propiconazole
#' concentrations in \eqn{\mu mol.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.}}
#' @references Nyman, A.-M., Schirmer, K., Ashauer, R., (2012) Toxicokinetic-toxicodynamic
#' modelling of survival of \emph{Gammarus pulex} in multiple pulse exposures to
#' propiconazole: model assumptions, calibration data requirements and predictive
#' power, \emph{Ecotoxicology}, (21), 1828-1840.
#'
#' @keywords data set
NULL


#' Survival data set for \emph{Gammarus pulex} exposed to propiconazole
#' during 10 days with time-variable
#' exposure concentration (non-standard pulsed toxicity experiments)
#'
#' Survival data set of laboratory toxicity tests with \emph{Gammarus pulex}
#' freshwater invertebrates exposed to several profiles of concentrations
#' (time-variable concentration for each time series)
#' of one fungicide (propiconazole) during 10 days.
#'
#' @name propiconazole_pulse_exposure
#' @docType data
#' @usage data(propiconazole_pulse_exposure)
#' @format A data frame with 74 observations on the following four variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{factor} with the
#' replicate code (\code{varControl}, \code{varA}, \code{varB} and \code{varC}).}
#' \item{\code{conc}}{A vector of class \code{numeric} with propiconazole
#' concentrations in \eqn{\mu mol.L^{-1}}.}
#' \item{\code{time}}{A vector of class \code{integer} with the time points
#' (in days from the beginning of the experiment \eqn{t = 0}).}
#' \item{\code{Nsurv}}{A vector of class \code{integer} with the number of
#' alive individuals at each time point for each concentration and each
#' replicate.}}
#' @references Nyman, A.-M., Schirmer, K., Ashauer, R., (2012) Toxicokinetic-toxicodynamic
#' modelling of survival of \emph{Gammarus pulex} in multiple pulse exposures to
#' propiconazole: model assumptions, calibration data requirements and predictive
#' power, \emph{Ecotoxicology}, (21), 1828-1840.
#'
#' @keywords data set
NULL


#' Reproduction and survival data sets for \emph{Daphnia magna} exposed to zinc
#' during 21 days
#'
#' Reproduction and survival data sets of a chronic laboratory toxicity tests with
#' \emph{Daphnia magna} freshwater invertebrate exposed to four concentrations
#' of zinc during 21 days. Four concentrations were
#' tested with three replicates per concentration. Each replicate contained 20
#' organisms. Reproduction and survival were monitored at 15 time points.
#'
#'
#' @name zinc
#' @docType data
#' @usage data(zinc)
#' @format A data frame with 180 observations on the following five variables:
#' \describe{ \item{\code{replicate}}{A vector of class \code{numeric} with the
#' replicate code (\code{1} to \code{12}).} \item{\code{conc}}{A vector of
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
#' @keywords data set
NULL
