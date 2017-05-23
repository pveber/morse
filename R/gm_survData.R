#' Creates a dataset for survival analysis
#'
#' This function creates a \code{gm_survData} object from experimental data
#' provided as two \code{data.frame}s. The resulting object can then be used for
#' plotting and model fitting. It can also be used to generate
#' \emph{individual-time} estimates.
#'
#' The \code{data} argument describes experimental results from a survival
#' assay. Each line of the \code{data.frame} corresponds to one experimental
#' measurement, that is a number of alive individuals for a given concentration
#' of pollutant at a certain time during the assay in a certain replicate. The
#' function fails if \code{data} does not meet the expected requirements. Please
#' run \code{\link{gm_survDataCheck}} to ensure \code{data} is well-formed.
#'
#' @param data_surv a \code{data.frame} containing the following four columns:
#'   \itemize{ #' \item \code{profile}: a vector of class \code{integer},
#'   \code{character} or \code{factor} for profile identification \item
#'   \code{replicate}: a vector of class \code{integer} or \code{factor} for
#'   replicate identification \item \code{time}: a vector of class
#'   \code{integer} with time points, min value must be 0 \item \code{Nsurv}: a
#'   vector of class \code{integer} providing the number of alive individuals at
#'   each time point for each profile and each replicate }
#' @param data_conc a \code{data.frame} containing the following four columns:
#'   \itemize{ #' \item \code{profile}: a vector of class \code{integer},
#'   \code{character} or \code{factor} for profile identification \item
#'   \code{replicate}: a vector of class \code{integer} or \code{factor} for
#'   replicate identification \item \code{conc}: a vector of class
#'   \code{numeric} with tested concentrations (positive values) at each time
#'   point for each profile and each replicate \item \code{time}: a vector of
#'   class \code{integer} with time points, min value must be 0 }
#'
#'
#' @return A dataframe of class \code{gm_survData}.
#'
#' @seealso \code{\link{gm_survDataCheck}}
#'
#' @keywords transformation
#'
#' @examples
#'
#' # (1) Load the survival dataset
#' data(zinc)
#'
#' # (2) Create an objet of class 'gm_survData'
#' dat <- gm_survData(zinc)
#' class(dat)
#'
#'
#' @export

gm_survData = function(data_surv, data_conc){

  ##
  ## Test the integrity of the data with gm_survDataCheck
  ##
  gm_survDataCheck(data_surv, data_conc)

  ##
  ## Sum replicate in 'data_surv'
  ## + arrange in 'profile' and 'time'
  ##

  data_surv_sumReplicate = data_surv %>%
    dplyr::group_by(profile, time) %>%
    dplyr::summarise(Nsurv = sum(Nsurv)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(profile, time)%>%
    # 'lag' function copy values lagged by 1 (see 'dplyr' package)
    dplyr::mutate( tprec = ifelse( time == 0, time, dplyr::lag(time) ) ) %>%
    dplyr::mutate( Nprec = ifelse( time == 0, Nsurv, dplyr::lag(Nsurv) ) ) %>%
    ungroup()

  ##
  ## Join data_surv and data_conc
  ##

  df.survData = dplyr::full_join(data_surv_sumReplicate,
                          data_conc,
                          by=c("profile","time")) %>%
    dplyr::arrange(profile, time) %>%
    dplyr::group_by(profile) %>%
    # Extrapolate last conc until end time in data_surv if necessary.
    # But we do not interpolate otherwise.
    # From package zoo : 'na.locf()' carry the last observation forward to replace your NA values.
    # We create a variable 'conc.interpol' to test if variable is interpolated or extrapolated
    #mutate(conc.origin = conc)%>%
    mutate(conc.interpol = zoo::na.approx(conc,time, na.rm = FALSE))%>%
    dplyr::mutate(conc = ifelse(is.na(conc.interpol), zoo::na.locf(conc, na.rm = FALSE), conc)) %>%
    select(-conc.interpol) %>%
    ungroup(profile, time)

  class(df.survData) <- c("gm_survData", "data.frame")

  return(df.survData)
  # return(list(df.survData = df.survData,
  #             df.survOnly = data_surv_sumReplicate,
  #             df.concOnly = data_conc_arrange))
}
