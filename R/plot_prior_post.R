#' Generic method to plot priors and posteriors.
#' 
#' Plot priors and posteriors of a \code{survFit} object
#' 
#' 
#' @param x an object used to select a method \code{plot_prior_post}
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
plot_prior_post <- function(x, ...){
  UseMethod("plot_prior_post")
}
#' Plot posteriors vs priors
#' 
#' Plot posteriors vs priors of a \code{survFit} object
#' 
#' 
#' @param x an object of class \code{survFit} used to select a method \code{plot_prior_post}
#' @param size_sample Size of the random generation of the distribution.
#' Default is \code{1e3}.
#' @param EFSA_name If \code{TRUE}, replace the current terminology by
#'  the one used in the recent EFSA PPR Scientific Opinion (2018).
#' 
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @references
#' EFSA PPR Scientific Opinion (2018)
#' \emph{Scientific Opinion on the state of the art of Toxicokinetic/Toxicodynamic (TKTD) effect models for regulatory risk assessment of pesticides for aquatic organisms}
#' \url{https://www.efsa.europa.eu/en/efsajournal/pub/5377}
#' 
#' @export
plot_prior_post.survFit <- function(x, size_sample = 1e3, EFSA_name = FALSE, ...){

  priors_distr <- priors_distribution(object = x, size_sample = size_sample) %>%
    select(contains("log10")) 

  mctot <- do.call("rbind", x$mcmc)
  
  if(!("hb_log10" %in% colnames(mctot))){
    posteriors_distr <- data.frame(
      kd_log10 = mctot[, "kd_log10"]
    )
    priors_distr <- dplyr::select(priors_distr, - "hb_log10")
  } else{
    posteriors_distr <- data.frame(
      kd_log10 = mctot[, "kd_log10"],
      hb_log10 = mctot[, "hb_log10"]
    )
  }
  
  if(x$model_type == "SD"){
    posteriors_distr$z_log10 = mctot[, "z_log10"]
    posteriors_distr$kk_log10 = mctot[, "kk_log10"]
  }
  if(x$model_type == "IT"){
    posteriors_distr$alpha_log10 = mctot[, "alpha_log10"]
    posteriors_distr$beta_log10 = mctot[, "beta_log10"]
  }

  if(EFSA_name == TRUE){
    if(x$model_type == "IT"){
      posteriors_distr <- dplyr::rename(posteriors_distr,
                    kD_log10 = kd_log10,
                    mw_log10 = alpha_log10)
      priors_distr <- dplyr::rename(priors_distr,
                    kD_log10 = kd_log10,
                    mw_log10 = alpha_log10)
    }
    if(x$model_type == "SD"){
      posteriors_distr <- dplyr::rename(posteriors_distr,
                    kD_log10 = kd_log10,
                    bw_log10 = kk_log10,
                    zw_log10 = z_log10)
      priors_distr <- dplyr::rename(priors_distr,
                    kD_log10 = kd_log10,
                    bw_log10 = kk_log10,
                    zw_log10 = z_log10)
    }
  }
  
  priors_distr <- priors_distr %>%
    gather(parameters, density)
 
  
  posteriors_distr <- posteriors_distr %>%
    gather(parameters, density)
  
  plt <- ggplot() + theme_minimal() + 
    geom_density(data = priors_distr,
                 aes(density, group = parameters), fill = "grey", color = NA) +
    geom_density(data = posteriors_distr,
                 aes(density, group = parameters), fill = "orange", color = NA) +
    facet_wrap(~ parameters, nrow = 1, scales = "free")
  
  return(plt)
}