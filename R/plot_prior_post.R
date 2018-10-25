#' Plot posteriors vs priors
#' 
#' Plot posteriors vs priors of a \code{survFit} object
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
#' @param x an object used to select a method \code{plot_prior_post}
#' @param size_sample Size of the random generation of the distribution.
#' Default is \code{1e3}.
#' @param \dots Further arguments to be passed to generic methods
#' 
#' @export
plot_prior_post.survFit <- function(x, size_sample = 1e3, ...){

  priors_distr <- priors_distribution(object = x, size_sample = size_sample) %>%
    select(contains("log10")) %>%
    gather(parameters, density)
  
  mctot <- do.call("rbind", x$mcmc)
  
  posteriors_distr <- data.frame(
    kd_log10 = mctot[, "kd_log10"],
    hb_log10 = mctot[, "hb_log10"]
  )
  if(x$model_type == "SD"){
    posteriors_distr$z_log10 = mctot[, "z_log10"]
    posteriors_distr$kk_log10 = mctot[, "kk_log10"]
  }
  if(x$model_type == "IT"){
    posteriors_distr$alpha_log10 = mctot[, "alpha_log10"]
    posteriors_distr$beta_log10 = mctot[, "beta_log10"]
  }
  
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