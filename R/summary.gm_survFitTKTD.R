#' Summary for \code{survFitTKTD} objects
#'
#' This is the generic \code{summary} S3 method for the \code{survFitTKTD} class.
#' It shows the quantiles of priors and posteriors on parameters.
#'
#' @param object an object of class \code{survFitTKTD}
#' @param quiet when \code{FALSE}, prints summary on standard output
#' @param \dots Further arguments to be passed to generic methods.
#'
#' @return The function returns a list with the following fields:
#' \item{Qpriors}{quantiles for the model's prior}
#' \item{Qposteriors}{quantiles for the model's posteriors}
#'
#' @examples
#' # (1) Load the data
#' data(propiconazole)
#'
#' # (2) Create a survData object
#' dat <- survData(propiconazole)
#'
#' \dontrun{
#' # (3) Run the survFitTKTD function
#' out <- survFitTKTD(dat)
#'
#' # (4) summarize the survFitTKTD object
#' summary(out)
#' }
#'
#' @keywords summary
#'
#' @importFrom stats qnorm qunif
#'
#' @export
summary.gm_survFitTKTD <- function(object, quiet = FALSE, ...) {

  # quantiles of priors parameters
  n.iter <- object$n.iter$end - object$n.iter$start

  modelData = object$modelData

  # kd
  kd_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
              mean = modelData$kd_meanlog10,
              sd = modelData$kd_sdlog10)

  kd <- 10^kd_log10

  # hb
  hb_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                    mean = modelData$hb_meanlog10,
                    sd = modelData$hb_sdlog10)

  hb <- 10^hb_log10


  if(object$model_type == "SD" ){
    # kk
    kk_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                      mean = modelData$kk_meanlog10,
                      sd = modelData$kk_sdlog10)

    kk <- 10^kk_log10

    # z
    z_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                     mean = modelData$z_meanlog10,
                     sd = modelData$z_sdlog10)

    z <- 10^z_log10

    res <- rbind(kd, hb, kk, z)

  }
  if(object$model_type == "IT" ){
    # alpha
    alpha_log10 <- qnorm(p = c(0.5, 0.025, 0.975),
                      mean = modelData$alpha_meanlog10,
                      sd = modelData$alpha_sdlog10)

    alpha <- 10^alpha_log10

    # z
    beta_log10 <- qunif(p = c(0.5, 0.025, 0.975),
                        min = modelData$beta_minlog10,
                        max = modelData$beta_maxlog10)

    beta <- 10^beta_log10

    res <- rbind(kd, hb, alpha, beta)
  }

  ans1 <- format(data.frame(res), scientific = TRUE, digits = 4)
  colnames(ans1) <- c("50%", "2.5%", "97.5%")

  # quantiles of estimated model parameters
  ans2 <- format(object$estim.par, scientific = TRUE, digits = 4)
  colnames(ans2) <- c("50%", "2.5%", "97.5%")

  # print
  if (! quiet) {
    cat("Summary: \n\n")
    cat("Model type:", object$model_type, "\n\n")
    cat("Priors on parameters (quantiles):\n\n")
    print(ans1)
    cat("\nPosterior of the parameters (quantiles):\n\n")
    print(ans2)
  }

  invisible(list(Qpriors = ans1,
                 Qposts = ans2))
}
