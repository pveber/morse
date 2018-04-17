priors_distribution <- function(object, size_sample = 1e3){
  
  param <- object$jags.data
  
  # kd
  kd_log10 <- rnorm(size_sample,
                    mean = param$kd_meanlog10,
                    sd = param$kd_sdlog10)
  
  kd <- 10^kd_log10
  
  
  # hb
  hb_log10 <- rnorm(size_sample,
                    mean = param$hb_meanlog10,
                    sd = param$hb_sdlog10)
  
  hb <- 10^hb_log10
  
  df = data.frame(hb = hb,
                  hb_log10 = hb_log10,
                  kd = kd,
                  kd_log10 = kd_log10)
  
  
  if(object$model_type == "SD"){
    
    # kk
    kk_log10 <- rnorm(size_sample,
                      mean = param$kk_meanlog10,
                      sd = param$kk_sdlog10)
    
    kk <- 10^kk_log10
    
    ## z
    z_log10 <- rnorm(size_sample,
                     mean = param$z_meanlog10,
                     sd = param$z_sdlog10)
    
    z <- 10^z_log10
    
    df$kk = kk
    df$kk_log10 = kk_log10
    df$z = z
    df$z_log10 = z_log10
    
  }
  if(object$model_type == "IT"){
    
    # alpha
    alpha_log10 <- rnorm(size_sample,
                         mean = param$alpha_meanlog10,
                         sd = param$alpha_sdlog10)
    
    alpha <- 10^alpha_log10
    
    # beta
    beta_log10 <- runif(size_sample,
                        min = param$beta_minlog10,
                        max = param$beta_maxlog10)
    
    beta <- 10^beta_log10

    df$alpha = alpha
    df$alpha_log10 = alpha_log10
    df$beta = beta
    df$beta_log10 = beta_log10
  }
  if(object$model_type == "PROPER"){
    
    # kk
    kk_log10 <- rnorm(size_sample,
                      mean = param$kk_meanlog10,
                      sd = param$kk_sdlog10)
    
    kk <- 10^kk_log10
    
    # alpha
    alpha_log10 <- rnorm(size_sample,
                         mean = param$alpha_meanlog10,
                         sd = param$alpha_sdlog10)
    
    alpha <- 10^alpha_log10
    
    # beta
    beta_log10 <- runif(size_sample,
                        min = param$beta_minlog10,
                        max = param$beta_maxlog10)
    
    beta <- 10^beta_log10
    
    df$kk = kk
    df$kk_log10 = kk_log10
    
    df$alpha = alpha
    df$alpha_log10 = alpha_log10
    df$beta = beta
    df$beta_log10 = beta_log10
  }

  return(df)
}
