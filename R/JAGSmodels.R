#
# Here, it's exactly the JAGS model used in `MORSE` package in the function `survFitTKTD()`
#

jags_TKTD_cstSD <-
"model {
  ######### priors

  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  kk_taulog10 <- 1 / kk_sdlog10^2
  z_taulog10 <- 1 / z_sdlog10^2

  kk_log10 ~ dnorm(kk_meanlog10, kk_taulog10)
  z_log10  ~ dnorm(z_meanlog10 , z_taulog10)
  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  ##### parameter transformation
  kk <- 10**kk_log10
  z  <- 10**z_log10
  ke <- 10**kd_log10
  hb <- 10**hb_log10

  ########## Computation of the likelihood

  for (i in 1:n_data){

    tz[i] <- ifelse(x[i] > z, -1/kd * log( 1- R[i]), bigtime)
    R[i] <- ifelse(x[i] > z, z/xcor[i], 0.1)
    xcor[i] <- ifelse(x[i] > 0, x[i], 10)
    tref[i] <- max(tprec[i], tz[i])

    psurv[i] <- exp(-m0 * (t[i] - tprec[i]) + ifelse(t[i] > tz[i], -ks * ((x[i] - NEC) * (t[i] - tref[i]) + x[i]/kd * ( exp(-kd * t[i]) - exp(-kd * tref[i]))), 0))

    Nsurv[i] ~ dbin(psurv[i] , Nprec[i])

    ## ---------------------- generated data

    Nsurv_ppc[i] ~ dbin(psurv[i] , Nprec[i])
  }

  ### initialization is requires to use 'Nsurv_sim[i-1]' require in JAGS language (avoid auto-loop issue).
  Nsurv_sim[1] ~ dbin(psurv[1], Nprec[1])
  for( i in 2:n_data){
    Nsurv_sim[i] ~ dbin(psurv[i], ifelse( i == i_prec[i], Nprec[i], Nsurv_sim[i-1]))
  }
}"


#
# IT model with log-normal function.
#

jags_TKTD_cstIT <-"
model {

  #-------------------------------------------------------- priors

  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  alpha_taulog10 <- 1 / alpha_sdlog10^2

  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  alpha_log10 ~ dnorm(alpha_meanlog10, alpha_taulog10)
  beta_log10 ~ dunif(beta_minlog10, beta_maxlog10)

  #------------------------------------------ parameter transformation

  kd <- 10**kd_log10
  hb <- 10**hb_log10

  alpha <- 10**alpha_log10
  beta <- 10**beta_log10

  ##------------------------------------ model

  for( i in 1:n_data){

    D[profile_ID[i], time_ID[i]] <- concentration[i] * ( 1 - exp( - kd * time[i] ))

    D_max[profile_ID[i], time_ID[i]] <- max(D[profile_ID[i],1:time_ID[i]])

    F[i]  <- D_max[profile_ID[i], time_ID[i]]**beta / ( D_max[profile_ID[i], time_ID[i]]**beta + alpha**beta )

    psurv[i] <-  exp(-hb * time[i]) * (1- F[i])

    Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

    ## ---------------------- generated data

    Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

  }

  ### initialization is requires to use 'Nsurv_sim[i-1]' require in JAGS language (avoid auto-loop issue).
  Nsurv_sim[1] ~ dbin(psurv[1]/psurv[1], Nprec[1])
  for( i in 2:n_data){
    Nsurv_sim[i] ~ dbin(psurv[i]/psurv[i_prec[i]], ifelse( i == i_prec[i], Nprec[i], Nsurv_sim[i-1]))
  }

}"


#
# SD with non-constant concentration
#


jags_TKTD_varSD <-
"model {
  ######### priors
  kk_log10 ~ dnorm(kk_meanlog10, kk_taulog10)
  z_log10  ~ dnorm(z_meanlog10 , z_taulog10)
  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  ##### parameter transformation
  kk <- 10**kk_log10
  z  <- 10**z_log10
  kd <- 10**kd_log10
  hb <- 10**hb_log10

  ########## Computation of the likelihood

  for( i in 1:n_dataLong){

    #### midpoint method:
    conc_midpoint[profile_ID_long[i], time_ID_long[i]] <-  (exp(kd * time_long[i]) * conc_long[i] + exp(kd * tprec_long[i])
    * concprec_long[i]) / 2 * (time_long[i] - tprec_long[i])

    D_int[profile_ID_long[i], time_ID_long[i]] <- kd * exp(-kd * time_long[i]) *
    sum( conc_midpoint[profile_ID_long[i], 1:time_ID_long[i]] )

    h[profile_ID_long[i], time_ID_long[i]] = kk * max(0, D_int[profile_ID_long[i], time_ID_long[i]] - z) + hb

    h_midPoint[profile_ID_long[i], time_ID_long[i]] = ((h[profile_ID_long[i], time_ID_long[i]-1] +
            h[profile_ID_long[i], time_ID_long[i]-1])/2) * (time_long[i] - tprec_long[i])

    H_int[profile_ID_long[i], time_ID_long[i]] = sum( h_midPoint[profile_ID_long[i], 1:time_ID_long[i]] )

  }
  for( i in 1:n_dataRed){

     H[profile_ID[i], time_ID[i]]  <- H_int[profile_ID[i], time_ID_red[i]]

     H[i] = sum(h_midPoint[gr, 1:time_interp[gr,i+1]])

     psurv[i] = exp( - H[i])

     Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

    ## ---------------------- generated data

    Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

  }

  ### initialization is requires to use 'Nsurv_sim[i-1]' require in JAGS language (avoid auto-loop issue).
  Nsurv_sim[1] ~ dbin(psurv[1]/psurv[1], Nprec[1])
  for( i in 2:n_dataRed){
    Nsurv_sim[i] ~ dbin(psurv[i]/psurv[i_prec[i]], ifelse( i == i_prec[i], Nprec[i], Nsurv_sim[i-1]))
  }

}"


#
# IT model with log-normal function.
#

jags_TKTD_varIT <-"model {
  #------------------------------------------ parameter transformation
  kd_taulog10 <- 1 / kd_sdlog10^2
  hb_taulog10 <- 1 / hb_sdlog10^2
  alpha_taulog10 <- 1 / alpha_sdlog10^2

  #-------------------------------------------------------- priors
  kd_log10 ~ dnorm(kd_meanlog10, kd_taulog10)
  hb_log10 ~ dnorm(hb_meanlog10, hb_taulog10)

  alpha_log10 ~ dnorm(alpha_meanlog10, alpha_taulog10)
  beta_log10 ~ dunif(beta_minlog10, beta_maxlog10)

  #------------------------------------------ parameter transformation

  kd <- 10**kd_log10
  hb <- 10**hb_log10

  alpha <- 10**alpha_log10
  beta <- 10**beta_log10

  ##------------------------------------ model

  for( i in 1:n_dataLong){

    # trapezoid method :
    diff.int[profile_ID_long[i], time_ID_long[i]] <-  (exp(kd * time_long[i]) * conc_long[i] + exp(kd * tprec_long[i])
    * concprec_long[i]) / 2 * (time_long[i] - tprec_long[i])

    D_int[profile_ID_long[i], time_ID_long[i]] <- kd * exp(-kd * time_long[i]) *
    sum( diff.int[profile_ID_long[i], 1:time_ID_long[i]] )

  }

  for( i in 1:n_dataRed){

    D[profile_ID[i], time_ID[i]]  <- D_int[profile_ID[i], time_ID_red[i]]

    D_max[profile_ID[i], time_ID[i]] <- max(D[profile_ID[i],1:time_ID[i]])

    F[i]  <- D_max[profile_ID[i], time_ID[i]]**beta / ( D_max[profile_ID[i], time_ID[i]]**beta + alpha**beta )

    psurv[i] <-  exp(-hb * time[i]) * (1- F[i])

    Nsurv[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

    ## ---------------------- generated data

    Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i_prec[i]] , Nprec[i])

  }

  ### initialization is requires to use 'Nsurv_sim[i-1]' require in JAGS language (avoid auto-loop issue).
  Nsurv_sim[1] ~ dbin(psurv[1]/psurv[1], Nprec[1])
  for( i in 2:n_dataRed){
    Nsurv_sim[i] ~ dbin(psurv[i]/psurv[i_prec[i]], ifelse( i == i_prec[i], Nprec[i], Nsurv_sim[i-1]))
  }

}"

