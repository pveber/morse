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

    ifelse(time[i] > 1,
           Nsurv_sim[i] ~ dbin(psurv[i] , Nsurv[i-1]),
           Nsurv_sim[i] ~ dbin(psurv[i] , Nprec[i]))
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

    ifelse(time[i] > 1,
           Nsurv[i] ~ dbin(psurv[i]/psurv[i-1] , Nprec[i]),
           Nsurv[i] ~ dbin(psurv[i]/1 , Nprec[i]) )

    ## ---------------------- generated data

    ifelse(time[i] > 1,
           Nsurv_ppc[i] ~ dbin(psurv[i]/psurv[i-1] , Nprec[i]),
           Nsurv_ppc[i] ~ dbin(psurv[i]/1 , Nprec[i]) )

    ifelse(time[i] > 1,
           Nsurv_sim[i] ~ dbin(psurv[i]/psurv[i-1] , Nsurv_sim[i-1]),
           Nsurv_sim[i] ~ dbin(psurv[i]/1 , Nprec[i])

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

  ########## Computation of the likdlihood

  for (gr in 1:n_data){

    #---- Integration for the internal concentration
    diff.int[gr,1]=0
    int.hazard[gr,1]=0

    for(i in 2:N.intC[gr]){

      ### midpoint method:
      conc_midPoint[gr,i] = ((exp(kd * time_interp[gr,i]) * conc_interp[gr,i] + exp(kd * time_interp[gr,i-1]) * conc_interp[gr,i-1]) / 2 ) * (time_interp[gr,i] - time_interp[gr,i-1])

    }
    for(i in 1:N_int[gr]){
      D[gr,i] = kd * exp(- kd * time_interp[gr,i]) * sum(conc_midPoint[gr,1:i])

      hazard[gr,i] = kk * max(0, D[gr,i] - z) + hb
    }
    for(i in 2:N.intC[gr]){

      #### midpoint method:
      hazard_midPoint[gr,i] = ((hazard[gr,i] + hazard[gr,i-1])/2) * (time_interp[gr,i] - time_interp[gr,i-1])

    }
    #---------------------------------------------

    hazard_integr[gr,1] = -sum(hazard_midPoint[gr, 1:id_time_interp[gr,2]])

    psurv[gr,1] = exp(hazard_integr[gr,1])

    Nsurv[gr,1] ~ dbin(psurv[gr,1] , Nsurv_prec[gr,1])

    for (t in 2:n.time[gr]){

      hazard_integr[gr,t] = -sum( hazard_midPoint[ gr, 1:id.intCtime[gr,t+1] ] )

      psurv[gr,t] = exp(hazard_integr[gr,t])

      Nsurv[gr,t] ~ dbin(psurv[gr,t] / psurv[gr,t-1] , Nsurv_prec[gr,t])

    }

  }
}"


#
# IT model with log-normal function.
#


