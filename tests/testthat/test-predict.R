test_that("Lag of column in matrix", {
  
  lambda = matrix(1:12, nrow = 3)
  expect_equal(cbind(lambda[,1], lambda[,1:(ncol(lambda)-1)]),
               matrix(c(1:3,1:9), nrow = 3) )
  
})


test_that("Surv.SD_Cext", {
  
  Cw=5:2
  time=1:4
  kk=0:2
  kd=0:2
  z=0:2
  hb=0:2
  
  # Here we build a line by line function to test the matrix-vector function
  Surv.SD_Cext_nonVectorize_ForTesting = function(Cw, time, kk, kd, z, hb){
    time.prec = c(time[1], time[1:(length(time)-1)])
    diff.int = (exp(time * kd) + exp(time.prec * kd) )*Cw/2 * (time-time.prec) 
    D = kd * exp(-kd * time) * cumsum(diff.int)
    lambda = kk * pmax(D-z,0) + hb 
    lambda.prec = c(lambda[1], lambda[1:(length(lambda)-1)])
    int.lambda =  (lambda + lambda.prec)/2 * (time-time.prec)
    return( exp(-cumsum(int.lambda)) )
  }
  
  matTest = matrix(NA,nrow = length(kk), ncol = length(time))
  for(i in 1:length(kk)){
    matTest[i,] =  Surv.SD_Cext_nonVectorize_ForTesting(Cw, time, kk[i], kd[i], z[i], hb[i])
  }
  
  expect_equal( Surv.SD_Cext(Cw, time, kk, kd, z, hb),  matTest,  tolerance=1e-5)
})

test_that("MCMC length one work", {

  Cw=1:2
  time=1:2
  replicate="A"
  kd=0.5
  hb=0.2
  z=16
  kk=1.9
  alpha=16
  beta=1.9
  mcmc_size=length(kd)

  # check No ERROR
  expect_error(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_error(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  # check No WARNING
  expect_warning(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_warning(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  
})

test_that("MCMC longer than one", {
  
  Cw=5:2
  time=1:4
  replicate="A"
  kd=c(0.5,2,1.3)
  hb=c(0.5,2,1.3)
  z=c(16,5,2)
  kk=c(0.5,2,1.3)
  alpha=c(16,5,2)
  beta=c(0.5,2,1.3)
  mcmc_size=length(kd)

  # check No ERROR
  expect_error(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_error(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  # check No WARNING
  expect_warning(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_warning(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  
})


test_that("predict_ode", {
  
  skip_on_cran()
  
  data("propiconazole")
  fit_cstSD <- survFit(survData(propiconazole), quiet = TRUE, model_type = "SD")
  
  data_4prediction <- data.frame(time = c(1:10, 1:10),
                                 conc = c(c(0,0,40,0,0,0,40,0,0,0),
                                          c(21,19,18,23,20,14,25,8,13,5)),
                                 replicate = c(rep("pulse", 10), rep("random", 10)))
  
  # check No ERROR
  expect_error(predict_ode(object = fit_cstSD, data_predict = data_4prediction), NA)


  data_4MFx <- data.frame(time = 1:10,
                          conc = c(0,0.5,8,3,0,0,0.5,8,3.5,0))
  # check No ERROR
  expect_error(MFx(object = fit_cstSD, data_predict = data_4MFx, ode = TRUE), NA)

})


test_that("predict_Nsurv_ode internal", {
  
  skip_on_cran()
  
  data("propiconazole")
  fit_cstSD <- survFit(survData(propiconazole), quiet = TRUE, model_type = "SD")
  fit_cstIT <- survFit(survData(propiconazole), quiet = TRUE, model_type = "IT")
  
  data("FOCUSprofile")  
  FOCUSprofile$Nsurv = sort(round(runif(nrow(FOCUSprofile), 0, 100)), decreasing = TRUE)
  
  # check No ERROR
  expect_error(predict_Nsurv_ode(object = fit_cstSD, data_predict = FOCUSprofile, mcmc_size = 10, interpolate_length  = NULL), NA)
  expect_error(predict_Nsurv_ode(object = fit_cstIT, data_predict = FOCUSprofile, mcmc_size = 10, interpolate_length  = NULL), NA)
  
  expect_error(predict_Nsurv_ode(object = fit_cstSD, data_predict = NULL, mcmc_size = 1000, interpolate_length = 10), NA)
  expect_error(predict_Nsurv(object = fit_cstSD, data_predict = NULL, mcmc_size = 1000, interpolate_length = 10), NA)

  data("propiconazole_pulse_exposure")
  expect_error(predict_Nsurv_ode(fit_cstSD, propiconazole_pulse_exposure, mcmc_size = NULL, interpolate_length = NULL), NA)
  expect_error(predict_Nsurv_ode(fit_cstSD, propiconazole_pulse_exposure, mcmc_size = NULL, interpolate_length = NULL), NA)

})

test_that("predict_interpolate", {
  
  data(FOCUSprofile)
  FC_predInterp = morse:::predict_interpolate(FOCUSprofile,  extend_time = 100)
  expect_true(nrow(FC_predInterp) >= nrow(FOCUSprofile))
  
})
