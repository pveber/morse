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
  interpolate_length = NULL
  interpolate_method = "linear"
  
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
  kd=0:2
  hb=0:2
  z=0:2
  kk=0:2
  alpha=0:2
  beta=1:3
  mcmc_size=length(kd)
  interpolate_length = NULL
  interpolate_method = "constant"
  
  morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size)
  morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size)          
  
  # check No ERROR
  expect_error(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_error(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  # check No WARNING
  expect_warning(morse:::SurvIT_ode(Cw, time, replicate, kd, hb, alpha, beta, mcmc_size = mcmc_size), NA)
  expect_warning(morse:::SurvSD_ode(Cw, time, replicate, kd, hb, z, kk, mcmc_size = mcmc_size), NA)
  
})




