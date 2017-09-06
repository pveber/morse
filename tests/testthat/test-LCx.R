# For constant concentration exposure

load(file = "tests/testdata/fit_cstSD.rda")
load(file = "tests/testdata/fit_cstIT.rda")

test_that("LCx_survFitCstExp", {
  
  LCx_cstIT <- LCx(fit_cstIT, X = 50)
  LCx_cstSD <- LCx(fit_cstSD, X = 50)

  expect_is(LCx_cstIT, "list")
  expect_is(LCx_cstSD, "list")

})

# For time varying concentration

load(file = "tests/testdata/fit_varSD.rda")
load(file = "tests/testdata/fit_varIT.rda")

test_that("LCx_survFitVarExp", {
  
  LCx_varIT <- LCx(fit_varIT, X = 50)
  LCx_varSD <- LCx(fit_varSD, X = 50)
  
  expect_is(LCx_varIT, "list")
  expect_is(LCx_varSD, "list")
  
})


