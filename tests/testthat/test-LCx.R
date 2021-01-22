test_that("test LCx", {
  
  skip_on_cran()
  
  data("propiconazole")
  fit_cstSD <- survFit(survData(propiconazole), quiet = TRUE, model_type = "SD")
  
  expect_error(LCx(fit_cstSD, X = 50), NA)
  
  LCx_cstSD <- LCx(fit_cstSD, X = 50)
  expect_is( plot(LCx_cstSD), "ggplot")

  expect_error( LCx(fit_cstSD, X = 0.00001), NA)
  
  LCx_cstSD <- LCx(fit_cstSD, X = 0.00001)
  expect_is( plot(LCx_cstSD), "ggplot")
})

