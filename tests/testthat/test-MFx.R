test_that("test LCx", {
  
  skip_on_cran()
  
  data("propiconazole")
  fit_cstSD <- survFit(survData(propiconazole), quiet = TRUE, model_type = "SD")
  data_4MFx <- data.frame(time = 1:10,
                          conc = c(0,0.5,8,3,0,0,0.5,8,3.5,0))
  expect_error(MFx(object = fit_cstSD, data_predict = data_4MFx, ode = TRUE), NA)
  
  MFx_PRZ_cstSD <- MFx(object = fit_cstSD, data_predict = data_4MFx, ode = TRUE)
  expect_is( plot(MFx_PRZ_cstSD), "ggplot")

})


