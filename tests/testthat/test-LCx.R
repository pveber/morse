test_that("test LCx", {
  
  skip_on_cran()
  
  data("propiconazole")
  fit_cstSD <- survFit(survData(propiconazole), quiet = TRUE, model_type = "SD")
  LCx_cstSD <- LCx(fit_cstSD, X = 50)
  plot(LCx_cstSD)

  LCx_cstSD <- LCx(fit_cstSD, X = 0.00001)
  plot(LCx_cstSD)
})

