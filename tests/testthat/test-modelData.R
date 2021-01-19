test_that("modelData", {
  
  skip_on_cran()
  
  data("propiconazole")
  data("propiconazole_pulse_exposure")
  
  df_MD_Cst <- morse:::addVariable_survDataCstExp(survData(propiconazole))
  expect_length(unique(df_MD_Cst$replicate), 8)
  expect_length(unique(df_MD_Cst$replicate_ID), 8)

  df_MD_Var <- morse:::survData_interpolate(survData(propiconazole_pulse_exposure))
  expect_length(unique(df_MD_Var$replicate), 4)
  expect_length(unique(df_MD_Var$replicate_ID), 4)
  
  # time points of the profile are within the time extended
  expect_true(all(propiconazole_pulse_exposure$time %in% df_MD_Var$time))
  
})
