data("cadmium2")

survDataset <- survData(cadmium2)

survfitDataset_SD <- survFit(survDataset, model_type = "SD")
survfitDataset_IT <- survFit(survDataset, model_type = "IT")


survfitDataset_IT

test_that("LCx_survFitCstExp", {
  
  LCx_IT <- LCx(survfitDataset_IT, X = 50)
  LCx_SD <- LCx(survfitDataset_SD, X = 50)

  expect_is(LCx_IT, "list")
  expect_is(LCx_SD, "list")

})




