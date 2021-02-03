test_that("survData replicate is factor", {
  
  data("propiconazole")
  class(propiconazole$replicate)
  expect_s3_class(survData(propiconazole), "data.frame")
  
  expect_s3_class(survData(propiconazole)$replicate, "factor")
  propiconazole$replicate = as.character(propiconazole$replicate)
  expect_type(propiconazole$replicate, "character")
  
  expect_s3_class(survData(propiconazole)$replicate, "factor")
  expect_type(survData(propiconazole)$replicate, "integer")
  
  data("cadmium1")
  class(cadmium1$replicate)
  expect_s3_class(survData(cadmium1)$replicate, "factor")
  cadmium1$replicate = as.character(cadmium1$replicate)
  expect_type(cadmium1$replicate, "character")
  
  expect_s3_class(survData(cadmium1)$replicate, "factor")
  expect_type(survData(cadmium1)$replicate, "integer")

})
