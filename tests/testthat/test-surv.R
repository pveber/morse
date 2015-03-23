## data
data(cadmium1)
data(cadmium2)
data(copper)
data(chlordan)
data(zinc)

# no error dataset
d <- list(cadmium1 = cadmium1,
          cadmium2 = cadmium2,
          copper = copper,
          chlordan = chlordan,
          zinc = zinc)

# error dataset

zinc0 <- as.list(zinc)

zinc1 <- zinc
colnames(zinc1) <- c("replica", "con", "time", "Nsur", "Nrepro")

zinc2 <- zinc
zinc2[46, "time"] <- 1
zinc2$time <- as.integer(zinc2$time)

zinc3 <- zinc
zinc3$conc <- as.character(zinc3$conc)

zinc4 <- zinc
zinc4$Nsurv <- as.numeric(zinc4$Nsurv)

zinc5 <- zinc
zinc5[69, "Nsurv"] <- -248
zinc5$Nsurv <- as.integer(zinc5$Nsurv)

zinc6 <- zinc
zinc6[1, "Nsurv"] <- 0
zinc6$Nsurv <- as.integer(zinc6$Nsurv)

zinc7 <- zinc
zinc7[107, "replicate"] <- "A"

zinc8 <- zinc
zinc8[25, "Nsurv"] <- 20
zinc8$Nsurv <- as.integer(zinc8$Nsurv)

zinc9 <- zinc
zinc9[, "replicate"] <- as.character(zinc9[, "replicate"])
zinc9[12, "replicate"] <- "D"
zinc9[, "replicate"] <- as.factor(zinc9[, "replicate"])

cadmium19 <- cadmium1
cadmium19[12, "replicate"] <- 5

d2 <- list(zinc0 = zinc0,
           zinc1 = zinc1,
           zinc2 = zinc2,
           zinc3 = zinc3,
           zinc4 = zinc4,
           zinc5 = zinc5,
           zinc6 = zinc6,
           zinc7 = zinc7,
           zinc8 = zinc8,
           zinc9 = zinc9,
           cadmium19 = cadmium19)

## tests

test_that("survDataCheck", {
  skip_on_cran()
  
  # no error dataset
  lapply(d, function(x) {
    dat <- survDataCheck(x)
    expect_equal(dim(dat)[1], 0)
    expect_equal(dim(dat)[2], 0)
    expect_is(dat, c("errorTable",
                     "data.frame"))
    expect_null(dat$id)
    expect_null(dat$msg)
    })
  
  # error dataset
  expect_named(survDataCheck(d2[["zinc0"]],
                             diagnosis.plot = FALSE), c("id", "msg"))
  
  expect_equal(survDataCheck(d2[["zinc0"]],
                             diagnosis.plot = FALSE)$id,
               "dataframeExpected")
  
  expect_equal(survDataCheck(d2[["zinc1"]],
                              diagnosis.plot = FALSE)$id,
               rep("missingColumn", 3))
  
  expect_equal(survDataCheck(d2[["zinc2"]],
                              diagnosis.plot = FALSE)$id[1],
               "firstTime0")
  
  expect_equal(survDataCheck(d2[["zinc3"]], diagnosis.plot = FALSE)$id,
               "concNumeric")
  
  expect_equal(reproDataCheck(d2[["zinc4"]], diagnosis.plot = FALSE)$id,
               "NsurvInteger")
  
  expect_equal(reproDataCheck(d2[["zinc5"]], diagnosis.plot = FALSE)$id[1],
               "tablePositive")
  
  expect_equal(survDataCheck(d2[["zinc6"]], diagnosis.plot = FALSE)$id[1],
               "Nsurv0T0")
  
  expect_equal(survDataCheck(d2[["zinc7"]], diagnosis.plot = FALSE)$id[1:2],
               c("duplicatedID", "missingReplicate"))
  
  expect_equal(survDataCheck(d2[["zinc8"]], diagnosis.plot = FALSE)$id,
               "NsurvIncrease")
  
  expect_equal(survDataCheck(d2[["zinc9"]], diagnosis.plot = FALSE)$id[4],
               "ReplicateLabel")
  
  expect_equal(survDataCheck(d2[["cadmium19"]], diagnosis.plot = FALSE)$id[4],
               "ReplicateLabel")
})