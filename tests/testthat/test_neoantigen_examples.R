
library(testthat)
require("NeoantigenProfiler") || stop("unable to load NeoantigenProfiler package")
#GenomicRanges:::.test()

version <- 1

test_that("Correct Version", {
  expect_equal(version, 1)
})
