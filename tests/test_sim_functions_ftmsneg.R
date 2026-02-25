library(testthat)
source("C:/Users/amentag/Desktop/these/These/clones/RDynlib/RDynLib/Similarity/sim_functions_with_options.R")

a <- cbind(mz = c(13.3, 15.4, 15.6, 23.4, 26.2),
           intensity = 1:5)
b <- cbind(mz = c(15.3, 33.4),
           intensity = 1:2)

test_that("remove_duplicates works", {
  res <- remove_duplicates(a[, 1], a[, 2])
  expect_true(is.list(res))
  expect_equal(names(res), c("mz", "intensity"))
  ## Here I'm not sure. Would maybe have preferred that 15.4 and 15.6 get
  ## merged into the same? i.e. result would be rather:
  ## c(13, 15, 23, 26), c(1, 3, 4, 5)?
  expect_equal(res$mz, c(13, 15, 16, 23, 26))
  expect_equal(res$intensity, c(1, 2, 3, 4, 5))
})

test_that("dynlib_map works", {
  
  ## Empty input. Expect empty result, but in a specific format
  res <- dynlib_map(matrix(numeric(), nrow = 0), matrix(numeric(), nrow = 0))
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  expect_true(nrow(res$x) == 0)
  expect_true(nrow(res$y) == 0)
  
  ## Only direct matches, no neutral loss matches
  res <- dynlib_map(a, b, 12, 12)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = 15, intensity = 2))
  expect_equal(y, cbind(mz = 15, intensity = 1))
  
  ## Direct and neutral loss matches
  res <- dynlib_map(a, b, 0, 10)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = c(15, 23), intensity = c(2, 4)))
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(y, cbind(mz = c(15, 33), intensity = c(1, 2)))
  
  ## No matches
  d <- cbind(mz = c(33.3, 46.1), intensity = c(1, 2))
  res <- dynlib_map(a, d, 0, 0)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = numeric(), intensity = numeric()))
  expect_equal(y, cbind(mz = numeric(), intensity = numeric()))
  
  ## Only neutral loss matches
  res <- dynlib_map(a, d, 0, 20)
  expect_true(is.list(res))
  expect_equal(names(res), c("x", "y"))
  x <- res$x
  expect_true(attributes(x)$wintensity_sum > 0)
  attributes(x)$wintensity_sum <- NULL
  y <- res$y
  expect_true(attributes(y)$wintensity_sum > 0)
  attributes(y)$wintensity_sum <- NULL
  expect_equal(x, cbind(mz = c(13, 26), intensity = c(1, 5)))
  expect_equal(y, cbind(mz = c(33, 46), intensity = c(1, 2)))
})





  
  

  
  
  
 
  
  
  
  
  
