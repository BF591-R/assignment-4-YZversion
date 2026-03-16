#!/usr/bin/Rscript
source("main.R")
library(testthat)

#generate consistent sample data for all tests using it
set.seed(42)
fx_intensity_df <- data.frame(GSM971958 = runif(1000, 0, 15),
                              GSM971959 = runif(1000, 0, 15),
                              GSM971967 = runif(1000, 0, 15),
                              GSM971968 = runif(1000, 0, 15))

row.names(fx_intensity_df) <- paste0(1:1000, "_s_at")
fx_intensity_mat <- as.matrix(fx_intensity_df)

describe("read_data()", {
  res <- read_data("data/example_intensity_data.csv", " ")
  
  it("loads data as a dataframe, not tibble", {
    expect_true(is.data.frame(res))
  })
  it("returns the correct number of dimensions", {
    expect_equal(dim(res), c(54675, 35))
  })
})

describe("calculate_variance()", {
  test_obj <- list(sdev = c(8,4,2))
  test_res <- c(0.76190476, 0.19047619, .04761905)
  
  it("correctly accesses the sdev values, squares them, and divides each by the sum", {
    expect_true(all(dplyr::near(calculate_variance_explained(test_obj), test_res, tol=0.001)))
  })
})

describe("make_variance_tibble()", {
  test_pca_res <- prcomp(scale(t(fx_intensity_mat)), scale=FALSE, center=FALSE)
  test_ve <- c(1, 1, 1, 1)
  test_tib <- make_variance_tibble(test_ve, test_pca_res)
  expected_cumsum <- cumsum(test_ve)
  
  it("returns a tibble of the concatenated results", {
    expect_true(is_tibble(test_tib))
  })
  it("returns a tibble with the correct dimensions from the input", {
    expect_equal(dim(test_tib), c(length(test_ve), 3))
  })
  it("returns a column entitled cumulative with the cumulative sum", {
    expect_equal(test_tib$cumulative, expected_cumsum)
  })
})

describe("list_significant_probes()", {
  test_tibble <- tibble(
    probeid = paste0(1:4, "_s_at"),
    padj = c(.005, .015, .01, .008)
  )
  
  it("correctly filters out values according to padj < threshold", {
    expect_equal(list_significant_probes(test_tibble,.01), c("1_s_at", "4_s_at"))
  })
})

describe("return_de_intensity()", {
  #tests is.matrix() and since function subsets the matrix should be identical
  #no need to worry about rounding, dimensions, etc. ?
  #
  test_mat <- as.matrix(rbind(fx_intensity_df[1, ], fx_intensity_df[4, ]))
  function_mat <- return_de_intensity(fx_intensity_df, c("1_s_at", "4_s_at"))
  
  it("returns a matrix", {
    expect_true(is.matrix(function_mat))
  })
  it("returns the correct values based on the filtering", {
    expect_true(identical(function_mat, test_mat))
  })
})

