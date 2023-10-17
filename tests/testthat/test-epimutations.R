source("init_data.R") 

context("epimutations::output")



test_that("epimutation output class", {
  #length
  expect_length(class(result_manova), 3)
  expect_length(class(result_mlm), 3)
  expect_length(class(result_iForest), 3)
  expect_length(class(result_mahdist), 3)
  expect_length(class(result_quantile), 3)
  expect_length(class(result_beta), 3)
  #class
  expect_true(all(class(result_manova) %in% c("tbl_df", "tbl", "data.frame")))
  expect_equal(class(result_mlm), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_iForest), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_mahdist), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_quantile), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_beta), c("tbl_df", "tbl", "data.frame"))
})

test_that("epimutation column names", {
  #length
  expect_length(colnames(result_manova), 16)
  expect_length(colnames(result_mlm), 16)
  expect_length(colnames(result_iForest), 16)
  expect_length(colnames(result_mahdist), 16)
  expect_length(colnames(result_quantile), 16)
  expect_length(colnames(result_beta), 16)
  #names
  expect_true(all(colnames(result_manova) %in% col_names))
  expect_true(all(colnames(result_mlm) %in% col_names))
  expect_true(all(colnames(result_iForest) %in% col_names))
  expect_true(all(colnames(result_mahdist) %in% col_names))
  expect_true(all(colnames(result_quantile) %in% col_names))
  expect_true(all(colnames(result_beta) %in% col_names))
})