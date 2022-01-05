source("init_data.R") 

context("epimutations::output")



test_that("epimutation output class", {
  #length
  expect_length(class(result_manova), 3)
  expect_length(class(result_mlm), 3)
  expect_length(class(result_isoforest), 3)
  expect_length(class(result_mahadistmcd), 3)
  expect_length(class(result_quantile), 3)
  expect_length(class(result_beta), 3)
  #class
  expect_true(all(class(result_manova) %in% c("tbl_df", "tbl", "data.frame")))
  expect_equal(class(result_mlm), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_isoforest), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_mahadistmcd), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_quantile), c("tbl_df", "tbl", "data.frame"))
  expect_equal(class(result_beta), c("tbl_df", "tbl", "data.frame"))
})

test_that("epimutation column names", {
  #length
  expect_length(colnames(result_manova), 15)
  expect_length(colnames(result_mlm), 15)
  expect_length(colnames(result_isoforest), 15)
  expect_length(colnames(result_mahadistmcd), 15)
  expect_length(colnames(result_quantile), 15)
  expect_length(colnames(result_beta), 15)
  #names
  expect_true(all(colnames(result_manova) %in% col_names))
  expect_true(all(colnames(result_mlm) %in% col_names))
  expect_true(all(colnames(result_isoforest) %in% col_names))
  expect_true(all(colnames(result_mahadistmcd) %in% col_names))
  expect_true(all(colnames(result_quantile) %in% col_names))
  expect_true(all(colnames(result_beta) %in% col_names))
})