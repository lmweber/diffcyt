context("Prepare data")
library(SummarizedExperiment)

test_that("'d_se' object has correct structure", {
  # note: using example data from 'helper_create_example_data.R'
  
  # Prepare data
  d_se <- prepareData(d_input, experiment_info, marker_info)
  
  expect_is(d_se, "SummarizedExperiment")
  
  expect_equal(length(assays(d_se)), 1)
  expect_match(names(assays(d_se)), "exprs")
  
  expect_equal(nrow(d_se), 4000)
  expect_equal(ncol(d_se), 20)
  
  expect_equal(colnames(rowData(d_se)), c("sample_id", "group_id"))
  expect_equal(colnames(colData(d_se)), c("channel_name", "marker_name", "marker_class"))
  
  expect_equal(names(metadata(d_se)), c("experiment_info", "n_cells"))
})
