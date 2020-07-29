
imputation_methods <- c("cc","pmm","rs","km","km_exp","km_wei","km_os","mrl")
tmp_formula <- formula(y~Surv(X,I)+z+(1|r))
data_sim <- simulate_singlecluster(
  n = 20,
  formula = tmp_formula,
  n_levels_fixeff = 2,
  type = "glmer",
  b = list(b=c(-4,-0.5,0.5)),
  censoring_dependent_on_covariate = FALSE,
  error_variance = 0,
  variance_fixeff = 1,
  variance_raneff = 1,
  number_of_clusters = 3,
  number_of_differential_clusters = 1,
  transform_fn = "log_positive"
  )
d_counts <- data_sim
experiment_info <- SummarizedExperiment::colData(data_sim)
da_formula <- createFormula(experiment_info,cols_fixed = c("X","z"), cols_random = "r", event_indicator = "I")
contrast <- diffcyt::createContrast(c(0, 1, 0))

test_testDA_censoredGLMM <- function(method){
  outs <- testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
                              contrast = contrast, imputation_method = method,
                              verbose = FALSE, mi_reps = 2, BPPARAM = BiocParallel::SerialParam())
  
  test_that(paste("class testDA_censoredGLMM correct for",method),{
    expect_true(is(outs, "SummarizedExperiment"))
    expect_equal(dim(SummarizedExperiment::rowData(outs)),c(3,4))
    expect_equal(dim(SummarizedExperiment::assay(outs)),c(3,20))
  })
  
  test_that(paste("testDA_censoredGLMM keeps entries",method),{
    expect_equal(SummarizedExperiment::assay(outs),SummarizedExperiment::assay(d_counts))
    expect_equal(sort(SummarizedExperiment::rowData(outs)[1]),sort(SummarizedExperiment::rowData(d_counts)[1]))
    expect_equal(sort(SummarizedExperiment::rowData(outs)[4]),sort(SummarizedExperiment::rowData(d_counts)[2]))
  })
  
  test_that(paste("testDA_censoredGLMM valid pvalues",method),{
    expect_true(all(SummarizedExperiment::rowData(outs)[ ,2] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,2] <= 1) )
    expect_true(all(SummarizedExperiment::rowData(outs)[ ,3] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,3] <= 1) )
  })
}

for(method in imputation_methods){
  test_testDA_censoredGLMM(method)
}
