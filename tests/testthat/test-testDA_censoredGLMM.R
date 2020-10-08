set.seed(123)
imputation_methods <- c("cc","pmm","rs","km","km_exp","km_wei","km_os","mrl")
# tmp_formula <- formula(y~Surv(X,I)+z+(1|r))
# create small data set with 2 differential clusters with 10 samples.
d_counts <- simulate_multicluster(alphas = runif(10,1e4,1e5),
                                  sizes = runif(10,1e4,1e5),
                                  nr_diff = 2,
                                  group=2,
                                  return_summarized_experiment = TRUE)$counts
experiment_info <- SummarizedExperiment::colData(d_counts)
experiment_info$status <- sample(c(0,1),size=10,replace = TRUE,prob = c(0.3,0.7))
experiment_info$covariate[experiment_info$status == 0] <-
  runif(10-sum(experiment_info$status),
        min=0,
        max=experiment_info$covariate[experiment_info$status == 0])
da_formula <- createFormula(experiment_info,
                            cols_fixed = c("covariate", "group_covariate"),
                            cols_random = "sample",event_indicator = "status")
contrast <- createContrast(c(0, 1, 0))
# outs <- testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
#                             contrast = contrast, mi_reps = 2, imputation_method = "km")

test_testDA_censoredGLMM <- function(method){
  outs <- suppressWarnings(suppressMessages(testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
                              contrast = contrast, imputation_method = method,
                              verbose = FALSE, mi_reps = 2, BPPARAM = BiocParallel::SerialParam())))
  
  test_that(paste("class testDA_censoredGLMM correct for",method),{
    expect_true(is(outs, "SummarizedExperiment"))
    expect_equal(dim(SummarizedExperiment::rowData(outs)),c(10,7))
    expect_equal(dim(SummarizedExperiment::assay(outs)),c(10,10))
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
