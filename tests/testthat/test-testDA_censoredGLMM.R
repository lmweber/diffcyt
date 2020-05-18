tmp_formula <- formula(y~Surv(X,I)+z+(1|r))
data_sim <- simulate_data(
  n = 20,
  formula = tmp_formula,
  n_levels_fixeff = 2,
  # n_levels_raneff = 20,
  type = "glmer",
  b = list(b=c(-4,-0.5,0.5)),
  # weibull_params = list(X = list(shape = 0.75, scale = 0.5),
  #                       C = list(shape = 0.5, scale = 0.75)),
  censoring_dependent_on_covariate = FALSE,
  error_variance = 0,
  variance_fixeff = 1,
  variance_raneff = 1,
  number_of_clusters = 10,
  number_of_differential_clusters = 2,
  transform_fn = "log_positive"
  )
d_counts <- data_sim
experiment_info <- SummarizedExperiment::colData(data_sim)
# data_sim$z <-
#   factor(data_sim$z,labels = 0:(length(unique(data_sim$z))-1))
# data_sim$r <-
#   factor(data_sim$r,labels = 0:(length(unique(data_sim$r))-1))
da_formula <- createFormula(experiment_info,cols_fixed = c("X","z"), cols_random = "r", event_indicator = "I")
# da_formula <- list(formula = tmp_formula,
#                    data = dplyr::select(experiment_info,"X","I","z","r"),
#                    random_terms = TRUE)
contrast <- diffcyt::createContrast(c(0, 1, 0))

outs <- testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
                            contrast = contrast, imputation_method = "rs",
                            verbose = FALSE, mi_reps = 10, BPPARAM = BiocParallel::MulticoreParam(workers=14))
SummarizedExperiment::rowData(outs)

test_that("class testDA_censoredGLMM correct",{
  expect_true(is(outs, "SummarizedExperiment"))
  expect_equal(dim(SummarizedExperiment::rowData(outs)),c(20,4))
  expect_equal(dim(SummarizedExperiment::assay(outs)),c(20,20))
})

test_that("testDA_censoredGLMM keeps entries",{
  expect_equal(SummarizedExperiment::assay(outs),SummarizedExperiment::assay(d_counts))
  expect_equal(SummarizedExperiment::rowData(outs)[1],SummarizedExperiment::rowData(d_counts)[1])
  expect_equal(SummarizedExperiment::rowData(outs)[4],SummarizedExperiment::rowData(d_counts)[2])
})

test_that(" testDA_censoredGLMM valid pvalues",{
  expect_true(all(SummarizedExperiment::rowData(outs)[ ,2] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,2] <= 1) )
  expect_true(all(SummarizedExperiment::rowData(outs)[ ,3] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,3] <= 1) )
})
