tmp_formula <- formula(y~Surv(X,I)+z+(1|r))
data_sim <- simulate_data(
  n = 20,
  formula = tmp_formula,
  n_levels_fixeff = 2,
  # n_levels_raneff = 20,
  type = "glmer",
  b = list(b=c(-5,-2,0.2)),
  # weibull_params = list(X = list(shape = 0.75, scale = 0.25),
  #                       C = list(shape = 0.5, scale = 0.5)),
  censoring_dependent_on_covariate = TRUE,
  error_variance = 0,
  # variance_fixeff = 0.1,
  # variance_raneff = 0.1,
  number_of_clusters = 10,
  number_of_differential_clusters = 1)

d_counts <- data_sim[["d_counts"]]
data_sim <- data_sim[["out"]]
data_sim$z <-
  factor(data_sim$z,labels = 0:(length(unique(data_sim$z))-1))
data_sim$r <-
  factor(data_sim$r,labels = 0:(length(unique(data_sim$r))-1))
da_formula <- list(formula = tmp_formula,
                   data = dplyr::select(data_sim,"X","I","z","r"),
                   random_terms = TRUE)
contrast <- diffcyt::createContrast(c(0, 1, 0))

outs <- testDA_censoredGLMM(d_counts = d_counts, formula = da_formula,
                            contrast = contrast, method_est = "cc",
                            verbose = FALSE, m = 10)

test_that("class testDA_censoredGLMM correct",{
  expect_equal(class(outs)[1],"SummarizedExperiment")
  expect_equal(dim(SummarizedExperiment::rowData(outs)),c(10,3))
  expect_equal(dim(SummarizedExperiment::assay(outs)),c(10,20))
})

test_that(" testDA_censoredGLMM valid pvalues",{
  expect_true(all(SummarizedExperiment::rowData(outs)[ ,2] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,2] <= 1) )
  expect_true(all(SummarizedExperiment::rowData(outs)[ ,3] >= 0) & all(SummarizedExperiment::rowData(outs)[ ,3] <= 1) )
})
