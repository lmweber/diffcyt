censored_variable <-  "X"
estimate <- paste0(censored_variable,"_est")
censoring_indicator <- "I"
covariates <-  c("Z")
response <-  "Y"
methods <- c("rs","km","mrl")

set.seed(123)
conditions <- list(formulas = list(lm_z = formula(Y ~ X + Z),
                                   glm_z = formula(Y ~ X + Z),
                                   glmer_z = formula(Y ~ X + Z + (1|R)),
                                   lm = formula(Y ~ X),
                                   glm = formula(Y ~ X),
                                   glmer = formula(Y ~ X + (1|R))),
                   formulas_surv = list(lm_z = formula(Y ~ Surv(X,I) + Z),
                                        glm_z = formula(Y ~ Surv(X,I) + Z),
                                        glmer_z = formula(Y ~ Surv(X,I) + Z + (1|R)),
                                        lm = formula(Y ~ Surv(X,I) ),
                                        glm = formula(Y ~ Surv(X,I)),
                                        glmer = formula(Y ~ Surv(X,I)  + (1|R))),
                   regression_type = list(lm_z = "lm",
                                          glm_z = "glm",
                                          glmer_z = "glmer",
                                          lm = "lm",
                                          glm = "glm",
                                          glmer = "glmer"))

test_data_sets <- list(lm_z = dplyr::arrange(simulate_data(50, conditions$formulas_surv$lm_z, type = "lm"),ID),
                       glm_z = dplyr::arrange(simulate_data(50, conditions$formulas_surv$glm_z, type = "glm"),ID),
                       glmer_z = dplyr::arrange(simulate_data(50, conditions$formulas_surv$glmer_z, type = "glmer"),ID),
                       lm = dplyr::arrange(simulate_data(50, conditions$formulas_surv$lm, type = "lm"),ID),
                       glm = dplyr::arrange(simulate_data(50, conditions$formulas_surv$glm, type = "glm"),ID),
                       glmer = dplyr::arrange(simulate_data(50, conditions$formulas_surv$glmer, type = "glmer"),ID))


test_for_equal_input_output <- function(output_to_test, output_name, test_data_sets){
  test_that(paste(output_name, "keeps entries"),{
    purrr::map2(output_to_test,
                test_data_sets,
                ~ expect_equal(c(.x[c(.x[,censoring_indicator] == 1),censored_variable])[[1]],
                               c(.y[c(.y[,censoring_indicator] == 1),censored_variable])[[1]]))
    purrr::map2(output_to_test,
                test_data_sets,
                ~ expect_equal(c(.x[ ,response])[[1]], c(.y[ ,response])[[1]]))
    purrr::map2(output_to_test,
                test_data_sets,
                ~ expect_equal(c(.x[ ,censoring_indicator])[[1]], c(.y[ ,censoring_indicator])[[1]]))
  })
}

test_for_valid_estimates <- function(output_to_test, output_name, test_data_sets){
  test_that(paste(output_name,"estimations are valid"),{
    purrr::map2(output_to_test,
                test_data_sets,
                ~ expect_true(all(.x[[estimate]] >= .y[[estimate]])))
  })
}



##### cond.single.imputation :
combinations <- expand.grid(dataset=4:6,imputation_method=methods,stringsAsFactors = FALSE)
combinations$imputation_method <- as.character(combinations$imputation_method)
conditional_single_imputation_result_without_z <- purrr::pmap(
  combinations,
  function(dataset,imputation_method) {
    dplyr::arrange(conditional_single_imputation(
      data = test_data_sets[[dataset]],
      censored_variable = censored_variable,
      censoring_indicator = censoring_indicator,
      response = response, id = "ID",
      imputation_method = imputation_method
      ), ID)})

combinations <- expand.grid(dataset=1:3,imputation_method=methods,stringsAsFactors = FALSE)
combinations$imputation_method <- as.character(combinations$imputation_method)
conditional_single_imputation_result_with_z <- purrr::pmap(
  combinations,
  function(dataset,imputation_method) {
    dplyr::arrange(conditional_single_imputation(
    data = test_data_sets[[dataset]],
    censored_variable = censored_variable,
    censoring_indicator = censoring_indicator,
    response = response,
    covariates = covariates, id = "ID",
    imputation_method = imputation_method
  ), ID)})

test_for_equal_input_output(conditional_single_imputation_result_without_z,"conditional_single_imputation",rep(test_data_sets[4:6],length(methods)))
test_for_equal_input_output(conditional_single_imputation_result_with_z,"conditional_single_imputation",rep(test_data_sets[1:3],length(methods)))

test_that("conditional_single_imputation keeps dimensionality, without z",{
  purrr::map2(conditional_single_imputation_result_without_z,
              rep(test_data_sets[4:6],length(methods)),
              ~ expect_equal(length(.x), length(.y)+1))
})
test_that("conditional_single_imputation keeps dimensionality, with z",{
  purrr::map2(conditional_single_imputation_result_with_z,
              rep(test_data_sets[1:3],length(methods)),
              ~ expect_equal(length(.x), length(.y)+1))
})

test_for_valid_estimates(conditional_single_imputation_result_without_z,"conditional_single_imputation",rep(test_data_sets[4:6],length(methods)))
test_for_valid_estimates(conditional_single_imputation_result_with_z,"conditional_single_imputation",rep(test_data_sets[1:3],length(methods)))



