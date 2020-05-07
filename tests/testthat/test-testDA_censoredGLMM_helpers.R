test_formulas <- list(formula(y~Surv(x,i)),
                      formula(y~Surv(x,i)+z),
                      formula(r~Surv(c,p)+q),
                      formula(y~Surv(x,i)+z+zz),
                      formula(y~Surv(x,i)+z+zz+zzz),
                      formula(y~Surv(x,i)+(1|r)),
                      formula(y~Surv(x,i)+z+(1|r)),
                      formula(y~Surv(x,i)+z+zz+(1|r)))
test_formulas_glmm <- list(formula(y~x),
                      formula(y~x+z),
                      formula(r~c+q),
                      formula(y~x+z+zz),
                      formula(y~x+z+zz+zzz),
                      formula(y~x+(1|r)),
                      formula(y~x+z+(1|r)),
                      formula(y~x+z+zz+(1|r)))

expected <- list(list(censored_variable="x",censoring_indicator="i",response="y",covariates=NULL,random_covariates=NULL),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates="z",random_covariates=NULL),
                 list(censored_variable="c",censoring_indicator="p",response="r",covariates="q",random_covariates=NULL),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates=c("z","zz"),random_covariates=NULL),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates=c("z","zz","zzz"),random_covariates=NULL),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates=NULL,random_covariates="r"),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates="z",random_covariates="r"),
                 list(censored_variable="x",censoring_indicator="i",response="y",covariates=c("z","zz"),random_covariates="r"))

outs <- purrr::map(test_formulas, ~extract_variables_from_formula(.x))
purrr::walk(1:8 , function(i){
  test_that(paste("extract_variables_from_formula works for formula:",paste(as.character(test_formulas[[i]]),collapse = " ")), {
    expect_equal(!!outs[[i]], !!expected[[i]])
  })
})

outs <- purrr::map(test_formulas, ~create_glmm_formula(.x))
purrr::walk(1:8 , function(i){
  test_that(paste("create_glmm_formula works for formula:",paste(as.character(test_formulas[[i]]),collapse = " ")), {
    expect_equal(!!outs[[i]], !!test_formulas_glmm[[i]])
  })
})

test_that("extract_variables_from_formula throws error if incorrect",{
  expect_error(is_valid_censored_formula(~Surv(x,i),throw_error = TRUE))
  expect_error(is_valid_censored_formula(Surv(x,i)~y,throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~x,throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~x+z,throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~1,throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~Surv(,i),throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~Surv(x,),throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~Surv(,),throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~Surv( , ),throw_error = TRUE))
  expect_error(is_valid_censored_formula(y~Surv(x,i)+ Surv(a,j),throw_error = TRUE))
})

test_that("extract_variables_from_formula is FALSE if incorrect",{
  expect_false(is_valid_censored_formula(~Surv(x,i),throw_error = FALSE))
  expect_false(is_valid_censored_formula(Surv(x,i)~y,throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~x,throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~x+z,throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~1,throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~Surv(,i),throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~Surv(x,),throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~Surv(,),throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~Surv( , ),throw_error = FALSE))
  expect_false(is_valid_censored_formula(y~Surv(x,i)+ Surv(a,j),throw_error = FALSE))
})
