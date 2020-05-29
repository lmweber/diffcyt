syn_data <- tibble::tibble(X=as.double(1:6),censored=rep(c(0,1),3), Z = rep(c(0,1),each=3))
dim_cens <- c(3,1)
dim_cens_2reps <- c(3,2)

test_that("kaplan_meier_imputation with survival tail method 'lao' works when last value is observed",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z"))),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",mi_reps = 2)), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'exp' works when last value is observed",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp"))),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp",mi_reps = 2)), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'wei' works when last value is observed",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei"))),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei",mi_reps = 2)), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'os' works when last value is observed",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os"))),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os",mi_reps = 2)), dim_cens_2reps)
})


syn_data <- tibble::tibble(X=as.double(1:6),censored=rep(c(1,0),3), Z = rep(c(0,1),each=3))
dim_cens <- c(3,1)
dim_cens_2reps <- c(3,2)

test_that("kaplan_meier_imputation with survival tail method 'lao' works when last value is censored",{
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z"))), dim_cens)
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",mi_reps = 2))), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'exp' works when last value is censored",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp"))),"double")
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp"))), dim_cens)
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "exp",mi_reps = 2))), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'wei' works when last value is censored",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "wei",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei"))),"double")
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei"))), dim_cens)
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "wei",mi_reps = 2))), dim_cens_2reps)
})

test_that("kaplan_meier_imputation with survival tail method 'os' works when last value is censored",{
  expect_equal(typeof(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os")),"double")
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os")), dim_cens)
  expect_equal(dim(kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "os",mi_reps = 2)), dim_cens_2reps)
  
  expect_equal(typeof(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os"))),"double")
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os"))), dim_cens)
  expect_equal(dim(suppressWarnings(kaplan_meier_imputation(syn_data,"X","censored","Z",tail_approx_method = "os",mi_reps = 2))), dim_cens_2reps)
})








syn_data <- tibble::tibble(X=as.double(1:8),censored=rep(c(0,1),4), Z = rep(rep(c(0,1),each=2),2))
test_that("mean_residual_life_imputation outcome values greater or equal",{
  mrl_out <- mean_residual_life_imputation(syn_data,"X","censored")
  ref <- syn_data[syn_data$censored==0, ]$X
  purrr::map2(mrl_out, ref, ~ expect_gte(.x, .y))
  mrl_out <- mean_residual_life_imputation(syn_data,"X","censored","Z")
  purrr::map2(mrl_out, ref, ~ expect_gte(.x, .y))
  })


test_that("estimate_cens_vars 'rs' works",{
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "rs")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "rs")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- risk_set_imputation(syn_data,"X","censored")
  expect_equal(ecov,ecov)
  
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z",imputation_method = "rs")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- risk_set_imputation(syn_data,"X","censored", "Z")
  expect_equal(ecov,ecov)
})
test_that("estimate_cens_vars 'km' works",{
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "km")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "km")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored")
  expect_equal(ecov,ecov)
  
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z", imputation_method = "km")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", "Z")
  expect_equal(ecov,ecov)
})
test_that("estimate_cens_vars 'km_exp' works",{  
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_exp")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_exp")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored",tail_approx_method = "exp")
  expect_equal(ecov,ecov)
  
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z", imputation_method = "km_exp")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", "Z", tail_approx_method = "exp")
  expect_equal(ecov,ecov)
})
test_that("estimate_cens_vars 'km_wei' works",{  
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_wei")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_wei")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", tail_approx_method = "wei")
  expect_equal(ecov,ecov)
  
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z", imputation_method = "km_wei")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", "Z", tail_approx_method = "wei")
  expect_equal(ecov,ecov)
})
test_that("estimate_cens_vars 'km_os' works",{
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_os")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "km_os")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", tail_approx_method = "os")
  expect_equal(ecov,ecov)
  
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z", imputation_method = "km_os")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- kaplan_meier_imputation(syn_data,"X","censored", "Z", tail_approx_method = "os")
  expect_equal(ecov,ecov)
})
test_that("estimate_cens_vars 'mrl' works",{
  expect_equal(class(estimate_cens_vars(syn_data,"X","censored",imputation_method = "mrl")),c("tbl_df","tbl","data.frame"))
  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored",imputation_method = "mrl")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- mean_residual_life_imputation(syn_data,"X","censored")
  expect_equal(ecov,ecov)

  set.seed(123)
  ecov <- estimate_cens_vars(syn_data,"X","censored", "Z" ,imputation_method = "mrl")
  ecov <- ecov[ecov$"censored"==0, ]$"X_est"
  set.seed(123)
  comp <- mean_residual_life_imputation(syn_data,"X","censored", "Z")
  expect_equal(ecov,ecov)
})

testtib_to_test  <- tibble::tibble(x=factor(rep(c(1,2),3)),y=factor(rep(c(1,2,3),2)),z=1:6,u=rep(c(1,2),3))
testtib_expected <- tibble::tibble(x=rep(c(0,1),3),y=factor(rep(c(1,2,3),2)),z=1:6,u=rep(c(1,2),3))

test_that("covariates_from_factor_to_numeric works",{
  expect_equal(covariates_from_factor_to_numeric(testtib_to_test,c("x")),
               testtib_expected)
  expect_equal(covariates_from_factor_to_numeric(testtib_to_test,c("")),
               testtib_to_test)
  expect_equal(covariates_from_factor_to_numeric(testtib_to_test,c("wrong_name")),
               testtib_to_test)
  expect_equal(covariates_from_factor_to_numeric(testtib_to_test,c("x","z","u")),
               testtib_expected)
  expect_error(covariates_from_factor_to_numeric(testtib_to_test,c("x","y","z","u")))
})
syn_data_zfac <- dplyr::mutate(syn_data,Zf=factor(Z),Zfm=factor(Z*10),Zfe=factor(X))
test_that("covariates_from_factor_to_numeric works",{
  expect_equal(covariates_from_factor_to_numeric(syn_data_zfac,"Zf")$Zf, syn_data_zfac$Z)
  expect_equal(covariates_from_factor_to_numeric(syn_data_zfac,"Zfm")$Zfm, syn_data_zfac$Z)
  expect_equal(covariates_from_factor_to_numeric(syn_data_zfac,c("Zf","Zfm"))$Zf, syn_data_zfac$Z)
  expect_equal(covariates_from_factor_to_numeric(syn_data_zfac,c("Zf","Zfm"))$Zfm, syn_data_zfac$Z)
  expect_error(covariates_from_factor_to_numeric(syn_data_zfac,"Zfe"))
})

