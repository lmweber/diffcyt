################################################################################
### linear_sum
n <- 10
covariates_list <- list(z1=rnorm(n),z2=rnorm(n),z3=rnorm(n))
covariates_list_err <- list(z1=rnorm(n+1),z2=rnorm(n),z3=rnorm(n-1))
random_covariates_list <- list(r1=rnorm(n),r2=rnorm(n),r3=rnorm(n))
random_covariates_list_err <- list(r1=rnorm(n+1),r2=rnorm(n),r3=rnorm(n-1))
betas <- c(1,1,1)
linear_sum(covariates_list = covariates_list,betas = betas,random_covariates_list = random_covariates_list)

test_that("linear_sum correct",{
  expect_equal(length(linear_sum(covariates_list,betas)), n)
  expect_equal(length(linear_sum(covariates_list,betas,random_covariates_list)), n)
  expect_equal(linear_sum(covariates_list,c(0,0,0)), matrix(0,ncol = 1,nrow = n))
  expect_equal(linear_sum(covariates_list,c(1,0,0)), matrix(covariates_list[[1]],ncol = 1,nrow = n))
  expect_equal(linear_sum(covariates_list,c(0,1,0)), matrix(covariates_list[[2]],ncol = 1,nrow = n))
  expect_equal(linear_sum(covariates_list,c(0,0,1)), matrix(covariates_list[[3]],ncol = 1,nrow = n))
  expect_equal(linear_sum(covariates_list,c(1,1,0)), matrix(covariates_list[[1]]+covariates_list[[2]],ncol = 1,nrow = n))
  expect_equal(linear_sum(covariates_list,c(0,0,0),random_covariates_list[1]), matrix(random_covariates_list[[1]],ncol = 1,nrow = n))
  expect_error(linear_sum(covariates_list,1))
  expect_error(linear_sum(covariates_list,rep(1,4)))
  expect_error(linear_sum(covariates_list_err,betas))
  expect_error(linear_sum(covariates_list,betas,random_covariates_list_err))
})

################################################################################
### simulate_singlecluster
lm_formula <- formula(Y ~ Surv(X,I) + Z)
glm_formula <- formula(Y ~ Surv(X,I) + Z)
glmer_formula <- formula(Y ~ Surv(X,I) + Z + (1|R))

test_that("simulate_singlecluster correct class",{
   expect_equal(class(simulate_singlecluster(100, lm_formula, type = "lm")),c("tbl_df","tbl","data.frame"))
   expect_equal(class(simulate_singlecluster(100, glm_formula, type = "glm")),c("tbl_df","tbl","data.frame"))
   expect_equal(class(simulate_singlecluster(100, glmer_formula, type = "glmer")),c("tbl_df","tbl","data.frame"))
   })

test_that("simulate_singlecluster correct size",{
   expect_error(simulate_singlecluster(-1, lm_formula, type = "lm"))
   expect_error(simulate_singlecluster(0, lm_formula, type = "lm"))
   expect_equal(dim(simulate_singlecluster(2, lm_formula, type = "lm")),c(2,8))
   expect_equal(dim(simulate_singlecluster(10, lm_formula, type = "lm")),c(10,8))
   expect_equal(dim(simulate_singlecluster(100, lm_formula, type = "lm")),c(100,8))
 })

test_that("simulate_singlecluster correct n_levels_fixeff for lm",{
   expect_error( length(unique(simulate_singlecluster(60, lm_formula, type = "lm",n_levels_fixeff = 0)$Z)))
   expect_equal( length(unique(simulate_singlecluster(60, lm_formula, type = "lm",n_levels_fixeff = 1)$Z)),1)
   expect_equal( length(unique(simulate_singlecluster(60, lm_formula, type = "lm",n_levels_fixeff = 2)$Z)),2)
   expect_equal( length(unique(simulate_singlecluster(60, lm_formula, type = "lm",n_levels_fixeff = 3)$Z)),3)
   expect_equal( length(unique(simulate_singlecluster(60, lm_formula, type = "lm",n_levels_fixeff = 60)$Z)),60)
 })
test_that("simulate_singlecluster correct n_levels_fixeff for glm",{
   expect_error( length(unique(simulate_singlecluster(60, glm_formula, type = "glm",n_levels_fixeff = 0)$Z)))
   expect_equal( length(unique(simulate_singlecluster(60, glm_formula, type = "glm",n_levels_fixeff = 1)$Z)),1)
   expect_equal( length(unique(simulate_singlecluster(60, glm_formula, type = "glm",n_levels_fixeff = 2)$Z)),2)
   expect_equal( length(unique(simulate_singlecluster(60, glm_formula, type = "glm",n_levels_fixeff = 3)$Z)),3)
   expect_equal( length(unique(simulate_singlecluster(60, glm_formula, type = "glm",n_levels_fixeff = 60)$Z)),60)
 })
test_that("simulate_singlecluster correct n_levels_fixeff for glmer",{
   expect_error( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_fixeff = 0)$Z)))
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_fixeff = 1)$Z)),1)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_fixeff = 2)$Z)),2)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_fixeff = 3)$Z)),3)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_fixeff = 60)$Z)),60)
 })

test_that("simulate_singlecluster correct n_levels_raneff for glmer",{
   expect_error( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 0)$R)))
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 1)$R)),1)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 2)$R)),2)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 3)$R)),3)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 60)$R)),60)
 })

test_that("simulate_singlecluster correct n_levels_fixeff and n_levels_raneff for glmer",{
   expect_error( simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 0,n_levels_fixeff = 0))
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 1,n_levels_fixeff = 1)$R)),1)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 1,n_levels_fixeff = 1)$Z)),1)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 2,n_levels_fixeff = 2)$R)),2)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 2,n_levels_fixeff = 2)$Z)),2)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 3,n_levels_fixeff = 3)$R)),3)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 3,n_levels_fixeff = 3)$Z)),3)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 60,n_levels_fixeff = 60)$R)),60)
   expect_equal( length(unique(simulate_singlecluster(60, glmer_formula, type = "glmer",n_levels_raneff = 60,n_levels_fixeff = 60)$Z)),60)
 })

test_that("simulate_singlecluster throw error for negative shapes or scales",{
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params = list(X = list(shape = -1, scale = 2),
                                                                                  C = list(shape = 0.5, scale = 10))))
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params = list(X = list(shape = 1, scale = -2),
                                                                                  C = list(shape = 0.5, scale = 10))))
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params = list(X = list(shape = 1, scale = 2),
                                                                                  C = list(shape = -0.5, scale = 10))))
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params = list(X = list(shape = 1, scale = 2),
                                                                                  C = list(shape = 0.5, scale = -10))))
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params_covariate_dependent_censoring = list(shape = -1, scale = 10)))
   expect_error(simulate_singlecluster(100, lm_formula, type = "lm", weibull_params_covariate_dependent_censoring = list(shape = 1, scale = -10)))
 })


test_that("simulate_singlecluster correct output censoring indicator",{
  expect_true(all(simulate_singlecluster(100, lm_formula, type = "lm",censoring_dependent_on_covariate = TRUE)$I %in% c(0,1)))
 })



simulate_singlecluster_params_beta <- function(b,glmer_formula){
   simulate_singlecluster(
       n = 10,
       formula = glmer_formula,
       n_levels_fixeff = c(2),
       n_levels_raneff = NULL,
       type = "glmer",
       b = b,
       weibull_params = list(X = list(shape = 0.75, scale = 0.25),
                             C = list(shape = 0.5, scale = 0.7)),
       censoring_dependent_on_covariate = FALSE,
       weibull_params_covariate_dependent_censoring = list(shape = 0.5, scale = 0.9),
       error_variance = 0.0,
       variance_raneff = 0.5,
       transform_fn = "boxcox_positive")
}

b_l2 <- c(-4,-1.5)
b_l3 <- c(-4,-1.5,0.5)
b_l4 <- c(-4,-1.5,0.5,0.5)
glmer_formula <- formula(Y ~ Surv(X,I) + (1|R))
test_that("error if no covariate",{
   expect_error(simulate_singlecluster_params_beta(list(b_l2),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l2,b_l2),2,glmer_formula))
   expect_error(simulate_singlecluster_params_beta(list(b_l3),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l3,b_l3),2,glmer_formula))
   expect_error(simulate_singlecluster_params_beta(list(b_l4),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l4,b_l4),2,glmer_formula))
})

glmer_formula <- formula(Y ~ Surv(X,I) + Z + (1|R))
test_that("list of betas matches with formula for 1 covariate",{
   expect_error(simulate_singlecluster_params_beta(list(b_l2),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l2,b_l2),2,glmer_formula))
   expect_true(is(simulate_singlecluster_params_beta(list(b_l3),glmer_formula), "data.frame"))
   # expect_true(is(simulate_singlecluster_params_beta(list(b_l3,b_l3),2,glmer_formula), "SummarizedExperiment"))
   expect_error(simulate_singlecluster_params_beta(list(b_l4),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l4,b_l4),2,glmer_formula))
})

glmer_formula <- formula(Y ~ Surv(X,I) + Z + ZZ + (1|R))
test_that("list of betas matches with formula for 2 covariates",{
   expect_error(simulate_singlecluster_params_beta(list(b_l2),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l2,b_l2),2,glmer_formula))
   expect_error(simulate_singlecluster_params_beta(list(b_l3),glmer_formula))
   # expect_error(simulate_singlecluster_params_beta(list(b_l3,b_l3),2,glmer_formula))
   expect_true(is(simulate_singlecluster_params_beta(list(b_l4),glmer_formula), "data.frame"))
   # expect_true(is(simulate_singlecluster_params_beta(list(b_l4,b_l4),2,glmer_formula), "SummarizedExperiment"))
})




################################################################################
### calculate_wanted_cluster_proportion
 test_that("calculate_wanted_cluster_proportion correct",{
    expect_equal(length(calculate_wanted_cluster_proportion(c(0.1))),1)
    expect_equal(length(calculate_wanted_cluster_proportion(c(0.1,0.1))),2)
    expect_equal(length(calculate_wanted_cluster_proportion(c(0.1,0.1,0.1))),3)
    expect_error(calculate_wanted_cluster_proportion(c(0.2,0.8,0.01)))
    expect_error(calculate_wanted_cluster_proportion(c(0.4,0.4)))
    expect_equal(calculate_wanted_cluster_proportion(c(0.1,0.1)),c(0.5,0.5))
    expect_equal(calculate_wanted_cluster_proportion(c(0.1)),1)
 })

################################################################################
### do_transformation
 test_that("do_transformation correct",{
    expect_equal(do_transformation(1:10,"identity"),1:10)
    expect_equal(do_transformation(1:10,function(x){x}),1:10)
    expect_equal(do_transformation(1:10,function(x){x+1}),2:11)
    expect_equal(do_transformation(1:10,"log"),log(1:10))
    expect_equal(do_transformation(1:10,log),log(1:10))
    expect_equal(do_transformation(1:10,"log_positive"),log(1:10)+abs(min(log(1:10))))
    expect_equal(do_transformation(seq(0,1,by=0.01),"log_positive"),log(seq(0,1,by=0.01))+abs(min(log(seq(0,1,by=0.01)))))
    expect_equal(do_transformation(1:10,sqrt),sqrt(1:10))

 })

