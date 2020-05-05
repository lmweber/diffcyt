syn_data <- tibble::tibble(X=1:6,censored=rep(c(0,1),3), Z = rep(c(0,1),each=3))
syn_data2 <- tibble::tibble(X=1:8,censored=rep(c(0,1),4), Z = rep(c(0,1),each=4))
syn_data3<- tibble::tibble(X=1:6,censored=rep(0,6), Z = rep(c(0,1),each=3))
syn_data4<- tibble::tibble(X=1:6,censored=rep(1,6), Z = rep(c(0,1),each=3))

test_that("risk_set_imputation correct output format",{
  expect_equal(typeof(risk_set_imputation(syn_data,"X","censored")),"double")
  expect_equal(length(risk_set_imputation(syn_data,"X","censored")),sum(syn_data$censored==0))
  expect_equal(length(risk_set_imputation(syn_data2,"X","censored")),sum(syn_data2$censored==0))
  expect_equal(length(risk_set_imputation(syn_data3,"X","censored")),sum(syn_data3$censored==0))
  expect_equal(length(risk_set_imputation(syn_data4,"X","censored")),sum(syn_data4$censored==0))
  expect_equal(typeof(suppressWarnings(risk_set_imputation(syn_data,"X","censored","Z"))),"double")
  expect_equal(length(suppressWarnings(risk_set_imputation(syn_data,"X","censored","Z"))),sum(syn_data$censored==0))
  expect_equal(length(suppressWarnings(risk_set_imputation(syn_data2,"X","censored","Z"))),sum(syn_data2$censored==0))
  expect_error(suppressWarnings(risk_set_imputation(syn_data3,"X","censored","Z")))
  expect_equal(length(suppressWarnings(risk_set_imputation(syn_data4,"X","censored","Z"))),sum(syn_data4$censored==0))
})

test_that("risk_set_imputation correct output values",{
  set.seed(123)
  expect_equal(risk_set_imputation(syn_data,"X","censored"),c(4,5,3))
  set.seed(123)
  expect_equal(suppressWarnings(risk_set_imputation(syn_data,"X","censored","Z")),c(4,4,3))
})


test_that("risk_set_cov_adjusted correct output format",{
  expect_equal(class(suppressWarnings(risk_set_cov_adjusted(syn_data,"X","censored", "Z"))),"list")
  expect_equal(class(suppressWarnings(risk_set_cov_adjusted(syn_data2,"X","censored", "Z"))),"list")
  expect_equal(class(suppressWarnings(risk_set_cov_adjusted(syn_data4,"X","censored", "Z"))),"list")
})
test_that("risk_set_cov_adjusted correct output values",{
  expect_equal(suppressWarnings(risk_set_cov_adjusted(syn_data,"X","censored", "Z")),list(c(2,3,4,5),c(4,5),6))
})
