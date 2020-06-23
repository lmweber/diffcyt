nr_samples <-  10
nr_cluster <- 10
nr_diff  <-  4
sizes  <-  runif(nr_samples,1e5,1e6)
theta  <-  0.001
  
# simulate non DA data
alphas_theo <- runif(nr_cluster,0,1)*(1-theta)/theta
probs <- alphas_theo/sum(alphas_theo)
counts <- t(dirmult::simPop(J=nr_samples,n=sizes,pi=probs,theta=theta)$data)


test_that("dimensionality is kept in simulate_multicluster",{
  expect_equal(dim(simulate_multicluster(alphas = alphas_theo,sizes = sizes)$counts),c(nr_cluster,nr_samples))
  expect_equal(dim(simulate_multicluster(alphas = alphas_theo,sizes = sample(sizes,replace = TRUE,20))$counts),c(nr_cluster,20))
  expect_equal(dim(simulate_multicluster(alphas = alphas_theo,sizes =  sample(sizes,replace = TRUE,8))$counts),c(nr_cluster,8))
})

test_that("output same as input in simulate_multicluster if from reference",{
  out <- simulate_multicluster(counts=counts,
                             nr_diff = nr_diff,
                             nr_samples = NULL,
                             # alphas = alphas_theo,
                             # sizes = sizes,
                             covariate = NULL,
                             slope = NULL,
                             diff_cluster = FALSE,
                             enforce_sum_alpha = TRUE)
  expect_equal(dim(out$counts),c(nr_cluster,nr_samples))
  expect_equal(sum(!is.na(out$row_data$paired)),nr_diff)
  expect_equal(sum(out$row_data$b1 != 0),4)
  expect_equal(length(out$col_data$covariate),nr_samples)
  expect_equal(dim(out$alphas),c(nr_cluster,nr_samples))
  expect_equal(length(out$theta),1)
  expect_equal(dim(out$var_counts),c(nr_cluster,nr_samples))
})

test_that("output same as input in simulate_multicluster if from parameters",{
  out <- simulate_multicluster(counts=NULL,
                                    nr_diff = nr_diff,
                                    nr_samples = NULL,
                                    alphas = alphas_theo,
                                    # theta = theta,
                                    sizes = runif(10,1e4,1e5),
                                    covariate = NULL,
                                    slope = NULL,
                                    diff_cluster = FALSE,
                                    enforce_sum_alpha = TRUE)
  expect_equal(dim(out$counts),c(nr_cluster,nr_samples))
  expect_equal(sum(!is.na(out$row_data$paired)),nr_diff)
  expect_equal(sum(out$row_data$b1 != 0),4)
  expect_equal(length(out$col_data$covariate),nr_samples)
  expect_equal(dim(out$alphas),c(nr_cluster,nr_samples))
  expect_equal(length(out$theta),1)
  expect_equal(dim(out$var_counts),c(nr_cluster,nr_samples))
})
test_that("list input in 'slope' and 'group_slope' parameter work in simulate_multicluster",{
  suppressWarnings({
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5),slope = 0))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5),slope = 0.1))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5),slope = list(-0.2,0.2)))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5),group_slope = list(-0.2,0.2)))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5), group_slope = 0))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5), group_slope = 0.2))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5), 
                                       slope= list(0.8,0.2), 
                                       group_slope = list(0.2,0.2)))
    expect_error(simulate_multicluster(nr_diff = nr_diff,alphas = alphas_theo,sizes = runif(10,1e4,1e5), 
                                       slope= list(0.7,0.2), 
                                       group_slope = list(0.2,0.9)))
  })
})
