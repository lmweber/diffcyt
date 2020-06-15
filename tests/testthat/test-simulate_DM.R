nr_samples <-  10
nr_cluster <- 100
nr_diff  <-  4
sizes  <-  runif(nr_samples,1e5,1e6)
theta  <-  0.001
  
# simulate non DA data
alphas_theo <- runif(nr_cluster,0,1)*(1-theta)/theta
probs <- alphas_theo/sum(alphas_theo)
counts <- t(dirmult::simPop(J=nr_samples,n=sizes,pi=probs,theta=theta)$data)



test_that("dimensionality is kept in simulate_multicluster",{
  expect_equal(dim(simulate_multicluster(counts,alphas = alphas_theo,theta = theta)$counts),c(nr_cluster,nr_samples))
  expect_equal(dim(simulate_multicluster(counts,nr_samples = 20,alphas = alphas_theo,theta = theta)$counts),c(100,20))
  expect_equal(dim(simulate_multicluster(counts,nr_samples = 8,alphas = alphas_theo,theta = theta)$counts),c(100,8))
})


test_that("output same as input in simulate_multicluster if from reference",{
  out <- simulate_multicluster(counts=counts,
                             nr_diff = nr_diff,
                             nr_samples = NULL,
                             alphas = alphas_theo,
                             theta = theta,
                             covariate = NULL,
                             slope = NULL,
                             random_cluster_choice = FALSE,
                             enforce_sum_alpha = TRUE)
  expect_equal(dim(out$counts),c(nr_cluster,nr_samples))
  expect_equal(length(out$diff_clus),nr_diff)
  expect_equal(dim(out$betas),c(2,nr_diff))
  expect_equal(length(out$z),nr_samples)
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
                                    random_cluster_choice = FALSE,
                                    enforce_sum_alpha = TRUE)
  expect_equal(dim(out$counts),c(nr_cluster,nr_samples))
  expect_equal(length(out$diff_clus),nr_diff)
  expect_equal(dim(out$betas),c(2,nr_diff))
  expect_equal(length(out$z),nr_samples)
  expect_equal(dim(out$alphas),c(nr_cluster,nr_samples))
  expect_equal(length(out$theta),1)
  expect_equal(dim(out$var_counts),c(nr_cluster,nr_samples))
})

