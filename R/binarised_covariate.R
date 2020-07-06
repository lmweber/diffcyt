
binarised_covariate <- function(covariate, method="median"){
  if (method=="median"){
    ind <- covariate > median(covariate)
  } else {
    quantiles <- quantile(covariate,probs=seq(0,1,length.out = length(covariate)))
    out <- lapply(quantiles,function(q){
      ind <- covariate>=q
      data.frame(m1=mean(covariate[ind]),m2=mean(covariate[!ind]))
    })
    out_df <- purrr::reduce(out,rbind)
    
    slope <- (out_df[-1,2]-out_df[-dim(out_df)[1],2])/((out_df[-1,1]-out_df[-dim(out_df)[1],1]))
    subsetting <- ceiling(dim(out_df)[1]*0.25):ceiling(dim(out_df)[1]*0.75)
    max_slope <- which.max(slope[subsetting])+min(subsetting)-1
    slope[max_slope]
    ind <- covariate>quantiles[max_slope]
  }
  return(factor(ind,levels = c(TRUE,FALSE),labels = c("high","low")))
}

