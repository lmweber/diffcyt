##########################################################################################
# Script to save small data set for examples (small subset of 'BCR-XL' data set from 
# Bodenmiller et al. 2012)
##########################################################################################


library(flowCore)
library(magrittr)


# filenames
files <- list.files("../../benchmark_data/BCR_XL/raw_data/experiment_15713_files", full.names = TRUE)

# groups: BCR-XL vs. reference
files_BCRXL <- files[grep("patient[1-8]_BCR-XL\\.fcs$", files)]
files_ref <- files[grep("patient[1-8]_Reference\\.fcs$", files)]

# load data
files_load <- c(files_BCRXL, files_ref)
data_BCRXL <- lapply(files_load, read.FCS, transformation = FALSE, truncate_max_range = FALSE)

# sample IDs
basename(files_load) %>% 
  gsub("^PBMC8_30min_|\\.fcs$", "", .) %>% 
  gsub("BCR-XL", "BCRXL", .) %>% 
  gsub("Reference", "ref", .) -> 
  sample_IDs

sample_IDs

names(data_BCRXL) <- sample_IDs


# subsample to create small example data set (6 samples, 100 cells per sample)
smp_keep <- sample_IDs[c(1:3, 9:11)]
smp_keep

n <- 100

for (i in 1:length(smp_keep)) {
  data_i <- flowFrame(exprs(data_BCRXL[[i]])[1:n, ])
  fn <- paste0("../inst/extdata/", smp_keep[i], ".fcs")
  write.FCS(data_i, filename = fn)
}


