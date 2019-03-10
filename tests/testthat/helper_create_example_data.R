# Helper functions to create example data for unit tests

# Function to create random data (one sample)
d_random <- function(n = 20000, mean = 0, sd = 1, ncol = 20, cofactor = 5) {
  d <- sinh(matrix(rnorm(n, mean, sd), ncol = ncol)) * cofactor
  colnames(d) <- paste0("marker", sprintf("%02d", 1:ncol))
  d
}

# Create random data (without differential signal)
set.seed(123)
d_input <- list(
  sample1 = d_random(), 
  sample2 = d_random(), 
  sample3 = d_random(), 
  sample4 = d_random()
)

experiment_info <- data.frame(
  sample_id = factor(paste0("sample", 1:4)), 
  group_id = factor(c("group1", "group1", "group2", "group2")), 
  stringsAsFactors = FALSE
)

marker_info <- data.frame(
  channel_name = paste0("channel", sprintf("%03d", 1:20)), 
  marker_name = paste0("marker", sprintf("%02d", 1:20)), 
  marker_class = factor(c(rep("type", 10), rep("state", 10)), 
                        levels = c("type", "state", "none")), 
  stringsAsFactors = FALSE
)
