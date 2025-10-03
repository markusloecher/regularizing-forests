library(dplyr)
library(tidyr)
library(ggplot2)


# # Load shrinkage + bagging results
# load("data/bagging_mse_concrete.RData")  # final_df for bag+shrink
# 
# # Load MSE grid search results
# load("data/MSE_out_concrete.RData")     # MSE_out, has d_RF
# 
# # Extract constants
# p <- MSE_out[[1]]$p
# n <- MSE_out[[1]]$n
# var_y <- var(MSE_out[[1]]$dat[, p + 1])
# nrep <- length(MSE_out)
# 
# # Initialize
# d_RF <- d_mtry <- d_depth <- d_bagging <- numeric(nrep)
# 
# rte_dft_rf <- numeric(nrep)
# 
# for (i in seq_len(nrep)) {
#   mse_mat <- MSE_out[[i]]$cv.sigma
#   
#   row_vals <- as.numeric(gsub("mtry = ", "", rownames(mse_mat)))
#   col_vals <- as.numeric(gsub("maxnodes = ", "", colnames(mse_mat)))
#   
#   mtry_dft <- which.min(abs(row_vals - 0.3 * p))
#   depth_dft <- which.min(abs(col_vals - 1.0 * n))  # full trees
#   
#   rte_mat <- mse_mat / var_y
#   
#   rte_def <- rte_mat[mtry_dft, depth_dft]
#   rte_opt <- min(rte_mat, na.rm = TRUE)
#   
#   d_RF[i]     <- rte_def - rte_opt #d_RF found
#   rte_dft_rf[i] <- rte_def
# }
# 
# 
# 
# #Use -tunability + d_RF to find tunability of bagging+shrinkage
# #diff_vec <- 100* (-final_df$tunability + d_RF)/final_df$rte_def
# diff_vec <- 100 * (-final_df$tunability + d_RF)/rte_dft_rf
# 
# final_df$rel_tunability_bagging <- diff_vec
# 
# save(final_df, file = "data/final_df_tuned_concrete.RData")



# List all bagging result files
bagging_files <- list.files("data/", pattern = "^bagging_mse_.*\\.RData$", full.names = TRUE)

for (bag_file in bagging_files) {
  dataname <- sub("bagging_mse_(.*)\\.RData", "\\1", basename(bag_file))
  
  # Corresponding MSE_out file
  mse_file <- file.path("data", paste0("MSE_out_", dataname, ".RData"))
  if (!file.exists(mse_file)) next  # Skip if MSE file not found
  
  # Load both RData files
  load(bag_file)    # loads final_df
  load(mse_file)    # loads MSE_out
  
  # Extract constants
  p <- MSE_out[[1]]$p
  n <- MSE_out[[1]]$n
  var_y <- var(MSE_out[[1]]$dat[, p + 1])
  nrep <- length(MSE_out)
  
  d_RF <- rte_dft_rf <- numeric(nrep)
  
  for (i in seq_len(nrep)) {
    mse_mat <- MSE_out[[i]]$cv.sigma
    row_vals <- as.numeric(gsub("mtry = ", "", rownames(mse_mat)))
    col_vals <- as.numeric(gsub("maxnodes = ", "", colnames(mse_mat)))
    
    mtry_dft <- which.min(abs(row_vals - 0.3 * p))
    depth_dft <- which.min(abs(col_vals - 1.0 * n))  # full trees
    
    rte_mat <- mse_mat / var_y
    rte_def <- rte_mat[mtry_dft, depth_dft]
    rte_opt <- min(rte_mat, na.rm = TRUE)
    
    d_RF[i] <- rte_def - rte_opt
    rte_dft_rf[i] <- rte_def
  }
  
  # Compute and store relative tunability gain
  diff_vec <- 100 * (-final_df$tunability + d_RF) / rte_dft_rf
  final_df$rel_tunability_bagging <- diff_vec
  
  # Save updated final_df
  save(final_df, file = file.path("data", paste0("final_df_tuned_", dataname, ".RData")))
}





# #SANITY CHECKS
# load("data/bagging_mse_concrete.RData")
# 
# nrep = 50
# 
# df <- data.frame(
#   rep = 1:nrep,
#   BagHS = final_df$rte_def,
#   RF = final_df$rte_opt
# )
# 
# df_long <- pivot_longer(df, cols = c("BagHS", "RF"), names_to = "Method", values_to = "RTE")
# 
# ggplot(df_long, aes(x = Method, y = RTE, color = Method)) +
#   geom_boxplot() +
#   ylab("Relative Test Error") +
#   ggtitle("Comparison of RTE: Bagging+Shrinkage vs Tuned RF")
# 
# #def here is bagging + hs and opt is the optimally tuned rf
# diff_vec <- final_df$rte_def - final_df$rte_opt
# 
# ggplot(data.frame(diff = diff_vec), aes(x = diff)) +
#   geom_histogram(bins = 30, fill = "steelblue", color = "white") +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
#   xlab("RTE(BagHS) - RTE(RF)") +
#   ggtitle("Distribution of RTE Differences")
# 
# #Stability across lambda
# #unimodal = better, multimodel = probably noisy data
# plot(final_df$best_lambda, final_df$rte_def)






##coefficient of variation

# Directory containing your files
file_list <- list.files(path = "data/", pattern = "^final_df_tuned_.*\\.RData$", full.names = TRUE)

# Initialize a list to collect results
cv_results <- list()

for (file in file_list) {
  load(file)  # loads final_df

  # Extract dataname
  dataname <- sub("final_df_tuned_(.*)\\.RData", "\\1", basename(file))

  # Compute coefficient of variation
  cv <- sd(final_df$mtry_opt) / mean(final_df$mtry_opt)

  # Store result
  cv_results[[dataname]] <- cv
}

# Convert to data frame
cv_df <- data.frame(
  dataname = names(cv_results),
  coefficient_of_variation = as.numeric(cv_results),
  row.names = NULL
)

#print(cv_df)
