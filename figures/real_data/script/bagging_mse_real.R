library(randomForest)
library(ranger)
library(future.apply)
plan(multisession, workers = 25)

source("functions/functions_bagging_mse.R")

simulate_rep <- function(rep, MSE_out) {
  dat <- MSE_out[[1]]$dat
  p <- MSE_out[[1]]$p
  n <- MSE_out[[1]]$n
  params <- generate_params(n)
  #print(params)
  mse_mat <- MSE_out[[rep]]$cv.sigma
  
  opt_idx <- which(mse_mat == min(mse_mat), arr.ind = TRUE)
  mtry_opt <- as.numeric(gsub("mtry = ", "", rownames(mse_mat)[opt_idx[,"row"]]))
  maxnodes_opt <- as.numeric(gsub("maxnodes = ", "", colnames(mse_mat)[opt_idx[,"col"]])) 
  
  test_idx <- sample(1:n, n / 2)
  n_train <- n - length(test_idx)
  maxnodes_used <- min(maxnodes_opt, n_train)
  
  #maybe not use this to begin with
  rf_opt <- randomForest(
    x = dat[-test_idx, 1:p],
    y = dat[-test_idx, p + 1],
    xtest = dat[test_idx, 1:p],
    ytest = dat[test_idx, p + 1],
    mtry = mtry_opt,
    nodesize = 1,
    maxnodes = maxnodes_used,
    ntree = 500
  )
  
  mse_opt <- mean((rf_opt$test$predicted - dat[test_idx, p + 1])^2)
  var_y <- var(dat[, p + 1])
  rte_opt <- mse_opt / var_y
  
  best_lambda <- NA
  mse_def <- Inf

  for (lambda in params) {
    rf_bag <- ranger(
      x = dat[-test_idx, 1:p],
      y = dat[-test_idx, p + 1],
      num.trees = 500,
      mtry = p,
      min.node.size = 1,
      max.depth = NULL,
      seed = rep,
      node.stats = TRUE,
      keep.inbag = TRUE
    )

    hshrink(rf_bag, lambda = lambda)
    preds <- predict(rf_bag, data = dat[test_idx, 1:p])$predictions
    mse <- mean((preds - dat[test_idx, p + 1])^2)

    if (mse < mse_def) {
      mse_def <- mse #bagging baseline mse from best performing lambda
      best_lambda <- lambda
    }
  }

  rte_def <- mse_def / var_y
  d_lambda <- rte_def - rte_opt
  
  data.frame(
    rep = rep,
    maxnodes_opt = maxnodes_opt,
    mtry_opt = mtry_opt,
    mse_opt = mse_opt,
    rte_opt = rte_opt,
    best_lambda = best_lambda,
    mse_def = mse_def,
    rte_def = rte_def,
    tunability = d_lambda
  )
}


files <- list.files("data", pattern = "^MSE_out_.*\\.RData$", full.names = TRUE)

#test
#files <- c("data/MSE_out_cpu.RData", "data/MSE_out_AquaticTox.RData")

for (file in files) {
  load(file)  # loads MSE_out
  dataname <- sub("^MSE_out_(.*)\\.RData$", "\\1", basename(file))
  
  reps <- 1:50
  results <- future_lapply(reps, simulate_rep, MSE_out = MSE_out, future.seed = TRUE)
  final_df <- do.call(rbind, results)
  
  save(final_df, file = file.path("data", paste0("bagging_mse_", dataname, ".RData")))
  cat("Saved results for:", dataname, "\n")
}



# load("data/MSE_out_.RData")
# length(colnames(MSE_out[[1]]$cv.sigma))


###Visual tool

library(ggplot2)
library(dplyr)
library(tools)
library(ggrepel)

#set the directory with the RData files
data_dir <- "data"  # change this to the actual path
files <- list.files(path = data_dir, pattern = "^final_df_tuned_.*\\.RData$", full.names = TRUE)

#initialize result collector
result_list <- list()

#loop over files
for (file in files) {
  load(file)  # loads final_df into environment
  
  dataset_name <- file_path_sans_ext(basename(file)) %>%
    sub("^final_df_tuned_", "", .)
  
  # Sanity check
  if (!exists("final_df")) next
  
  # Skip empty or invalid datasets
  if (nrow(final_df) == 0 || !"mtry_opt" %in% names(final_df)) next
  
  cov_maxnodes <- sd(final_df$best_lambda, na.rm = TRUE) / mean(final_df$best_lambda, na.rm = TRUE)
  median_tunability <- median(final_df$rel_tunability_bagging, na.rm = TRUE)
  
  result_list[[dataset_name]] <- data.frame(
    dataset = dataset_name,
    cov_maxnodes = cov_maxnodes,
    median_rel_tunability = median_tunability
  )
  
  rm(final_df)  # clean up for next load
}

#combine into one data frame
summary_df <- bind_rows(result_list)
lm_model <- lm(median_rel_tunability ~ cov_maxnodes, data = summary_df)
summary_lm <- summary(lm_model)
confint_lm <- confint(lm_model)

p_val <- signif(summary_lm$coefficients[2, 4], 3)
slope_ci <- confint_lm[2, ]  # lower and upper CI for slope
slope_ci_text <- sprintf("95%% CI: [%.3f, %.3f]", slope_ci[1], slope_ci[2])
p_val_text <- paste("p-value =", p_val)

#plot
p <- ggplot(summary_df, aes(x = cov_maxnodes, y = median_rel_tunability, label = dataset)) +
  geom_point(size = 3) +
  geom_text_repel(size = 6, max.overlaps = 100) +
  geom_smooth(method = "lm", se = TRUE, level = 0.95, color = "black", linetype = "dashed") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
           label = p_val_text,
           size = 5, fontface = "italic") +
  labs(
    #x = "CoV (maxnodes)",
    x = expression("CoV(" * lambda * ")"),
    y = "Median Relative Tunability (%)"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    plot.title = element_text(size = 20, hjust = 0.5)
  )

#ggsave("plots/trend_best_lambda.pdf", plot = p, width = 7, height = 5)
