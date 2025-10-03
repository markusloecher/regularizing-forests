library(future.apply)
library(randomForest)
library(ranger)

# Load grid search results
load("Rdata/MSE_high10.RData")

# Sim settings
n <- 100
nval <- 100
p <- 1000
s <- 10
rho <- 0.35
beta.type <- 2
ntrees <- 500
params <- c(0, 0.1, 1, 5, 10, 25, 50, 80, 95)
SNR <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6)
reps <- 1:100

# Plan parallel
plan(multisession, workers = 40)
#change workers as needed

results <- future_lapply(SNR, function(snr) {
  sim_res <- lapply(reps, function(rep) {
    set.seed(rep)
    xy <- sim.xy(n, nval, p, s, rho, beta.type, snr)
    x <- as.data.frame(xy$x)
    y <- xy$y
    xval <- as.data.frame(xy$xval)
    yval <- xy$yval
    sigma2 <- xy$sigma^2
    
    mse_mat <- MSE_out[[rep]][["MSE_test"]][[paste0("SNR ", snr)]] #MSE matrix
    mse_min <- min(mse_mat) #get the lowest MSE
    #Get the corresponding maxnodes, mtry
    opt_idx <- which(mse_mat == mse_min, arr.ind = TRUE)
    maxnodes_opt <- as.numeric(gsub("maxnode ", "", rownames(mse_mat)[opt_idx[,"row"]]))
    mtry_opt     <- as.numeric(gsub("mtry ", "", colnames(mse_mat)[opt_idx[,"col"]])) #returns proportion
    mtry_opt <- mtry_opt * p
    
    #fit a forest for those parameters
    rf_opt <- randomForest(x = x, y = y, 
                           xtest = xval, ytest = yval,
                           ntree = ntrees, 
                           seed = rep,
                           maxnodes = maxnodes_opt, mtry = mtry_opt)
    
    #set the obtained mse as optimal (best)
    mse_opt <- mean((rf_opt$test$predicted - yval)^2)
    rte_opt <- mse_opt / sigma2
    
    
    best_lambda <- NA
    mse_def <- Inf
    
    for (lambda in params) {
      #fresh unshrunk ranger rf
      rf_bag <- ranger(x = x, y = y, num.trees = ntrees, 
                       mtry = p, #bagging
                       min.node.size = 1, 
                       seed = rep, node.stats = TRUE, keep.inbag = TRUE)
      
      hshrink(rf_bag, lambda = lambda) #shrinkage
      preds <- predict(rf_bag, data = xval)$predictions
      mse <- mean((preds - yval)^2)
      if (mse < mse_def) {
        mse_def <- mse #baseline mse from best performing lambda
        best_lambda <- lambda
      }
    }
    
    rte_def <- mse_def / sigma2
    d_lambda <- rte_def - rte_opt #tunability
    
    data.frame(rep = rep, snr = snr,
               maxnodes_opt = maxnodes_opt,
               mtry_opt = mtry_opt,
               mse_opt = mse_opt,
               rte_opt = rte_opt,
               best_lambda = best_lambda,
               mse_def = mse_def,
               rte_def = rte_def,
               tunability = d_lambda)
  })
  cat("Done with SNR:", snr, "\n")
  do.call(rbind, sim_res)
}, future.seed = TRUE)

bagging_df <- do.call(rbind, results)
save(bagging_df, file = "Rdata/bagging_mse_parallel_high10.RData")
