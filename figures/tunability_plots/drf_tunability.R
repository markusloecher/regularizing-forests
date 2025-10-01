library(dplyr)
library(tidyr)
library(ggplot2)

###FOR -tunability + d_RF to get rte_default - rte_bagginghs###

#load shrinkage data
load("Rdata/bagging_mse_parallel_low.RData")  #bagging_df

#load simulation data MSE_out
load("Rdata/MSE_low.RData")     # MSE_out

# Extract simulation parameters
temp <- MSE_out[[1]]
n <- temp[["n"]]
sigma <- temp[["sigma"]]
SNR <- temp[["SNR"]]
maxnode <- temp[["maxnode"]]
mtry <- temp[["mtry"]]
nrep <- length(MSE_out)

mtry0 <- 0.3
maxnode0 <- 1
mtry_dft <- which.min(abs(mtry - mtry0))
depth_dft <- which.min(abs(maxnode - maxnode0 * n))

# Compute d_RF
# From Zhou's function
d_RF <- matrix(NA, nrow = nrep, ncol = length(SNR))
for (k in seq_along(SNR)) {
  for (i in seq_len(nrep)) {
    RTE <- MSE_out[[i]]$MSE_test[[k]] / sigma[k]^2
    RTE0 <- min(RTE)
    d_RF[i, k] <- RTE[depth_dft, mtry_dft] - RTE0
  }
}

#using it as a dataframe 
df_d_rf <- as.data.frame(d_RF)
colnames(df_d_rf) <- paste0("snr_", SNR)
df_d_rf$rep <- seq_len(nrep)

df_d_rf_long <- pivot_longer(df_d_rf, 
                             cols = starts_with("snr_"),
                             names_to = "snr",
                             names_prefix = "snr_",
                             values_to = "d_rf")
df_d_rf_long$snr <- as.numeric(df_d_rf_long$snr)

#bagging_df has shrinkage+bagging mse data
df_combined <- left_join(bagging_df, df_d_rf_long, by = c("rep", "snr")) %>%
  mutate(tunability_total = - tunability + d_rf)

snr_levels <- c(0.05, 0.09, 0.14, 0.25, 0.42, 0.71, 1.22, 2.07, 3.52, 6.0)

df_summary <- df_combined %>%
  group_by(snr) %>%
  summarise(
    d = mean(tunability_total, na.rm = TRUE),
    se = sd(tunability_total, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

df_summary$snr <- factor(df_summary$snr, levels = snr_levels)

ggplot(df_summary, aes(x = snr, y = d)) +
  geom_point(color = "black", size = 2) +
  geom_line(aes(group = 1), color = "black", linewidth = 1) +
  geom_errorbar(aes(ymin = d - se, ymax = d + se), width = 0.1) +
  xlab("SNR") +
  ylab("Tunability (Low)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_text(size = rel(1.2)),
    axis.text = element_text(size = rel(1.1))
  )

save(df_combined, file = "Rdata/tunability_rf_low.RData")
