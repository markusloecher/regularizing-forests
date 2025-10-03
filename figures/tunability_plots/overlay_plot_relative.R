library(reticulate)
library(dplyr)
library(ggplot2)

#Load RData

load("Rdata/MSE_low.RData")  #MSE_out
temp      <- MSE_out[[1]]
n         <- temp[['n']]
sigma     <- temp[['sigma']]
SNR       <- temp[['SNR']]
maxnode   <- temp[['maxnode']]
mtry      <- temp[['mtry']]
nrep      <- length(MSE_out)
l_SNR     <- length(SNR)
mtry0     <- 0.3
maxnode0  <- 1
mtry_dft  <- which.min(abs(mtry - mtry0))
depth_dft <- which.min(abs(maxnode - maxnode0 * n))
baseline  <- "max"


#compute relative tunability from RData
d_RF <- d_mtry <- d_depth <- d_additional <- matrix(NA, nrep, l_SNR)
opt_RF <- matrix(NA, nrep, 2 * l_SNR)
colnames(opt_RF) <- rep(c("maxnodes", "mtry"), l_SNR)
#store default rf values
rte_default_rf <- matrix(NA, nrow = nrep, ncol = l_SNR)

for (k in 1:l_SNR) {
  for (i in 1:nrep) {
    RTE <- MSE_out[[i]]$MSE_test[[k]] / sigma[k]^2
    RTE0 <- min(RTE)
    ind_opt_RF <- which(RTE == RTE0, arr.ind = TRUE)
    opt_RF[i, (2 * k - 1):(2 * k)] <- c(maxnode[ind_opt_RF[1]] / n, mtry[ind_opt_RF[2]])
    rte_default_rf[i, k] <- RTE[depth_dft, mtry_dft]
    
    if (baseline == "max") {
      d_RF[i, k]    <- 100 * (RTE[depth_dft, mtry_dft] - RTE0) / RTE[depth_dft, mtry_dft]
      d_mtry[i, k]  <- 100 * (RTE[depth_dft, mtry_dft] - min(RTE[depth_dft, ])) / RTE[depth_dft, mtry_dft]
      d_depth[i, k] <- 100 * (RTE[depth_dft, mtry_dft] - min(RTE[, mtry_dft])) / RTE[depth_dft, mtry_dft]
    } else {
      stop("Only 'max' baseline implemented in this block.")
    }
  }
}

se <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

d_RF_ave    <- colMeans(d_RF)
d_RF_se     <- apply(d_RF, 2, se)
d_mtry_ave  <- colMeans(d_mtry)
d_mtry_se   <- apply(d_mtry, 2, se)
d_depth_ave <- colMeans(d_depth)
d_depth_se  <- apply(d_depth, 2, se)
#d_additional_ave <- colMeans(d_additional)
#d_additional_se  <- apply(d_additional, 2, se)

dat_d_ave <- data.frame(
  d = c(d_RF_ave, 
        d_mtry_ave, d_depth_ave 
        #d_additional_ave
        ),
  se = c(d_RF_se, d_mtry_se, d_depth_se
         #d_additional_se
         ),
  SNR = rep(SNR, times = 3),
  Type = rep(c("RF", "mtry", "maxnodes"
               #"Additional"
               ), each = l_SNR),
  mtry_dft = mtry[mtry_dft],
  depth_dft = maxnode[depth_dft],
  alpha = 0.4
)


#default rf RTE values as a df
df_rte_dft <- as.data.frame(rte_default_rf)
colnames(df_rte_dft) <- paste0("snr_", SNR)
df_rte_dft$rep <- seq_len(nrep)

df_rte_dft_long <- pivot_longer(df_rte_dft,
                                cols = starts_with("snr_"),
                                names_to = "snr",
                                names_prefix = "snr_",
                                values_to = "rte_dft_rf") %>%
  mutate(snr = as.numeric(snr))


#load shrinkage pickle file and compute relative tunability
pd <- import("pandas")
df_py <- py_to_r(pd$read_pickle("Rdata/hs_mse_low_parallel.pkl"))

#join default rf rte values
df_py <- left_join(df_py, df_rte_dft_long, by = c("rep", "snr"))

shrinkage_summary <- df_py %>%
  arrange(rte_best) %>%
  group_by(snr, rep) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(rel_tunability = 100 * (rte_def - rte_best) / rte_def) %>%
  group_by(snr) %>%
  summarise(d = mean(rel_tunability, na.rm = TRUE),
            se = sd(rel_tunability, na.rm = TRUE)/sqrt(n()),
            .groups = "drop") %>%
  mutate(Type = "Shrinkage",
         mtry_dft = NA,
         depth_dft = NA,
         SNR = snr,
         alpha = 1) %>%
  select(SNR, d, se, Type, mtry_dft, depth_dft, alpha)


#load saved Bagging + Shrinkage RData

load("Rdata/tunability_rf_low.RData")  # loads df_combined

#join rf_dft values
df_combined <- left_join(df_combined, df_rte_dft_long, by = c("rep", "snr"))

df_combined <- df_combined %>%
  #mutate(rel_tunability = 100 * (tunability_total / rte_def)) %>% 
  mutate(rel_tunability = 100 * (tunability_total / rte_dft_rf)) %>%
  group_by(snr) %>%
  summarise(d = mean(rel_tunability, na.rm = TRUE),
            se = sd(rel_tunability, na.rm = TRUE)/sqrt(n()),
            .groups = "drop") %>%
  mutate(Type = "Bagging + Shrinkage",
         mtry_dft = NA,
         depth_dft = NA,
         SNR = snr,
         alpha = 1) %>%
  select(SNR, d, se, Type, mtry_dft, depth_dft, alpha)


#combine everything and plot
combined_df <- rbind(dat_d_ave, 
                     shrinkage_summary, 
                     df_combined
                     #df_bag_shrink_summary
                     )

gp_rel_combined <- ggplot(combined_df, aes(x = SNR, y = d, group = Type)) +
  geom_point(aes(color = Type, alpha = alpha)) +
  geom_line(aes(color = Type, alpha = alpha), linewidth = 1) +
  geom_errorbar(aes(ymin = d - se, ymax = d + se, color = Type, alpha = alpha), width = 0.1) +
  scale_color_manual(values = c("RF" = "#C77CFF",
                                "mtry" = "#00FFFF",
                                "maxnodes" = "#00BA38",
                                #"Additional" = "#F8766D",
                                "Shrinkage" = "black",
                                "Bagging + Shrinkage" = "red")) +
  scale_alpha_identity() +
  scale_x_continuous(breaks = unique(combined_df$SNR), trans = "log") +
  xlab("Signal-to-Noise Ratio") +
  ylab("Relative Tunability (%)") +
  ggtitle("Low") + #change as needed
  theme_bw() +
  theme(
    #legend.position = c(0.02, 0.98),
    legend.position = "none",
    legend.justification = c(0, 1),
    plot.title = element_text(hjust = 0.5, size = rel(2.7), face = 'bold'),
    legend.key.size = unit(0.04, "npc"),
    legend.spacing.y = unit(.01, "npc"),
    legend.text = element_text(size = rel(2.3)),
    legend.title = element_blank(),
    axis.title = element_text(hjust = 0.5, size = rel(2.2)),
    axis.text = element_text(size = rel(1.5))
  )

#Save as needed
#ggsave("plots/overlay_plot_relative_low.pdf", gp_rel_combined, width = 8, height = 6)