enlist <- function (...) 
{
  result <- list(...)
  if ((nargs() == 1) & is.character(n <- result[[1]])) {
    result <- as.list(seq(n))
    names(result) <- n
    for (i in n) result[[i]] <- get(i)
  }
  else {
    n <- sys.call()
    n <- as.character(n)[-1]
    if (!is.null(n2 <- names(result))) {
      which <- n2 != ""
      n[which] <- n2[which]
    }
    names(result) <- n
  }
  result
}


se <- function(x, margin =2){
  
  if (margin ==2) {
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(colSums(!is.na(x)))
  }else if(margin ==1){
    out <- apply(x, margin, sd, na.rm = T)/
      sqrt(rowSums(!is.na(x)))
  }
  
  
  return(out)
}


#library(randomForest)
library(tidyverse)
library(ggplot2)
library(grid)
#library(gridExtra)

# load('data/final_df_tuned_pah1.RData')
# bagging_result <- final_df


summary_and_boxplot_real <- function(FinalResult, 
                                     dt_name,
                                     mtry0=0.3, maxnodes0=1,
                                     lgd_pos = c(0.02, 0.98),
                                     ymin_box = NULL,
                                     ymax_box = NULL,
                                     width = 5.5,
                                     height = 7,
                                     device = "pdf",
                                     save = T){
  
  n <- FinalResult[[1]]$n
  p <- FinalResult[[1]]$p
  K <- FinalResult[[1]]$K
  dat <- FinalResult[[1]]$dat
  
  mtry <- FinalResult[[1]]$mtry
  maxnodes <- FinalResult[[1]]$maxnodes
  var_y <- FinalResult[[1]]$dat[,(p+1)] %>% var()
  nrep <- length(FinalResult)
  
  l_maxnodes <- length(maxnodes)
  l_mtry <- length(mtry)
  
  # Calculate the normalizing constants for tunability
  normal_cst <- var_y
  
  # Initialize lists for tunability d
  d_RF <- d_mtry <- d_depth <- d_add <- numeric(nrep)
  
  
  opt_RF  <- matrix(NA, nrow = nrep, ncol = 2)
  colnames(opt_RF) <- c("maxnodes", "mtry")
  
  # Index of default values of parameters
  mtry_dft <- which.min(abs(mtry - ceiling(mtry0*p)))
  depth_dft <- which.min(abs(maxnodes - maxnodes0*n))
  
  for (i in 1:nrep) {
    RTE     <- FinalResult[[i]]$cv.sigma/normal_cst
    RTE_dft <- RTE[mtry_dft, depth_dft]
    
    # optimal combination of mtry and maxnodes
    ind_opt_RF <- which(RTE == min(RTE), arr.ind = T)
    opt_RF[i, ] <- c(maxnodes[ind_opt_RF[2]]/n, mtry[ind_opt_RF[1]]/p)
    
    # d_RF[i]    <- RTE_dft - min(RTE)
    # d_mtry[i]  <- RTE_dft - min(RTE[, depth_dft])
    # d_depth[i] <- RTE_dft - min(RTE[mtry_dft, ])
    # d_add[i]   <- min(c(min(RTE[, depth_dft]),
    #                     min(RTE[mtry_dft, ]))) - min(RTE)
    
    RTE_best <- min(RTE)
    
    # Relative tunability as percent gain over optimal configuration
    d_RF[i]    <- 100 * (RTE_dft - min(RTE)) / RTE_dft
    d_mtry[i]  <- 100 * (RTE_dft - min(RTE[, depth_dft])) / RTE_dft
    d_depth[i] <- 100 * (RTE_dft - min(RTE[mtry_dft, ])) / RTE_dft
  }    
  
  
  opt_RF_med  <- apply(opt_RF, 2, median, na.rm=TRUE)%>% matrix(nrow = 2)
  opt_RF_mad   <- (apply(opt_RF, 2, mad, na.rm=TRUE) /
                     sqrt(colSums(!is.na(opt_RF)))) %>% matrix(nrow = 2)
  rownames(opt_RF_med) <- rownames(opt_RF_mad) <- c("maxnodes","mtry")
  
  # # Generate boxplots of tunability d
  # dat_d <- data.frame(d         = c(d_RF    %>% as.numeric(),
  #                                   d_mtry  %>% as.numeric(),
  #                                   d_depth %>% as.numeric(),
  #                                   d_add   %>% as.numeric()),
  #                     Type      = c("RF", "mtry", "maxnodes", "Additional") %>% 
  #                       rep(each = nrep) %>% as.factor(),
  #                     mtry_dft  = mtry[mtry_dft] %>% as.factor(),
  #                     depth_dft = maxnodes[depth_dft] %>% as.factor(),
  #                     dataset   = dt_name %>% as.factor()
  # )
  
  d_bagging <- bagging_result$rel_tunability_bagging

  # After computing d_RF, d_mtry, d_depth, d_bagging:
  if (is.null(ymin_box)) ymin_box <- min(c(d_RF, d_mtry, d_depth, d_bagging), na.rm = TRUE)
  if (is.null(ymax_box)) ymax_box <- max(c(d_RF, d_mtry, d_depth, d_bagging), na.rm = TRUE)
  
  dat_d <- data.frame(d         = c(d_RF    %>% as.numeric(),
                                    d_mtry  %>% as.numeric(),
                                    d_depth %>% as.numeric(),
                                    d_bagging),
                      Type      = c("RF", "mtry", "maxnodes", "BaggingHS") %>% 
                        rep(each = nrep) %>% as.factor(),
                      mtry_dft  = mtry[mtry_dft] %>% as.factor(),
                      depth_dft = maxnodes[depth_dft] %>% as.factor(),
                      dataset   = dt_name %>% as.factor()
  )
  
  
  gp_d <- ggplot(dat_d, aes(x = Type, y = d)) +
    geom_boxplot(aes(color = Type)) +
    xlab("Type") +
    ylab("Relative Tunability(%)") + 
    #ggtitle(dt_name)
    ggtitle(sub("^MSE_out_", "", tools::file_path_sans_ext(basename(dt_name))))
  
  
  
  if (!is.null(ymax_box)) {
    #gp_d <- gp_d + coord_cartesian(ylim = c(0,ymax_box)) 
    gp_d <- gp_d + scale_y_continuous(limits = c(ymin_box, ymax_box))
  }
  
  gp_d <- gp_d + 
    theme_bw() +
    theme(legend.position = "None",
          legend.justification = c(0,1),
          plot.title = element_text(hjust = 0.5, size = 26, face = 'bold'),
          legend.key.size = unit(0.04, "npc"),
          legend.spacing.y = unit(.01, "npc"),
          legend.text = element_text(size = rel(1.5)),
          legend.title = element_blank(),
          axis.title = element_text(hjust = 0.5,size = 20),
          axis.text  = element_text(size = 15)
    )
  # gp_d
  dataname <- sub("^MSE_out_", "", tools::file_path_sans_ext(basename(dt_name)))
  filenm <- file.path("plots", paste0("RelTune_", basename(dataname), ".pdf"))
  #filenm <- paste0("Tunability_", dt_name,"_boxplot.png")
  if (save) {
    ggsave(filenm, plot =gp_d, width = width, height = height, device = device)
  }
  
  
  dat_opt <- data.frame(opt = opt_RF_med %>% as.numeric(),
                        mad = opt_RF_mad %>% as.numeric(),
                        dataset = dt_name,
                        Type = c("maxnodes","mtry") %>% as.factor()
  )
  
  
  
  result <- enlist(dt_name, mtry0, maxnodes0, 
                   dat_d, dat_opt, 
                   gp_d
  )
  
  return(result)
  
}

# load("data/MSE_out_fb.RData")
# 
# summary_and_boxplot_real(
#   FinalResult = MSE_out,
#   dt_name = "data/MSE_out_fb.RData",
#   mtry0 = 0.3,
#   maxnodes0 = 1,
#   save = TRUE
# )

data_files <- c("bike", "boston", "concrete", "cpu", "csm", "fb", "servo", "solar",
                 "AquaticTox", "pah", "pdgfr", "phen")

#data_files <- c("concrete")

for (data_name in data_files) {
  
  # Construct file paths
  path_final_df <- paste0("data/final_df_tuned_", data_name, ".RData")
  path_mse_out  <- paste0("data/MSE_out_", data_name, ".RData")
  
  # Load bagging + shrinkage results
  load(path_final_df)
  bagging_result <- final_df  # assign inside global environment for function use
  
  # Load grid search results
  load(path_mse_out)
  
  # Run the summary + boxplot function
  summary_and_boxplot_real(
    FinalResult = MSE_out,
    dt_name = path_mse_out,
    mtry0 = 0.3,
    maxnodes0 = 1,
    save = TRUE
  )
}



#
#Histogram distributions of parameters 
library(ggplot2)

# Load data
load("data/final_df_tuned_solar.RData") #change as needed

# Create ggplot histogram
p <- ggplot(final_df, aes(x = maxnodes_opt)) +
  geom_histogram(fill = "#00BA38", color = "black", bins = 9) +
  ggtitle("Optimal Tree Depth (solar)") +
  xlab("maxnodes") +
  ylab("Frequency") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 26),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text  = element_text(size = 15)
  )

# Save as PDF
ggsave("plots/hist_maxnodes_opt_solar.pdf", plot = p, width = 7, height = 5)
