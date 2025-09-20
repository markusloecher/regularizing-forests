generate_boxplots <- function(melted_df,
                              shrink_modes = c("hs", "hs_entropy","hs_permutation", "no_shrinkage"), 
                              yLimTop = 1,
                              yLimBottom = 0.0175,
                              to_facet = TRUE,
                              ylab='MDI'
) {
  # Load the necessary libraries
  library(ggplot2)
  library(reshape2)
  
  # Define a custom color palette with attractive colors
  custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6")
  #df = subset(df, df$shrink_mode %in% shrink_modes)
  #df$shrink_mode = gsub("hs_permutation", "hs_permute", df$shrink_mode)
  melted_df = subset(melted_df, melted_df$shrink_mode %in% shrink_modes)
  melted_df$shrink_mode = gsub("hs_permutation", "hs_permute", melted_df$shrink_mode)
  # Get the unique 'shrink_mode' values for the entire DataFrame
  unique_shrink_modes <- unique(melted_df$shrink_mode)
  
  # Define custom x-axis labels with LaTeX mathematical expressions
  x_labels <- c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]), expression(X[5]))
  
  # Define the height of the FacetGrid plots (moderately larger)
  plot_height <- 6
  p = list()
  #print(unique(melted_df$relevance))
  # Iterate through unique 'relevance' values and create attractive diagrams
  if (to_facet){
    p <- ggplot(melted_df, aes(x=Xvar, y=value, fill=shrink_mode)) +
      geom_boxplot(outlier.size=0.25, width = 0.95,lwd=0.2) +
      #scale_fill_brewer(palette = "Set2") +
      scale_fill_manual(values=custom_palette) +
      labs(x='', y=ylab) +
      theme_minimal() +
      coord_cartesian(ylim = c(yLimBottom, yLimTop)) +
      #theme(legend.position='right') +
      guides(fill=guide_legend(title='')) +
      # guides(fill="none") + 
      scale_x_discrete(labels=x_labels) +
      theme(
        legend.position = c(0.05, 1.075),  # top-left corner (x, y)
        legend.justification = c(0, 1),   # anchor top-left of the legend box
        legend.background = element_rect(fill = alpha('white', 0), color = NA),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 9)
       # legend.background = element_rect(fill = "transparent", color = NA)
      ) +
      theme(
        panel.spacing.x = unit(0., "lines"),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
      ) +
      scale_y_sqrt() +
      #ggtitle(paste('r=', relevance_value/100)) +
      #theme(plot.title = element_text(hjust = 0.5)) +
      facet_wrap(~ relevance, nrow = 1, labeller = label_bquote(r == .(relevance/100))) +
      theme(
        strip.text = element_text(size = 14),
        plot.title = element_text(hjust = 0.5)
      )
  } else {
    for (relevance_value in unique(melted_df$relevance)) {
      # Filter the data for the current 'relevance' value
      subset_df <- subset(melted_df, relevance == relevance_value)
      #return(subset_df)
      #cat(relevance_value, as.character(relevance_value), "\n")
      # Create a ggplot object for boxplots
      p[[as.character(relevance_value)]] <- ggplot(subset_df, aes(x=Xvar, y=value, fill=shrink_mode)) +
        geom_boxplot(outlier.size=0.25, width = 0.75) +
        scale_fill_manual(values=custom_palette) +
        labs(x='', y=ylab) +
        theme_minimal() +
        coord_cartesian(ylim = c(0, yLimTop)) +
        #theme(legend.position='right') +
        #guides(fill=guide_legend(title='')) +
        guides(fill="none") + 
        scale_x_discrete(labels=x_labels) +
        ggtitle(paste('r=', relevance_value/100)) +
        theme(plot.title = element_text(hjust = 0.5)) # Center the title
    }
    
    # Save the plot to a file with a relevant filename
    #ggsave(filename=paste('Relevance_', relevance_value, '_boxplot.png', sep=''), plot=p, width=12, height=plot_height)
  }
  return(p)
}





perf_boxplots <- function(df, shrink_modes = c("hs", "hs_entropy","hs_permutation", "no_shrinkage"),
                          fname = NULL, 
                          maximizeOverLambda = FALSE,
                          perf = "ROC.AUC",
                          title = "ROC",
                          verbose = 0
) {
  # Load the necessary libraries
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  # Define a custom color palette with attractive colors
  custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6")
  
  for (j in 1:length(df)){
    df[[j]] = subset(df[[j]], df[[j]]$shrink_mode %in% shrink_modes)
    df[[j]]$shrink_mode = gsub("hs_permutation", "hs_permute", df[[j]]$shrink_mode)
    df[[j]]$shrink_mode = gsub("no_shrinkage", "none", df[[j]]$shrink_mode)
    
    if (maximizeOverLambda){
      #browser()
      maxAUC = group_by(df[[j]],shrink_mode,lambda) %>% 
        summarise(avgAUC =mean(!!! rlang::syms(perf))) %>% 
        group_by(shrink_mode) %>%  
        summarise(maxAUC =max(avgAUC), opLamIndex=which.max(avgAUC), opLam = lambda[which.max(avgAUC)])
      #print(maxAUC)
      K=nrow(maxAUC)
      for (i in 1:K){
        x =subset(df[[j]], shrink_mode == maxAUC$shrink_mode[i] & lambda == maxAUC$opLam[i])
        if (i ==1) {
          bestROC = x
        } else {
          bestROC = rbind.data.frame(bestROC,x)
        }
      }
    } else {
      bestROC = df[[j]]
    }
    bestROC$dataset = names(df)[j]
    
    if (j ==1) {
      dfAll = bestROC 
    } else {
      dfAll = rbind.data.frame(dfAll,bestROC)
    }
  }#### end of loop over data sets
  if (verbose>1) return(dfAll)
  
  #maxAUC = group_by(dfAll,dataset, shrink_mode) %>% 
  #  summarise(avgAUC =mean(!!! rlang::syms(perf)))
  #print(maxAUC)
  
  if (perf == "R2"){
    p=ggplot(dfAll, aes(x=dataset, y=R2, fill=shrink_mode))
    ylab = expression(R^2)
  } else if (perf == "ROC.AUC") {
    p=ggplot(dfAll, aes(x=dataset, y=ROC.AUC, fill=shrink_mode))
    ylab = "AUC"
  }
  p=p + geom_boxplot(outlier.size=0.25) +
    scale_fill_manual(values=custom_palette) +
    labs(x='', y=ylab) +
    theme_minimal() +
    #coord_cartesian(ylim = c(0, )) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position='right') +
    guides(fill=guide_legend(title='shrink mode')) +
    #guides(fill="none") 
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) # Center the title
  
  
  if (!is.null(fname)) {
    pdf_filename = fname#paste0("Performance1.pdf")
    ggsave(pdf_filename, plot = p, device = "pdf", width = 8, height = 4, units = "in")
  }
  return(p)
}

# Example usage with your DataFrame 'your_df'
#generate_boxplots(df)

perf_boxplots_trees <- function(df, 
           shrink_modes = c("hs", "hs_entropy","hs_permutation", "no_shrinkage"),
           minTrees=5,
           perf = "ROC.AUC",
           addDataInfo2Title = TRUE,
           fname = NULL) {
  # Load the necessary libraries
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  
  # Define a custom color palette with attractive colors
  custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f")#, "#9b59b6")
  
  panels = list()##list of ggplots 
  for (j in 1:length(df)){
    df[[j]] = subset(df[[j]], df[[j]]$shrink_mode %in% shrink_modes & df[[j]]$num_trees >= minTrees)
    df[[j]]$shrink_mode = gsub("hs_permutation", "hs_permute", df[[j]]$shrink_mode)
    df[[j]]$shrink_mode = gsub("no_shrinkage", "none", df[[j]]$shrink_mode)
    ##df[[j]]$num_trees = factor(df[[j]]$num_trees)
    #browser()
    #THE FOLLOWING ASSUMES that no more maximization over lambda has 
    #to be performed by this function!
    meanAUC = group_by(df[[j]], shrink_mode,num_trees)  %>%  
      summarise(avgAUC =mean(!!! rlang::syms(perf)), #opLamIndex=which.max(avgAUC), opLam = lambda[which.max(avgAUC)],
                medAUC =median(!!! rlang::syms(perf)),
                sdAUC =sd(!!! rlang::syms(perf)), n = n())
    if (perf == "ROC.AUC"){
      ylab = "AUC"
      maxAUC0 = group_by(df[[j]],shrink_mode,lambda, num_trees) %>% 
        summarise(avgAUC =mean(ROC.AUC), sdAUC =sd(ROC.AUC), n = n())
    } else if (perf == "R2"){
      ylab = expression(R^2)
      maxAUC0 = group_by(df[[j]],shrink_mode,lambda, num_trees) %>% 
        summarise(avgAUC =mean(R2), sdAUC =sd(R2), n = n())
    }
    
    maxAUC = maxAUC0 %>% group_by(shrink_mode,num_trees)  %>%
      summarise(maxavgAUC =max(avgAUC), opLamIndex=which.max(avgAUC), opLam = lambda[which.max(avgAUC)],
                sd_AUC = sdAUC[which.max(avgAUC)],
                n = n[which.max(avgAUC)])
    
    # #print(maxAUC)
    title=names(df)[j]
    if (addDataInfo2Title){
      #rc =dim(df[[j]])
      if (perf == "ROC.AUC"){
        rc=list("breast-cancer" =c(277,17), "diabetes" =c(768,8), "german-credit" =c(1000,20), 
                "haberman" =c(306,3), "heart" =c(270,15), "ionosphere" =c(351,34), "juvenile" =c(3640,286), "recidivism" =c(6172,20))
      } else if (perf == "R2"){
        rc=list("abalone" =c(), "california-housing" =c(), "diabetes-reg" =c(), 
                "friedman1" =c(), "friedman2" =c(), "friedman3" =c(), "red-wine" =c(), "satellite-image" =c())
      }
      
      n= rc[[title]][1];p= rc[[title]][2]
      title = paste0(title, "\n (n=",n, ", p=", p,")")
    }
    #browser()
    panels[[j]] =ggplot(meanAUC, aes(x=num_trees, y=avgAUC)) +
      geom_line(aes(col = shrink_mode)) +
      scale_color_manual(values=custom_palette) +
      labs(x='', y='') +
      theme_minimal() +
      #coord_cartesian(ylim = c(0, )) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position="none") +scale_x_sqrt() +
      #guides(fill=guide_legend(title='shrink mode')) +
      guides(fill="none") +
      ggtitle(title) +
      theme(plot.title = element_text(size = 12, hjust = 0.5)) # Center the title
    
    names(panels)[j] = names(df)[j]
    cat("done with", names(df)[j] , "(", j, ")\n")
  }### end of loop over data sets
  
  
  return(panels)
}

# Example usage with your DataFrame 'your_df'
#perf_boxplots_trees(df)
