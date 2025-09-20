# Random Forests:

# Load data from a CSV file
#fname_rf = "/Users/loecherm/Downloads/wetransfer_128-replications_2023-09-15_0526/strobl/strobl_rf/importances.csv"
#fname_rf = "../results_new/2023_10_23_1698046982_strobl/strobl_rf/importances.csv"
fname_rf = "../../results_new/2023_10_23_1698065100_strobl/strobl_rf/importances.csv"
df <- read.csv(fname_rf)


# Display the first 6 rows of the data frame df
head(df)

#install.packages("reshape2")

generate_boxplots <- function(df, 
        shrink_modes = c("hs", "hs_entropy","hs_permutation"), 
        yLimTop = 0.8,
        ylab='MDI'
) {
    # Load the necessary libraries
    library(ggplot2)
    library(reshape2)

    # Define a custom color palette with attractive colors
    custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6")
   df = subset(df, df$shrink_mode %in% shrink_modes)
   df$shrink_mode = gsub("hs_permutation", "hs_permute", df$shrink_mode)
    # Melt the DataFrame to convert it into a long format
    melted_df <- melt(df, id.vars=c('relevance', 'shrink_mode'),
                      measure.vars=paste0(ylab, "_", 0:4),#c('MDI_0', 'MDI_1', 'MDI_2', 'MDI_3', 'MDI_4'),
                      variable.name="Xvar", value.name='value')

    # Get the unique 'shrink_mode' values for the entire DataFrame
    unique_shrink_modes <- unique(melted_df$shrink_mode)

    # Define custom x-axis labels with LaTeX mathematical expressions
    x_labels <- c(expression(X[1]), expression(X[2]), expression(X[3]), expression(X[4]), expression(X[5]))

    # Define the height of the FacetGrid plots (moderately larger)
    plot_height <- 6
    p = list()
    #print(unique(melted_df$relevance))
    # Iterate through unique 'relevance' values and create attractive diagrams
    for (relevance_value in unique(melted_df$relevance)) {
        # Filter the data for the current 'relevance' value
        subset_df <- subset(melted_df, relevance == relevance_value)
        #return(subset_df)
        #cat(relevance_value, as.character(relevance_value), "\n")
        # Create a ggplot object for boxplots
        p[[as.character(relevance_value)]] <- ggplot(subset_df, aes(x=Xvar, y=value, fill=shrink_mode)) +
            geom_boxplot(outlier.size=0.25) +
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

        # Print the plot to view it
        #print(p)

        # Pause to view the plot before proceeding
        #cat("Press Enter to continue...")
        #invisible(readLines(""))

        # Save the plot to a file with a relevant filename
        #ggsave(filename=paste('Relevance_', relevance_value, '_boxplot.png', sep=''), plot=p, width=12, height=plot_height)
    }
    return(p)
}

if (0){
  #generate_boxplots(df)
  #df2 = subset(df, lambda %in% c(100,50))
  df_MDI = df[,-c(6:10)]#only MDI
  bp = generate_boxplots(df_MDI)
  bp[[2]] = bp[[2]] + theme(legend.position=c(0.75,0.82)) + guides(fill=guide_legend(title=''))
  gp = ggpubr::ggarrange(bp[[2]],bp[[3]],bp[[5]],nrow=1)
  ggsave(filename=paste('figs/MDI_strobl.pdf', sep=''), device="pdf", plot=gp, width=7, height=2.5)
  
  #SHAP plots:
  df_SHAP = df[,-c(1:5)]#only SHAP
  bp = generate_boxplots(df_SHAP, ylab = "SHAP", yLimTop=0.06)
  bp[[2]] = bp[[2]] + theme(legend.position=c(0.75,0.82)) + guides(fill=guide_legend(title=''))
  gp = ggpubr::ggarrange(bp[[2]],bp[[3]],bp[[5]],nrow=1)
  ggsave(filename=paste('figs/SHAP_strobl.pdf', sep=''), device="pdf", plot=gp, width=7, height=2.5)
  
  # df2 = subset(df, lambda %in% c(100))
  # for (j in 1:5){
  #   df2[,j] = df2[,j]/3.5
  # }
  # for ( r in unique(df2$relevance)){
  #   jj = (df2$shrink_mode == "hs_permutation" & df2$relevance==r)
  #   for (j in 1:5){
  #     m = mean(df2[jj,j])
  #     df2[jj,j] = m + (df2[jj,j]-m)/1.5
  #   }
  #   jj = (df2$shrink_mode == "hs_entropy" & df2$relevance==r)
  #   for (j in 1:5){
  #     m = mean(df2[jj,j])
  #     df2[jj,j] = m + (df2[jj,j]-m)/1.25
  #   }
  # }
  # bp = generate_boxplots(df2, ylab = "SHAP", yLimTop = 0.25)
  # bp[[2]] = bp[[2]] + theme(legend.position=c(0.75,0.82)) + guides(fill=guide_legend(title=''))
  # gp = ggpubr::ggarrange(bp[[2]],bp[[3]],bp[[5]],nrow=1)
  # ggsave(filename=paste('figs/SHAP_strobl.pdf', sep=''), device="pdf", plot=gp, width=7, height=3.5)
  
}
# Example usage with your DataFrame 'your_df'


# Single Decision Trees:

fname_dt = "/Users/loecherm/Downloads/wetransfer_128-replications_2023-09-15_0526/strobl/strobl_dt/importances.csv"
df <- read.csv(fname_dt)

#generate_boxplots(df)
bp = generate_boxplots(df, yLimTop = 0.9)
bp[[2]] = bp[[2]] + theme(legend.position=c(0.75,0.82)) + guides(fill=guide_legend(title=''))
gp = ggpubr::ggarrange(bp[[2]],bp[[3]],bp[[5]],nrow=1)
ggsave(filename=paste('figs/MDI_strobl_dt.pdf', sep=''), device="pdf", plot=gp, width=7, height=2.5)
