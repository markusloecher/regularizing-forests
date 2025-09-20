# library
library(ggplot2)
library(dplyr)

# Random Forests

fname_rf = "/Users/loecherm/Downloads/wetransfer_128-replications_2023-09-15_0526/strobl/strobl_rf/scores.csv"
df2 <- read.csv(fname_rf)

nrow(df2)#22400
length(unique(df2$lambda))*length(unique(df2$shrink_mode))*length(unique(df2$replication))*length(unique(df2$relevance))

#find optimal lambda first:
maxAUC = group_by(df2,shrink_mode,lambda) %>% 
  summarise(avgAUC =mean(ROC.AUC)) #%>% #medAUC =median(ROC.AUC))
  group_by(shrink_mode) %>%  
  summarise(maxAUC =max(avgAUC), opLam = which.max(avgAUC))

df2 = subset(df2, lambda %in% c(100,50))
df2 = subset(df2, shrink_mode %in% c("hs", "hs_entropy", "hs_permutation"))

df2$relevance = factor(df2$relevance/100)
# grouped boxplot
# Define a custom color palette with attractive colors
custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6")

gp = ggplot(df2, aes(x=relevance, y=ROC.AUC, fill=shrink_mode)) + 
  geom_boxplot(outlier.size=0.5) +
  scale_fill_manual(values=custom_palette) +
  labs(x='relevance', y='AUC') +
  theme_minimal()
gp

ggsave(filename=paste('figs/AUC_strobl_boxplot.pdf', sep=''), device="pdf", plot=gp, width=5, height=3)


################## Single Decision Tree

fname_dt = "/Users/loecherm/Downloads/wetransfer_128-replications_2023-09-15_0526/strobl/strobl_dt/scores.csv"
df2 <- read.csv(fname_dt)

#find optimal lambda first:
maxAUC = group_by(df2,shrink_mode,lambda) %>% 
  summarise(avgAUC =mean(ROC.AUC)) %>% #medAUC =median(ROC.AUC))
group_by(shrink_mode) %>%  
  summarise(maxAUC =max(avgAUC), opLam = which.max(avgAUC))


df2 = subset(df2, lambda %in% c(100,50))
df2 = subset(df2, shrink_mode %in% c("hs", "hs_entropy", "hs_permutation"))

df2$relevance = factor(df2$relevance/100)
# grouped boxplot
# Define a custom color palette with attractive colors
custom_palette <- c("#3498db", "#e74c3c", "#2ecc71", "#f1c40f", "#9b59b6")

gp = ggplot(df2, aes(x=relevance, y=ROC.AUC, fill=shrink_mode)) + 
  geom_boxplot() +
  scale_fill_manual(values=custom_palette) +
  labs(x='relevance', y='AUC') +
  theme_minimal()
gp

ggsave(filename=paste('figs/AUC_strobl_boxplot_dt.pdf', sep=''), device="pdf", plot=gp, width=5, height=3)
