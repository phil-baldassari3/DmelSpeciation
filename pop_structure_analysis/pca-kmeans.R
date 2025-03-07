library(ggplot2)
library(tidyverse)

setwd("/Users/philipbaldassari/Desktop/dros_speciation/Pop_Structure_Analysis/PCA_csvs/PCA_kmeans")

# SET THESE
plot_title <- "PCA on Autosomal SNPs with K-means Clusters (minor allele frequency \u2265 5%, k=5)"
pca_k_csv <- "PCA_maf0.05_Autosomes_k5.csv"


# Read the data
data <- read.csv(pca_k_csv)

#color pallet
my_colors <- c(
  "blue", "cyan", "burlywood", "blueviolet", "indianred1", "green3", "maroon2", 
  "black", "darkgreen", "magenta", "ivory3", "khaki2", "gold4", "brown", 
  "yellow", "plum1", "chocolate", "chartreuse", "skyblue", "red", "orange")



# Combined visualization
ggplot(data, aes(EV2, EV1)) +
  geom_point(aes(color = as.factor(Population), shape = as.factor(cluster)), size = 3, alpha=I(0.7)) +
  labs(title = plot_title, color = "Population", shape = "Cluster") +
  scale_color_manual(values = my_colors) +
  theme_bw()


