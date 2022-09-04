#loading packages
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

#set directory
setwd("/Users/philipbaldassari/Desktop/p-value_heatmap")

#opening dfs
zim_cos_control_pvalues <- read.csv("zim_cos_control_p-values.csv")
zim_cos_control_sites_pvalues <- read.csv("zim_cos_control_sites_p-values.csv")
zim_cos_control_genes_pvalues <- read.csv("zim_cos_control_genes_p-values.csv")

#creating breaks for significance ranges
FB_zim_cos_control_pvalues <- zim_cos_control_pvalues
FB_zim_cos_control_pvalues$p_value <- cut(FB_zim_cos_control_pvalues$FB_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))
head(FB_zim_cos_control_pvalues)

FM_zim_cos_control_pvalues <- zim_cos_control_pvalues
FM_zim_cos_control_pvalues$p_value <- cut(FM_zim_cos_control_pvalues$FM_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))


FB_zim_cos_control_sites_pvalues <- zim_cos_control_sites_pvalues
FB_zim_cos_control_sites_pvalues$p_value <- cut(FB_zim_cos_control_sites_pvalues$FB_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))

FM_zim_cos_control_sites_pvalues <- zim_cos_control_sites_pvalues
FM_zim_cos_control_sites_pvalues$p_value <- cut(FM_zim_cos_control_sites_pvalues$FM_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))


FB_zim_cos_control_genes_pvalues <- zim_cos_control_genes_pvalues
FB_zim_cos_control_genes_pvalues$p_value <- cut(FB_zim_cos_control_genes_pvalues$FB_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))

FM_zim_cos_control_genes_pvalues <- zim_cos_control_genes_pvalues
FM_zim_cos_control_genes_pvalues$p_value <- cut(FM_zim_cos_control_genes_pvalues$FM_p, breaks = c(0, 0.01, 0.05, 0.1, 0.5, 1))

#ordering plots
comparisons <- unique(FB_zim_cos_control_pvalues$comparison)
GO_terms <- rev(unique(FB_zim_cos_control_pvalues$GO_term))
comparisons_sites <- unique(FB_zim_cos_control_sites_pvalues$comparison)
GO_terms_sites <- rev(unique(FB_zim_cos_control_sites_pvalues$GO_term))
comparisons_genes <- unique(FB_zim_cos_control_genes_pvalues$comparison)
GO_terms_genes <- rev(unique(FB_zim_cos_control_genes_pvalues$GO_term))


#plotting
FB_zim_cos_control_pvalues_plot <- ggplot(FB_zim_cos_control_pvalues, aes(x=factor(comparison, level=comparisons), y=factor(GO_term, level=GO_terms), fill = p_value)) + 
  ggtitle("Flybase GO Enrichment p-values") + 
  geom_tile() + scale_fill_manual(breaks = levels(FB_zim_cos_control_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")

FM_zim_cos_control_pvalues_plot <- ggplot(FM_zim_cos_control_pvalues, aes(x=factor(comparison, level=comparisons), y=factor(GO_term, level=GO_terms), fill = p_value)) + 
  ggtitle("Flymine GO Enrichment p-values") + 
  geom_tile() + scale_fill_manual(breaks = levels(FM_zim_cos_control_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")



FB_zim_cos_control_sites_pvalues_plot <- ggplot(FB_zim_cos_control_sites_pvalues, aes(x=factor(comparison, level=comparisons_sites), y=factor(GO_term, level=GO_terms_sites), fill = p_value)) + 
  ggtitle("Flybase GO Enrichment p-values on sites") + 
  geom_tile() + scale_fill_manual(breaks = levels(FB_zim_cos_control_sites_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")

FM_zim_cos_control_sites_pvalues_plot <- ggplot(FM_zim_cos_control_sites_pvalues, aes(x=factor(comparison, level=comparisons_sites), y=factor(GO_term, level=GO_terms_sites), fill = p_value)) + 
  ggtitle("Flymine GO Enrichment p-values on sites") + 
  geom_tile() + scale_fill_manual(breaks = levels(FM_zim_cos_control_sites_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")



FB_zim_cos_control_genes_pvalues_plot <- ggplot(FB_zim_cos_control_genes_pvalues, aes(x=factor(comparison, level=comparisons_genes), y=factor(GO_term, level=GO_terms_genes), fill = p_value)) + 
  ggtitle("Flybase GO Enrichment p-values on genes") + 
  geom_tile() + scale_fill_manual(breaks = levels(FB_zim_cos_control_genes_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")

FM_zim_cos_control_genes_pvalues_plot <- ggplot(FM_zim_cos_control_genes_pvalues, aes(x=factor(comparison, level=comparisons_genes), y=factor(GO_term, level=GO_terms_genes), fill = p_value)) + 
  ggtitle("Flymine GO Enrichment p-values on genes") + 
  geom_tile() + scale_fill_manual(breaks = levels(FM_zim_cos_control_genes_pvalues$p_value), values = c("red", "indianred2", "pink4", "grey60", "grey40")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=14), axis.text.y = element_text(size=14)) +
  labs(x="Pairwise Fst Comparison", y="GO term")



#saving plots
png(filename = "FB_zim_cos_control_pvalues_plot.png", width = 1000, height = 700)
FB_zim_cos_control_pvalues_plot
dev.off()


png(filename = "FM_zim_cos_control_pvalues_plot.png", width = 1000, height = 700)
FM_zim_cos_control_pvalues_plot
dev.off()


png(filename = "FB_zim_cos_control_sites_pvalues_plot.png", width = 1000, height = 700)
FB_zim_cos_control_sites_pvalues_plot
dev.off()


png(filename = "FM_zim_cos_control_sites_pvalues_plot.png", width = 1000, height = 700)
FM_zim_cos_control_sites_pvalues_plot
dev.off()


png(filename = "FB_zim_cos_control_genes_pvalues_plot.png", width = 1000, height = 700)
FB_zim_cos_control_genes_pvalues_plot
dev.off()


png(filename = "FM_zim_cos_control_genes_pvalues_plot.png", width = 1000, height = 700)
FM_zim_cos_control_genes_pvalues_plot
dev.off()









