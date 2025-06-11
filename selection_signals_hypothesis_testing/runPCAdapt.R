library(pcadapt)
library(tidyverse)

#setting directory
setwd("~/Desktop/dros_speciation/selection_signals_hypothesis_testing/PCAdapt")

#reading files
bim <- read.table("bed_for_PCAdapt/biallelicSNPs_NoSingletons_missing0.05_dmel_r6_allChr_202samples.bim")
bed_file <- "bed_for_PCAdapt/biallelicSNPs_NoSingletons_missing0.05_dmel_r6_allChr_202samples.bed"
pcadapt_file <- read.pcadapt(bed_file, type = "bed")

#setting populations
poplist.names <- c(rep("FR", 40), rep("RAL", 84), rep("ZH", 14), rep("ZI", 51), rep("ZS", 13))

#running PCAdapt (Note: default maf threshold is 0.05)
x <- pcadapt(input = pcadapt_file, K = 20, LD.clumping = list(size = 500, thr = 0.1)) #, min.maf = 0.01


#finding K
plot(x, option = "screeplot")
plot(x, option = "screeplot", K = 10)
plot(x, option = "scores", pop = poplist.names)
plot(x, option = "scores", pop = poplist.names, i = 2, j = 3)
plot(x, option = "scores", pop = poplist.names, i = 3, j = 4)
plot(x, option = "scores", pop = poplist.names, i = 4, j = 5)

### K = 2 ###

#running pcadapt with k=4
x <- pcadapt(input = pcadapt_file, K = 2, LD.clumping = list(size = 500, thr = 0.1))  #, min.maf = 0.01
summary(x)
summary(x$pvalues)


#QQ
plot(x, option = "qqplot")

#assessing
par(mfrow = c(1, 2))
for (i in 1:2)
  plot(x$loadings[, i], pch = 19, cex = .3, ylab = paste0("Loadings PC", i))



#results table
results <- data.frame(
  CHROM = bim$V1[x$pass],
  POS   = bim$V4[x$pass],
  PVAL  = x$pvalues[x$pass],
  QVAL  = p.adjust(x$pvalues[x$pass], method = "BH")
)
results$CHROM <- as.character(results$CHROM)
results$CHROM[results$CHROM == "23"] <- "X"

#distribution of adjusted p-values
hist(results$QVAL, 
     breaks = 50, 
     col = "orange", 
     border = "black",
     main = "Distribution of Adjusted P-values",
     xlab = "P-values",
     ylab = "Number of SNPs")
abline(v = 0.05, col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = "p = 0.05", box.col = "white", bg = rgb(1, 1, 1, alpha = 0.85))


#significant snps
significant_snps <- results %>%
  filter(QVAL < 0.05)

write_csv(significant_snps, "PCAdapt_significant_SNPs.csv")




#top SNPs
top_snps <- results %>% 
  filter(QVAL <= 1e-20)

write_csv(top_snps, "PCAdapt_top_SNPs.csv")








#Mamhattan Plot
chrom_order <- c("2L", "2R", "3L", "3R", "X")
results$CHROM <- factor(results$CHROM, levels = chrom_order)

chr_info <- results %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS)) %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))

results <- results %>%
  left_join(chr_info, by = "CHROM") %>%
  mutate(BPcum = POS + chr_start)

results$minusLog10Q <- -log10(results$QVAL)


log_cutoff <- -log10(cutoff)


#genes
candidates <- read_csv("candidate_genes.csv")
head(candidates)

candidates <- candidates %>%
  mutate(CHROM = factor(CHROM, levels = chrom_order)) %>%
  arrange(CHROM, center_POS)
head(candidates)

chr_cand_info <- results %>%
  group_by(CHROM) %>%
  summarize(chr_len = max(POS)) %>%
  mutate(chr_start = cumsum(lag(chr_len, default = 0)))
candidates <- candidates %>%
  left_join(chr_cand_info, by = "CHROM") %>%
  mutate(BPcum = center_POS + chr_start)
head(candidates)

candidates <- candidates %>%
  mutate(y_end = c(30, 30, 65, 45, 37, 45, 40, 70, 40, 70, 40, 45, 37)) %>%
  mutate(y_start = y_end + 5)
head(candidates)


ggplot(results, aes(x = BPcum, y = minusLog10Q, color = CHROM)) +
  geom_point(alpha = 0.7, size = 0.6) +
  #geom_hline(yintercept = 20, linetype = "dashed", color = "red", linewidth = 0.8) +
  scale_color_manual(values = rep(c("gray20", "gray50"), length.out = length(chrom_order))) +
  scale_x_continuous(
    label = chrom_order,
    breaks = results %>% group_by(CHROM) %>% summarize(center = mean(BPcum)) %>% pull(center)
  ) +
  labs(x = "Chromosome", y = expression(-log[10](P))) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  
  geom_segment(data = candidates,
               aes(x = BPcum, xend = BPcum, y = y_start, yend = y_end),
               arrow = arrow(length = unit(0.1, "cm")),
               color = "blue", inherit.aes = FALSE) +
  
  # Gene labels
  geom_text(data = candidates,
            aes(x = BPcum, y = y_start + 1, label = Symbol),
            angle = 45, hjust = 0, vjust = 0.5,
            size = 2.5, color = "blue", fontface = "bold",
            inherit.aes = FALSE)







