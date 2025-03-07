library(ggplot2)
library(tidyverse)
library(RColorBrewer)

#set working directory
setwd("~/Desktop/dros_speciation/Pop_Structure_Analysis/dmel_ADMIXTURE/")

#opening CV error data
CV_chrX <- read.csv("ChrX_CV.txt", sep="\t")
CV_Autosomes <- read.csv("Autosomes_CV.txt", sep="\t")

#plotting CV error data
ggplot(data=CV_chrX, aes(x=k, y=CV_error)) +
  geom_line() + geom_point() +
  ggtitle("ADMIXTURE Cross-Validation Error on ChromX") +
  xlab("K") + ylab("Cross-Validation Error") +
  scale_x_continuous(breaks = seq(min(CV_chrX$k), max(CV_chrX$k), by = 1))

ggplot(data=CV_Autosomes, aes(x=k, y=CV_error)) +
  geom_line() + geom_point() +
  ggtitle("ADMIXTURE Cross-Validation Error on Autosomes") +
  xlab("K") + ylab("Cross-Validation Error") +
  scale_x_continuous(breaks = seq(min(CV_Autosomes$k), max(CV_Autosomes$k), by = 1))




#population data
samples_table <-  read.table("Dmel_samples_plink.fam")
samples <- samples_table[,1]
pops <- scan("pop.txt", what=character())

pop_order <- c('ZW','ZS','ZH','MC','LZV','DRM','SP','SD','ZI','RG','GH','EF','EA','EG','FR','N','B','T','I','RAL','w')

#opening Q matrices
QchrX5 <- read.table("thinned100_biallelic_maf0.01_missing0.05_SNPs_dmel_r6_ChrX.5.Q")
Qautosomes5 <- read.table("thinned100_biallelic_maf0.01_missing0.05_SNPs_dmel_r6_Autosomes.5.Q")
QchrX7 <- read.table("thinned100_biallelic_maf0.01_missing0.05_SNPs_dmel_r6_ChrX.7.Q")
Qautosomes7 <- read.table("thinned100_biallelic_maf0.01_missing0.05_SNPs_dmel_r6_Autosomes.7.Q")

#editing Q matrix dataframes
QchrX5_edit <- QchrX5 %>%
  mutate(Sample = samples) %>%
  mutate(Pop = pops) %>%
  relocate(Sample, .before=1) %>%
  relocate(Pop, .before=2) %>%
  arrange(match(Pop, pop_order))

Qautosomes5_edit <- Qautosomes5 %>%
  mutate(Sample = samples) %>%
  mutate(Pop = pops) %>%
  relocate(Sample, .before=1) %>%
  relocate(Pop, .before=2) %>%
  arrange(match(Pop, pop_order))

QchrX7_edit <- QchrX7 %>%
  mutate(Sample = samples) %>%
  mutate(Pop = pops) %>%
  relocate(Sample, .before=1) %>%
  relocate(Pop, .before=2) %>%
  arrange(match(Pop, pop_order))

Qautosomes7_edit <- Qautosomes7 %>%
  mutate(Sample = samples) %>%
  mutate(Pop = pops) %>%
  relocate(Sample, .before=1) %>%
  relocate(Pop, .before=2) %>%
  arrange(match(Pop, pop_order))


#population groups for plotting
unique_Pop <- unique(QchrX7_edit$Pop)
breaks <- match(unique_Pop, Qautosomes5_edit$Pop) 
breaks <- c((breaks - 1), 602)
averages <- rollapply(breaks, 2, mean, align = "right")



#Plotting
custom_palette <- brewer.pal(7, "Set1")

barplot(t(as.matrix(QchrX5_edit[, c("V1", "V2", "V3", "V4", "V5")])), col=custom_palette,
        ylab="Ancestry", border=NA, space = 0, main="ChrX, k=5")
segments(x0 = min(breaks), x1 = max(breaks), y0 = -0.1, xpd = TRUE)
segments(x0 = breaks, y0 = -0.1, y1 = -0.05, xpd = TRUE)
mtext(unique(QchrX5_edit$Pop), side = 1, padj = 2, at = averages, cex=0.7)


barplot(t(as.matrix(Qautosomes5_edit[, c("V1", "V2", "V3", "V4", "V5")])), col=custom_palette,
        ylab="Ancestry", border=NA, space = 0, main="Autosomes, k=5")
segments(x0 = min(breaks), x1 = max(breaks), y0 = -0.1, xpd = TRUE)
segments(x0 = breaks, y0 = -0.1, y1 = -0.05, xpd = TRUE)
mtext(unique(Qautosomes5_edit$Pop), side = 1, padj = 2, at = averages, cex=0.7)


barplot(t(as.matrix(QchrX7_edit[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7")])), col=custom_palette,
        ylab="Ancestry", border=NA, space = 0, main="ChrX, k=7")
segments(x0 = min(breaks), x1 = max(breaks), y0 = -0.1, xpd = TRUE)
segments(x0 = breaks, y0 = -0.1, y1 = -0.05, xpd = TRUE)
mtext(unique(QchrX7_edit$Pop), side = 1, padj = 2, at = averages, cex=0.7)


barplot(t(as.matrix(Qautosomes7_edit[, c("V1", "V2", "V3", "V4", "V5", "V6", "V7")])), col=custom_palette,
        ylab="Ancestry", border=NA, space = 0, main="Autosomes, k=7")
segments(x0 = min(breaks), x1 = max(breaks), y0 = -0.1, xpd = TRUE)
segments(x0 = breaks, y0 = -0.1, y1 = -0.05, xpd = TRUE)
mtext(unique(Qautosomes7_edit$Pop), side = 1, padj = 2, at = averages, cex=0.7)


