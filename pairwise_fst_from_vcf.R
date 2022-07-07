#Fst using the Hudson estimator from Hudson et al. 1992 and Bhatia et al. 2013
#Fst_hud = (((p1 - p2)**2) - ((p1*(1-p1))/(n1-1)) - ((p2*(1-p2))/(n2-1))) / ((p1*(1-p2)) + (p2*(1-p1)))

#loading packages
library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)

#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_vcfs/allele-counts_per-pop_downsampled")

#reading files
ChrX <- read.csv("vcf_ChrX_allele-counts_downsample.csv")
head(ChrX)
Chr2L <- read.csv("vcf_Chr2L_allele-counts_downsample.csv")
Chr2R <- read.csv("vcf_Chr2R_allele-counts_downsample.csv")
Chr3L <- read.csv("vcf_Chr3L_allele-counts_downsample.csv")
Chr3R <- read.csv("vcf_Chr3R_allele-counts_downsample.csv")


#pairwise Fst
##ZS vs RAL
ChrX_ZSvRAL <- ChrX %>%
  select(POS, n_ZS, n_RAL, maf_ZS, maf_RAL) %>%
  filter(n_ZS == 4) %>%
  filter(n_RAL == 125) %>%
  mutate(Fst_hud = ((((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS)))))

ChrX_ZSvRAL[is.na(ChrX_ZSvRAL)] = 0
mean_ChrX_ZSvRAL <- mean(ChrX_ZSvRAL$Fst_hud)
sd_ChrX_ZSvRAL <- sd(ChrX_ZSvRAL$Fst_hud)

ChrX_ZSvRAL <- ChrX_ZSvRAL %>%
  mutate(z = (Fst_hud - mean_ChrX_ZSvRAL)/sd_ChrX_ZSvRAL) %>%
  filter(z >= 3)

head(ChrX_ZSvRAL)

Chr2L_ZSvRAL <- Chr2L %>%
  select(POS, n_ZS, n_RAL, maf_ZS, maf_RAL) %>%
  filter(n_ZS == 4) %>%
  filter(n_RAL == 125) %>%
  mutate(Fst_hud = (((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))

Chr2L_ZSvRAL[is.na(Chr2L_ZSvRAL)] = 0
mean_Chr2L_ZSvRAL <- mean(Chr2L_ZSvRAL$Fst_hud)
sd_Chr2L_ZSvRAL <- sd(Chr2L_ZSvRAL$Fst_hud)

Chr2L_ZSvRAL <- Chr2L_ZSvRAL %>%
  mutate(z = (Fst_hud - mean_Chr2L_ZSvRAL)/sd_Chr2L_ZSvRAL) %>%
  filter(z >= 3)


Chr2R_ZSvRAL <- Chr2R %>%
  select(POS, n_ZS, n_RAL, maf_ZS, maf_RAL) %>%
  filter(n_ZS == 4) %>%
  filter(n_RAL == 125) %>%
  mutate(Fst_hud = (((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))

Chr2R_ZSvRAL[is.na(Chr2R_ZSvRAL)] = 0
mean_Chr2R_ZSvRAL <- mean(Chr2R_ZSvRAL$Fst_hud)
sd_Chr2R_ZSvRAL <- sd(Chr2R_ZSvRAL$Fst_hud)

Chr2R_ZSvRAL <- Chr2R_ZSvRAL %>%
  mutate(z = (Fst_hud - mean_Chr2R_ZSvRAL)/sd_Chr2R_ZSvRAL) %>%
  filter(z >= 3)


Chr3L_ZSvRAL <- Chr3L %>%
  select(POS, n_ZS, n_RAL, maf_ZS, maf_RAL) %>%
  filter(n_ZS == 4) %>%
  filter(n_RAL == 125) %>%
  mutate(Fst_hud = (((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))

Chr3L_ZSvRAL[is.na(Chr3L_ZSvRAL)] = 0
mean_Chr3L_ZSvRAL <- mean(Chr3L_ZSvRAL$Fst_hud)
sd_Chr3L_ZSvRAL <- sd(Chr3L_ZSvRAL$Fst_hud)

Chr3L_ZSvRAL <- Chr3L_ZSvRAL %>%
  mutate(z = (Fst_hud - mean_Chr3L_ZSvRAL)/sd_Chr3L_ZSvRAL) %>%
  filter(z >= 3)



Chr3R_ZSvRAL <- Chr3R %>%
  select(POS, n_ZS, n_RAL, maf_ZS, maf_RAL) %>%
  filter(n_ZS == 4) %>%
  filter(n_RAL == 125) %>%
  mutate(Fst_hud = (((maf_ZS - maf_RAL)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_RAL*(1-maf_RAL))/(n_RAL-1))) / ((maf_ZS*(1-maf_RAL)) + (maf_RAL*(1-maf_ZS))))

Chr3R_ZSvRAL[is.na(Chr3R_ZSvRAL)] = 0
mean_Chr3R_ZSvRAL <- mean(Chr3R_ZSvRAL$Fst_hud)
sd_Chr3R_ZSvRAL <- sd(Chr3R_ZSvRAL$Fst_hud)

Chr3R_ZSvRAL <- ChrX_ZSvRAL %>%
  mutate(z = (Fst_hud - mean_Chr3R_ZSvRAL)/sd_Chr3R_ZSvRAL) %>%
  filter(z >= 3)



##ZS vs ZI
ChrX_ZSvZI <- ChrX %>%
  select(POS, n_ZS, n_ZI, maf_ZS, maf_ZI) %>%
  filter(n_ZS == 4) %>%
  filter(n_ZI == 180) %>%
  mutate(Fst_hud = (((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))

ChrX_ZSvZI[is.na(ChrX_ZSvZI)] = 0
mean_ChrX_ZSvZI <- mean(ChrX_ZSvZI$Fst_hud)
sd_ChrX_ZSvZI <- sd(ChrX_ZSvZI$Fst_hud)

ChrX_ZSvZI <- ChrX_ZSvZI %>%
  mutate(z = (Fst_hud - mean_ChrX_ZSvZI)/sd_ChrX_ZSvZI) %>%
  filter(z >= 3)


Chr2L_ZSvZI <- Chr2L %>%
  select(POS, n_ZS, n_ZI, maf_ZS, maf_ZI) %>%
  filter(n_ZS == 4) %>%
  filter(n_ZI == 180) %>%
  mutate(Fst_hud = (((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))

Chr2L_ZSvZI[is.na(Chr2L_ZSvZI)] = 0
mean_Chr2L_ZSvZI <- mean(Chr2L_ZSvZI$Fst_hud)
sd_Chr2L_ZSvZI <- sd(Chr2L_ZSvZI$Fst_hud)

Chr2L_ZSvZI <- Chr2L_ZSvZI %>%
  mutate(z = (Fst_hud - mean_Chr2L_ZSvZI)/sd_Chr2L_ZSvZI) %>%
  filter(z >= 3)


Chr2R_ZSvZI <- Chr2R %>%
  select(POS, n_ZS, n_ZI, maf_ZS, maf_ZI) %>%
  filter(n_ZS == 4) %>%
  filter(n_ZI == 180) %>%
  mutate(Fst_hud = (((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))

Chr2R_ZSvZI[is.na(Chr2R_ZSvZI)] = 0
mean_Chr2R_ZSvZI <- mean(Chr2R_ZSvZI$Fst_hud)
sd_Chr2R_ZSvZI <- sd(Chr2R_ZSvZI$Fst_hud)

Chr2R_ZSvZI <- Chr2R_ZSvZI %>%
  mutate(z = (Fst_hud - mean_Chr2R_ZSvZI)/sd_Chr2R_ZSvZI) %>%
  filter(z >= 3)


Chr3L_ZSvZI <- Chr3L %>%
  select(POS, n_ZS, n_ZI, maf_ZS, maf_ZI) %>%
  filter(n_ZS == 4) %>%
  filter(n_ZI == 180) %>%
  mutate(Fst_hud = (((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))

Chr3L_ZSvZI[is.na(Chr3L_ZSvZI)] = 0
mean_Chr3L_ZSvZI <- mean(Chr3L_ZSvZI$Fst_hud)
sd_Chr3L_ZSvZI <- sd(Chr3L_ZSvZI$Fst_hud)

Chr3L_ZSvZI <- Chr3L_ZSvZI %>%
  mutate(z = (Fst_hud - mean_Chr3L_ZSvZI)/sd_Chr3L_ZSvZI) %>%
  filter(z >= 3)


Chr3R_ZSvZI <- Chr3R %>%
  select(POS, n_ZS, n_ZI, maf_ZS, maf_ZI) %>%
  filter(n_ZS == 4) %>%
  filter(n_ZI == 180) %>%
  mutate(Fst_hud = (((maf_ZS - maf_ZI)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_ZI*(1-maf_ZI))/(n_ZI-1))) / ((maf_ZS*(1-maf_ZI)) + (maf_ZI*(1-maf_ZS))))

Chr3R_ZSvZI[is.na(Chr3R_ZSvZI)] = 0
mean_Chr3R_ZSvZI <- mean(Chr3R_ZSvZI$Fst_hud)
sd_Chr3R_ZSvZI <- sd(Chr3R_ZSvZI$Fst_hud)

Chr3R_ZSvZI <- Chr3R_ZSvZI %>%
  mutate(z = (Fst_hud - mean_Chr3R_ZSvZI)/sd_Chr3R_ZSvZI) %>%
  filter(z >= 3)



##ZS vs SAfr
ChrX_ZSvSAfr <- ChrX %>%
  select(POS, n_ZS, n_SAfr, maf_ZS, maf_SAfr) %>%
  filter(n_ZS == 4) %>%
  filter(n_SAfr == 40) %>%
  mutate(Fst_hud = (((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))

ChrX_ZSvSAfr[is.na(ChrX_ZSvSAfr)] = 0
mean_ChrX_ZSvSAfr <- mean(ChrX_ZSvSAfr$Fst_hud)
sd_ChrX_ZSvSAfr <- sd(ChrX_ZSvSAfr$Fst_hud)

ChrX_ZSvSAfr <- ChrX_ZSvSAfr %>%
  mutate(z = (Fst_hud - mean_ChrX_ZSvSAfr)/sd_ChrX_ZSvSAfr) %>%
  filter(z >= 3)


Chr2L_ZSvSAfr <- Chr2L %>%
  select(POS, n_ZS, n_SAfr, maf_ZS, maf_SAfr) %>%
  filter(n_ZS == 4) %>%
  filter(n_SAfr == 40) %>%
  mutate(Fst_hud = (((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))

Chr2L_ZSvSAfr[is.na(Chr2L_ZSvSAfr)] = 0
mean_Chr2L_ZSvSAfr <- mean(Chr2L_ZSvSAfr$Fst_hud)
sd_Chr2L_ZSvSAfr <- sd(Chr2L_ZSvSAfr$Fst_hud)

Chr2L_ZSvSAfr <- Chr2L_ZSvSAfr %>%
  mutate(z = (Fst_hud - mean_Chr2L_ZSvSAfr)/sd_Chr2L_ZSvSAfr) %>%
  filter(z >= 3)


Chr2R_ZSvSAfr <- Chr2R %>%
  select(POS, n_ZS, n_SAfr, maf_ZS, maf_SAfr) %>%
  filter(n_ZS == 4) %>%
  filter(n_SAfr == 40) %>%
  mutate(Fst_hud = (((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))

Chr2R_ZSvSAfr[is.na(Chr2R_ZSvSAfr)] = 0
mean_Chr2R_ZSvSAfr <- mean(Chr2R_ZSvSAfr$Fst_hud)
sd_Chr2R_ZSvSAfr <- sd(Chr2R_ZSvSAfr$Fst_hud)

Chr2R_ZSvZI <- Chr2R_ZSvZI %>%
  mutate(z = (Fst_hud - mean_Chr2R_ZSvSAfr)/sd_Chr2R_ZSvSAfr) %>%
  filter(z >= 3)


Chr3L_ZSvSAfr <- Chr3L %>%
  select(POS, n_ZS, n_SAfr, maf_ZS, maf_SAfr) %>%
  filter(n_ZS == 4) %>%
  filter(n_SAfr == 40) %>%
  mutate(Fst_hud = (((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))

Chr3L_ZSvSAfr[is.na(Chr3L_ZSvSAfr)] = 0
mean_Chr3L_ZSvSAfr <- mean(Chr3L_ZSvSAfr$Fst_hud)
sd_Chr3L_ZSvSAfr <- sd(Chr3L_ZSvSAfr$Fst_hud)

Chr3L_ZSvSAfr <- Chr3L_ZSvSAfr %>%
  mutate(z = (Fst_hud - mean_Chr3L_ZSvSAfr)/sd_Chr3L_ZSvSAfr) %>%
  filter(z >= 3)


Chr3R_ZSvSAfr <- Chr3R %>%
  select(POS, n_ZS, n_SAfr, maf_ZS, maf_SAfr) %>%
  filter(n_ZS == 4) %>%
  filter(n_SAfr == 40) %>%
  mutate(Fst_hud = (((maf_ZS - maf_SAfr)**2) - ((maf_ZS*(1-maf_ZS))/(n_ZS-1)) - ((maf_SAfr*(1-maf_SAfr))/(n_SAfr-1))) / ((maf_ZS*(1-maf_SAfr)) + (maf_SAfr*(1-maf_ZS))))

Chr3R_ZSvSAfr[is.na(Chr3R_ZSvSAfr)] = 0
mean_Chr3R_ZSvSAfr <- mean(Chr3R_ZSvSAfr$Fst_hud)
sd_Chr3R_ZSvSAfr <- sd(Chr3R_ZSvSAfr$Fst_hud)

Chr3R_ZSvSAfr <- Chr3R_ZSvSAfr %>%
  mutate(z = (Fst_hud - mean_Chr3R_ZSvSAfr)/sd_Chr3R_ZSvSAfr) %>%
  filter(z >= 3)



#plotting
##ZS vs. RAL
png(filename="fst_plots/ChrX_fst_hud_ZSvsRAL.png", width = 1000, height = 300)
ggplot(ChrX_ZSvRAL, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("ChrX Fst (Hudson 1992) ZS (n=4) vs. RAL (n=125) z>=3")
dev.off()

png(filename="fst_plots/Chr2L_fst_hud_ZSvsRAL.png", width = 1000, height = 300)
ggplot(Chr2L_ZSvRAL, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))   +xlab("Position") +ylab("Fst") +ggtitle("Chr2L Fst (Hudson 1992) ZS (n=4) vs. RAL (n=125) z>=3")
dev.off()

png(filename="fst_plots/Chr2R_fst_hud_ZSvsRAL.png", width = 1000, height = 300)
ggplot(Chr2R_ZSvRAL, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr2R Fst (Hudson 1992) ZS (n=4) vs. RAL (n=125) z>=3")
dev.off()

png(filename="fst_plots/Chr3L_fst_hud_ZSvsRAL.png", width = 1000, height = 300)
ggplot(Chr3L_ZSvRAL, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3L Fst (Hudson 1992) ZS (n=4) vs. RAL (n=125) z>=3")
dev.off()

png(filename="fst_plots/Chr3R_fst_hud_ZSvsRAL.png", width = 1000, height = 300)
ggplot(Chr3R_ZSvRAL, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3R Fst (Hudson 1992) ZS (n=4) vs. RAL (n=125) z>=3")
dev.off()





##ZS vs. ZI
png(filename="fst_plots/ChrX_fst_hud_ZSvsZI.png", width = 1000, height = 300)
ggplot(ChrX_ZSvZI, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("ChrX Fst (Hudson 1992) ZS (n=4) vs. ZI (n=180) z>=3")
dev.off()

png(filename="fst_plots/Chr2L_fst_hud_ZSvsZI.png", width = 1000, height = 300)
ggplot(Chr2L_ZSvZI, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr2L Fst (Hudson 1992) ZS (n=4) vs. ZI (n=180) z>=3")
dev.off()

png(filename="fst_plots/Chr2R_fst_hud_ZSvsZI.png", width = 1000, height = 300)
ggplot(Chr2R_ZSvZI, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr2R Fst (Hudson 1992) ZS (n=4) vs. ZI (n=180) z>=3")
dev.off()

png(filename="fst_plots/Chr3L_fst_hud_ZSvsZI.png", width = 1000, height = 300)
ggplot(Chr3L_ZSvZI, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3L Fst (Hudson 1992) ZS (n=4) vs. ZI (n=180) z>=3")
dev.off()

png(filename="fst_plots/Chr3R_fst_hud_ZSvsZI.png", width = 1000, height = 300)
ggplot(Chr3R_ZSvZI, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3R Fst (Hudson 1992) ZS vs. ZI (n=180) z>=3")
dev.off()



##ZS vs SAfr
png(filename="fst_plots/ChrX_fst_hud_ZSvsSAfr.png", width = 1000, height = 300)
ggplot(ChrX_ZSvSAfr, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("ChrX Fst (Hudson 1992) ZS (n=4) vs. SAfr (n=40) z>=3")
dev.off()

png(filename="fst_plots/Chr2L_fst_hud_ZSvsSAfr.png", width = 1000, height = 300)
ggplot(Chr2L_ZSvSAfr, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr2L Fst (Hudson 1992) ZS (n=4) vs. SAfr (n=40) z>=3")
dev.off()

png(filename="fst_plots/Chr2R_fst_hud_ZSvsSAfr.png", width = 1000, height = 300)
ggplot(Chr2R_ZSvSAfr, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr2R Fst (Hudson 1992) ZS (n=4) vs. SAfr (n=40) z>=3")
dev.off()

png(filename="fst_plots/Chr3L_fst_hud_ZSvsSAfr.png", width = 1000, height = 300)
ggplot(Chr3L_ZSvSAfr, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3L Fst (Hudson 1992) ZS (n=4) vs. SAfr (n=40) z>=3")
dev.off()

png(filename="fst_plots/Chr3R_fst_hud_ZSvsSAfr.png", width = 1000, height = 300)
ggplot(Chr3R_ZSvSAfr, aes(x=POS, y=Fst_hud,)) +geom_point(shape=19, size=0.5)  +scale_x_continuous(labels = comma, limits = c(0,NA))  +xlab("Position") +ylab("Fst") +ggtitle("Chr3R Fst (Hudson 1992) ZS vs. SAfr (n=40) z>=3")
dev.off()



desc_Chr2R_ZSvZI <- Chr2R_ZSvZI %>%
  arrange(desc(Fst_hud))

desc_Chr2R_ZSvZI



