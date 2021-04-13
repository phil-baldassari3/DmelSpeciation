#packages
library(SNPRelate)
library(ggplot2)
library(dplyr)

#setting directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_pca-plots")

##################################################################################
#convert vcf to gds (ONLY DO THIS ONCE!)
#full genome
snpgdsVCF2GDS('zim-cos.vcf', 'zim-cos.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
#ChrX
snpgdsVCF2GDS('zim-cos_ChrX.vcf', 'zim-cos_ChrX.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
#Chr2L
snpgdsVCF2GDS('zim-cos_Chr2L.vcf', 'zim-cos_Chr2L.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
#Chr2R
snpgdsVCF2GDS('zim-cos_Chr2R.vcf', 'zim-cos_Chr2R.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
#Chr3L
snpgdsVCF2GDS('zim-cos_Chr3L.vcf', 'zim-cos_Chr3L.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
#Chr3R
snpgdsVCF2GDS('zim-cos_Chr3R.vcf', 'zim-cos_Chr3R.gds', method=("copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)
##################################################################################

#Full genome
gds <- snpgdsOpen('zim-cos.gds')
snpgdsSummary(gds, show=TRUE)
#ChrX
gdsX <- snpgdsOpen('zim-cos_ChrX.gds')
snpgdsSummary(gdsX, show=TRUE)
#Chr2L
gds2L <- snpgdsOpen('zim-cos_Chr2L.gds')
snpgdsSummary(gds2L, show=TRUE)
#Chr2R
gds2R <- snpgdsOpen('zim-cos_Chr2R.gds')
snpgdsSummary(gds2R, show=TRUE)
#Chr3L
gds3L <- snpgdsOpen('zim-cos_Chr3L.gds')
snpgdsSummary(gds3L, show=TRUE)
#Chr3R
gds3R <- snpgdsOpen('zim-cos_Chr3R.gds')
snpgdsSummary(gds3R, show=TRUE)


#pca (ONLY DO ONCE?)
#Full Genome
pca <- snpgdsPCA(gds, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)
#ChrX
pcaX <- snpgdsPCA(gdsX, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)
#Chr2L
pca2L <- snpgdsPCA(gds2L, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)
#Chr2R
pca2R <- snpgdsPCA(gds2R, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)
#Chr3L
pca3L <- snpgdsPCA(gds3L, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)
#Chr3R
pca3R <- snpgdsPCA(gds3R, sample.id=NULL, snp.id=NULL, autosome.only=FALSE)


#%variance for axes
#Full Genome
var <- pca$varprop*100
rounded <- round(var, 2)
head(var)
head(rounded)
#ChrX
varX <- pcaX$varprop*100
roundedX <- round(varX, 2)
head(varX)
head(roundedX)
#Chr2L
var2L <- pca2L$varprop*100
rounded2L <- round(var2L, 2)
head(var2L)
head(rounded2L)
#Chr2R
var2R <- pca2R$varprop*100
rounded2R <- round(var2R, 2)
head(var2R)
head(rounded2R)
#Chr3L
var3L <- pca3L$varprop*100
rounded3L <- round(var3L, 2)
head(var3L)
head(rounded3L)
#Chr3R
var3R <- pca3R$varprop*100
rounded3R <- round(var3R, 2)
head(var3R)
head(rounded3R)


#population info
poptxt <- scan("pop.txt", what=character())
head(poptxt)


#make pca data frames
#Full Genome
tab <- data.frame(sample.id = pca$sample.id,
                  Population = poptxt,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab)
#ChrX
tabX <- data.frame(sample.id = pcaX$sample.id,
                  Population = poptxt,
                  EV1 = pcaX$eigenvect[,1],    # the first eigenvector
                  EV2 = pcaX$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tabX)
#Chr2L
tab2L <- data.frame(sample.id = pca2L$sample.id,
                  Population = poptxt,
                  EV1 = pca2L$eigenvect[,1],    # the first eigenvector
                  EV2 = pca2L$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab2L)
#Chr2R
tab2R <- data.frame(sample.id = pca2R$sample.id,
                  Population = poptxt,
                  EV1 = pca2R$eigenvect[,1],    # the first eigenvector
                  EV2 = pca2R$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab2R)
#Chr3L
tab3L <- data.frame(sample.id = pca3L$sample.id,
                  Population = poptxt,
                  EV1 = pca3L$eigenvect[,1],    # the first eigenvector
                  EV2 = pca3L$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab3L)
#Chr3R
tab3R <- data.frame(sample.id = pca3R$sample.id,
                  Population = poptxt,
                  EV1 = pca3R$eigenvect[,1],    # the first eigenvector
                  EV2 = pca3R$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab3R)



###plotting and formatting

#axis labels
#Full Genome
xlab <- paste("EV2 (", toString(rounded[2]), "% Variance)", sep="")
ylab <- paste("EV1 (", toString(rounded[1]), "% Variance)", sep="")
xlab
ylab
#ChrX
xlabX <- paste("EV2 (", toString(roundedX[2]), "% Variance)", sep="")
ylabX <- paste("EV1 (", toString(roundedX[1]), "% Variance)", sep="")
xlabX
ylabX
#Chr2L
xlab2L <- paste("EV2 (", toString(rounded2L[2]), "% Variance)", sep="")
ylab2L <- paste("EV1 (", toString(rounded2L[1]), "% Variance)", sep="")
xlab2L
ylab2L
#Chr2R
xlab2R <- paste("EV2 (", toString(rounded2R[2]), "% Variance)", sep="")
ylab2R <- paste("EV1 (", toString(rounded2R[1]), "% Variance)", sep="")
xlab2R
ylab2R
#Chr3L
xlab3L <- paste("EV2 (", toString(rounded3L[2]), "% Variance)", sep="")
ylab3L <- paste("EV1 (", toString(rounded3L[1]), "% Variance)", sep="")
xlab3L
ylab3L
#Chr3R
xlab3R <- paste("EV2 (", toString(rounded3R[2]), "% Variance)", sep="")
ylab3R <- paste("EV1 (", toString(rounded3R[1]), "% Variance)", sep="")
xlab3R
ylab3R



#scatterplot
#Full Genome
ggplot(tab, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab) +ylab(ylab) +ggtitle("Whole Genome (minor alleles filtered )")
#ChrX
ggplot(tabX, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlabX) +ylab(ylabX) +ggtitle("ChrX (minor alleles filtered )")
#Chr2L
ggplot(tab2L, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab2L) +ylab(ylab2L) +ggtitle("Chr2L (minor alleles filtered )")
#Chr2R
ggplot(tab2R, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab2R) +ylab(ylab2R) +ggtitle("Chr2R (minor alleles filtered )")
#Chr3L
ggplot(tab3L, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab3L) +ylab(ylab3L) +ggtitle("Chr3L (minor alleles filtered )")
#Chr3R
ggplot(tab3R, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab3R) +ylab(ylab3R) +ggtitle("Chr3R (minor alleles filtered )")




#if you need to filter to only see specific populations
#Full Genome
just_african <- tab %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_african, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab) +ylab(ylab)
#ChrX
just_africanX <- tabX %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_africanX, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlabX) +ylab(ylabX)
#Chr2L
just_african2L <- tab2L %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_african2L, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab2L) +ylab(ylab2L)
#Chr2R
just_african2R <- tab2R %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_african2R, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab2R) +ylab(ylab2R)
#Chr3L
just_african3L <- tab3L %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_african3L, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab3L) +ylab(ylab3L)
#Chr3R
just_african3R <- tab3R %>%
  filter(!(Population == "FR"), !(Population == "RAL"))

ggplot(just_african3R, aes(EV2, EV1, colour=Population)) +geom_point(shape=19, size=2) +xlab(xlab3R) +ylab(ylab3R)








