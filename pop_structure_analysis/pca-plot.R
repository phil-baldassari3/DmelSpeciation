#packages
library(SNPRelate)
library(ggplot2)
library(dplyr)

#setting directory
setwd("/Users/philipbaldassari/Desktop/dros_speciation/Pop_Structure_Analysis")

##########################################################################

#convert vcf to gds (ONLY DO ONCE)
snpgdsVCF2GDS('biallelic_maf0.1_missing0.05_SNPs_dmel_r6_Autosomes.vcf', 'biallelic_maf0.1_missing0.05_SNPs_dmel_r6_Autosomes.gds', method=c("biallelic.only", "copy.num.of.ref"), snpfirstdim=FALSE, compress.annotation="ZIP_RA.max", compress.geno="", ref.allele=NULL, ignore.chr.prefix="chr", verbose=TRUE)

##########################################################################




# SET THESE
plot_title <- "PCA on Autosomal SNPs (minor allele frequency \u2265 10%)"
gds_file <- "biallelic_maf0.1_missing0.05_SNPs_dmel_r6_Autosomes.gds"
pca_csv <- "PCA_maf0.1_Autosomes.csv"
gds <- snpgdsOpen(gds_file)




snpgdsSummary(gds, show=TRUE)

#pca
pca <- snpgdsPCA(gds, sample.id=NULL, snp.id=NULL, autosome.only=FALSE, num.thread=5)



#%variance for axes
variance <- pca$varprop*100
rounded_percent <- round(variance, 2)
head(variance)
head(rounded_percent)

#population info
poptxt <- scan("pop.txt", what=character())
head(poptxt)

#make pca data frame
tab <- data.frame(sample.id = pca$sample.id,
                  Population = poptxt,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)

head(tab)

write.csv(tab, file = pca_csv)


###plotting and formatting

#axis labels
xlab <- paste("EV2 (", toString(rounded_percent[2]), "% Variance)", sep="")
ylab <- paste("EV1 (", toString(rounded_percent[1]), "% Variance)", sep="")

xlab
ylab



#filtering
#tab <- tab %>%
  #filter(Population == "T")


#color pallet
my_colors <- c(
  "blue", "cyan", "burlywood", "blueviolet", "indianred1", "green3", "maroon2", 
  "black", "darkgreen", "magenta", "ivory3", "khaki2", "gold4", "brown", 
  "yellow", "plum1", "chocolate", "chartreuse", "skyblue", "red", "orange")

#scatterplot
ggplot(tab, aes(EV2, EV1, colour=Population)) + 
  geom_point(shape=20, size=3, alpha=I(0.7)) + 
  xlab(xlab) + ylab(ylab) +
  ggtitle(plot_title) +
  scale_color_manual(values = my_colors) +
  theme_bw()





