#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(grid)



#fst and pi plotting function

fstpiplot <- function(plotfilename, chromosome, fst_file, pifile1, pifile2, pifile3, pifile4, pifile5) {
  ##input a plot filename, a chromosome (e.g. "2L"), and the files as strings. if you do not have all 4 pifiles, input NA for pifile3 and/or pifile4

  if (chromosome == "X") {
    substr_num <- 13
  } else {
    substr_num <- 18
  }
  
  
  
  if (chromosome == "X") {
    len <- 22422827
  } else if (chromosome == "2L") {
    len <- 23011544
  } else if (chromosome == "2R") {
    len <- 21146708
  } else if (chromosome == "3L") {
    len <- 24543557
  } else if (chromosome == "3R") {
    len <- 27905053
  }
  
  
  #opening dataframes
  setwd("/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann/10kbp_windowed_GO_ann")
  
  fst <- read.csv(fst_file)
  fst <- fst %>%
    filter(Chrom == chromosome)
  
  setwd("/Users/philipbaldassari/Desktop/zim-cos_vcfs/r6_filtered_vcfs")
  
  pi1 <- read.csv(pifile1, sep='\t')
  pi2 <- read.csv(pifile2, sep='\t')
  pi3 <- read.csv(pifile3, sep='\t')
  
  if (is.na(pifile4)) {
    pi4 <- NA
  } else {
    pi4 <- read.csv(pifile4, sep='\t')
  }
  
  if (is.na(pifile5)) {
    pi5 <- NA
  } else {
    pi5 <- read.csv(pifile5, sep='\t')
  }
  
  
  setwd("/Users/philipbaldassari/Desktop/pi_plots")
  
  
  #plotting
  if (is.na(pifile4)) {
    print('plot')
    
    pfst <- ggplot() + geom_line(data=fst, aes(x=window_start, y = Avg_Fst), linetype="solid", color = "black") + xlab(" ") + ylab("Mean Fst") + xlim(0,len) + ylim(0,0.25) + ggtitle(paste("Chrom:", chromosome, "Fst (10kbp sliding window)", substring(fst_file, 31, (nchar(fst_file)-substr_num)), sep=" ")) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
  
    
    p1 <- ggplot() + geom_line(data=pi1, aes(x=BIN_START, y = PI), linetype="solid", color = "red") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile1, 13, 15), "pi (10kbp sliding window)", sep=" ")) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p2 <- ggplot() + geom_line(data=pi2, aes(x=BIN_START, y = PI), linetype="solid", color = "blue") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile2, 13, 15), "pi (10kbp sliding window)", sep=" ")) + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    p3 <- ggplot() + geom_line(data=pi3, aes(x=BIN_START, y = PI), linetype="solid", color = "green") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile3, 13, 15), "pi (10kbp sliding window)", sep=" "))
    
    
    png(filename = plotfilename, width = 3000, height = 900)
    grid.newpage()
    grid.draw(rbind(ggplotGrob(pfst), ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), size = "last"))
    dev.off()
    
    
    
    
    
  } else if (is.na(pifile5)) {
    print('plot')
    
    pfst <- ggplot() + geom_line(data=fst, aes(x=window_start, y = Avg_Fst), linetype="solid", color = "black") + xlab(" ") + xlim(0,len) + ylab("Mean Fst") + ggtitle(paste("Chrom:", chromosome, "Fst (10kbp sliding window)", substring(fst_file, 31, (nchar(fst_file)-substr_num)), sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p1 <- ggplot() + geom_line(data=pi1, aes(x=BIN_START, y = PI), linetype="solid", color = "red") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.25) + ggtitle(paste(substring(pifile1, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p2 <- ggplot() + geom_line(data=pi2, aes(x=BIN_START, y = PI), linetype="solid", color = "blue") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile2, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p3 <- ggplot() + geom_line(data=pi3, aes(x=BIN_START, y = PI), linetype="solid", color = "green") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile3, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
   
    p4 <- ggplot() + geom_line(data=pi4, aes(x=BIN_START, y = PI), linetype="solid", color = "brown") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile4, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
     
    png(filename = plotfilename, width = 3000, height = 900)
    grid.newpage()
    grid.draw(rbind(ggplotGrob(pfst), ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), size = "last"))
    dev.off()
    
    
  } else {
    print('plot')
    
    
    pfst <- ggplot() + geom_line(data=fst, aes(x=window_start, y = Avg_Fst), linetype="solid", color = "black") + xlab(" ") + ylab("Mean Fst") + xlim(0,len) + ylim(0,0.25) + ggtitle(paste("Chrom:", chromosome, "Fst (10kbp sliding window)", substring(fst_file, 31, (nchar(fst_file)-substr_num)), sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p1 <- ggplot() + geom_line(data=pi1, aes(x=BIN_START, y = PI), linetype="solid", color = "red") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile1, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p2 <- ggplot() + geom_line(data=pi2, aes(x=BIN_START, y = PI), linetype="solid", color = "blue") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile2, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p3 <- ggplot() + geom_line(data=pi3, aes(x=BIN_START, y = PI), linetype="solid", color = "green") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile3, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    p4 <- ggplot() + geom_line(data=pi4, aes(x=BIN_START, y = PI), linetype="solid", color = "brown") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile4, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())

    p5 <- ggplot() + geom_line(data=pi5, aes(x=BIN_START, y = PI), linetype="solid", color = "purple") + xlab(" ") + ylab("pi") + xlim(0,len) + ylim(0,0.025) + ggtitle(paste(substring(pifile5, 13, 15), "pi (10kbp sliding window)", sep=" "))  + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
    
    
    png(filename = plotfilename, width = 3000, height = 900)
    grid.newpage()
    grid.draw(rbind(ggplotGrob(pfst), ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), ggplotGrob(p5), size = "last"))
    dev.off()
    
  }
  
}





#plotting

#fstpiplot("testplot.png", "X", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv", "win10000_pi_ZS_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", NA, NA)


fstpiplot("ChrX_ZS_RAL_ZI.png", "X", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv", "win10000_pi_ZS_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2L_ZS_RAL_ZI.png", "2L", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2R_ZS_RAL_ZI.png", "2R", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3L_ZS_RAL_ZI.png", "3L", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3R_ZS_RAL_ZI.png", "3R", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", NA, NA)




fstpiplot("ChrX_ZS_RAL_ZI_FR_SAfr.png", "X", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv", "win10000_pi_ZS_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_ChrX.vcf.windowed.pi")
fstpiplot("Chr2L_ZS_RAL_ZI_FR_SAfr.png", "2L", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2L.vcf.windowed.pi")
fstpiplot("Chr2R_ZS_RAL_ZI_FR_SAfr.png", "2R", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2R.vcf.windowed.pi")
fstpiplot("Chr3L_ZS_RAL_ZI_FR_SAfr.png", "3L", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3L.vcf.windowed.pi")
fstpiplot("Chr3R_ZS_RAL_ZI_FR_SAfr.png", "3R", "windowed_10kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3R.vcf.windowed.pi")




fstpiplot("ChrX_ZH_RAL_ZI.png", "X", "windowed_10kbp_GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv", "win10000_pi_ZH_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2L_ZH_RAL_ZI.png", "2L", "windowed_10kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZH_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2R_ZH_RAL_ZI.png", "2R", "windowed_10kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZH_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3L_ZH_RAL_ZI.png", "3L", "windowed_10kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZH_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3R_ZH_RAL_ZI.png", "3R", "windowed_10kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZH_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", NA, NA)




fstpiplot("ChrX_ZW_RAL_ZI.png", "X", "windowed_10kbp_GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv", "win10000_pi_ZW_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2L_ZW_RAL_ZI.png", "2L", "windowed_10kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZW_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2R_ZW_RAL_ZI.png", "2R", "windowed_10kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZW_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3L_ZW_RAL_ZI.png", "3L", "windowed_10kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZW_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3R_ZW_RAL_ZI.png", "3R", "windowed_10kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv", "win10000_pi_ZW_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", NA, NA)





fstpiplot("ChrX_ZS_ZH_ZW.png", "X", "windowed_10kbp_GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv", "win10000_pi_ZS_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZH_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZW_r6_zim-cos_ChrX.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2L_ZS_ZH_ZW.png", "2L", "windowed_10kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZH_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZW_r6_zim-cos_Chr2L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr2R_ZS_ZH_ZW.png", "2R", "windowed_10kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZH_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZW_r6_zim-cos_Chr2R.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3L_ZS_ZH_ZW.png", "3L", "windowed_10kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZH_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZW_r6_zim-cos_Chr3L.vcf.windowed.pi", NA, NA)
fstpiplot("Chr3R_ZS_ZH_ZW.png", "3R", "windowed_10kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv", "win10000_pi_ZS_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZH_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZW_r6_zim-cos_Chr3R.vcf.windowed.pi", NA, NA)






fstpiplot("ChrX_FR_vs_RAL_ZI_SAfr.png", "X", "windowed_10kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_ChrX.vcf.windowed.pi", NA)
fstpiplot("Chr2L_FR_vs_RAL_ZI_SAfr.png", "2L", "windowed_10kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2L.vcf.windowed.pi", NA)
fstpiplot("Chr2R_FR_vs_RAL_ZI_SAfr.png", "2R", "windowed_10kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2R.vcf.windowed.pi", NA)
fstpiplot("Chr3L_FR_vs_RAL_ZI_SAfr.png", "3L", "windowed_10kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3L.vcf.windowed.pi", NA)
fstpiplot("Chr3R_FR_vs_RAL_ZI_SAfr.png", "3R", "windowed_10kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3R.vcf.windowed.pi", NA)





fstpiplot("ChrX_RAL_vs_FR_ZI_SAfr.png", "X", "windowed_10kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_ChrX.vcf.windowed.pi", NA)
fstpiplot("Chr2L_RAL_vs_FR_ZI_SAfr.png", "2L", "windowed_10kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2L.vcf.windowed.pi", NA)
fstpiplot("Chr2R_RAL_vs_FR_ZI_SAfr.png", "2R", "windowed_10kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2R.vcf.windowed.pi", NA)
fstpiplot("Chr3L_RAL_vs_FR_ZI_SAfr.png", "3L", "windowed_10kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3L.vcf.windowed.pi", NA)
fstpiplot("Chr3R_RAL_vs_FR_ZI_SAfr.png", "3R", "windowed_10kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3R.vcf.windowed.pi", NA)






fstpiplot("ChrX_ZI_vs_RAL_FR_SAfr.png", "X", "windowed_10kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_ChrX.vcf.windowed.pi", NA)
fstpiplot("Chr2L_ZI_vs_RAL_FR_SAfr.png", "2L", "windowed_10kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2L.vcf.windowed.pi", NA)
fstpiplot("Chr2R_ZI_vs_RAL_FR_SAfr.png", "2R", "windowed_10kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2R.vcf.windowed.pi", NA)
fstpiplot("Chr3L_ZI_vs_RAL_FR_SAfr.png", "3L", "windowed_10kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3L.vcf.windowed.pi", NA)
fstpiplot("Chr3R_ZI_vs_RAL_FR_SAfr.png", "3R", "windowed_10kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3R.vcf.windowed.pi", NA)





fstpiplot("ChrX_SAfr_vs_RAL_FR_ZI.png", "X", "windowed_10kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv", "win10000_pi_RAL_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_ChrX.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_ChrX.vcf.windowed.pi", NA)
fstpiplot("Chr2L_SAfr_vs_RAL_FR_ZI.png", "2L", "windowed_10kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2L.vcf.windowed.pi", NA)
fstpiplot("Chr2R_SAfr_vs_RAL_FR_ZI.png", "2R", "windowed_10kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr2R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr2R.vcf.windowed.pi", NA)
fstpiplot("Chr3L_SAfr_vs_RAL_FR_ZI.png", "3L", "windowed_10kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3L.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3L.vcf.windowed.pi", NA)
fstpiplot("Chr3R_SAfr_vs_RAL_FR_ZI.png", "3R", "windowed_10kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv", "win10000_pi_RAL_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_FR_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_ZI_r6_zim-cos_Chr3R.vcf.windowed.pi", "win10000_pi_SAfr_r6_zim-cos_Chr3R.vcf.windowed.pi", NA)












