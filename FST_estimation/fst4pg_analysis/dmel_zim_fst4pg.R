library(tidyverse)
library(fst4pg)

#set directory
setwd("~/Desktop/dros_speciation/Fst/fst4pg")



#opening files
pop_n_maf_df <- read.csv("~/Desktop/dros_speciation/Fst/Hudson_FST/population_allele_frequencies.csv")
info_df <- read.csv("~/Desktop/dros_speciation/Fst/Hudson_FST/REF_and_ALT_alleles.csv")



#data cleaning
info_df$MARKER <- NA
names(info_df) <- gsub('CHROM', 'CHR', names(info_df))
names(info_df) <- gsub('REF', 'Ref', names(info_df))
names(info_df) <- gsub('ALT', 'Alt', names(info_df))

n_df <- pop_n_maf_df %>%
  select(CHROM, POS, Zim_N, ZS_N, ZH_N, ZR_N, ZC_N, WCAfr_N, EAfr_N, Zam_N, OOA.OW_N, OOA.NW_N, SAfr_N)

freq_df <- pop_n_maf_df %>%
  select(CHROM, POS, Zim_maf, ZS_maf, ZH_maf, ZR_maf, ZC_maf, WCAfr_maf, EAfr_maf, Zam_maf, OOA.OW_maf, OOA.NW_maf, SAfr_maf)

names(n_df) <- gsub('_N', '', names(n_df))
names(freq_df) <- gsub('_maf', '', names(freq_df))




#separating chromosomes
info_Chr2L <- info_df %>%
  filter(CHR == "2L")
info_Chr2R <- info_df %>%
  filter(CHR == "2R")
info_Chr3L <- info_df %>%
  filter(CHR == "3L")
info_Chr3R <- info_df %>%
  filter(CHR == "3R")
info_ChrX <- info_df %>%
  filter(CHR == "X")

n_Chr2L <- n_df %>%
  filter(CHROM == "2L") %>%
  select(!c(CHROM, POS))
n_Chr2R <- n_df %>%
  filter(CHROM == "2R") %>%
  select(!c(CHROM, POS))
n_Chr3L <- n_df %>%
  filter(CHROM == "3L") %>%
  select(!c(CHROM, POS))
n_Chr3R <- n_df %>%
  filter(CHROM == "3R") %>%
  select(!c(CHROM, POS))
n_ChrX <- n_df %>%
  filter(CHROM == "X") %>%
  select(!c(CHROM, POS))

freq_Chr2L <- freq_df %>%
  filter(CHROM == "2L") %>%
  select(!c(CHROM, POS))
freq_Chr2R <- freq_df %>%
  filter(CHROM == "2R") %>%
  select(!c(CHROM, POS))
freq_Chr3L <- freq_df %>%
  filter(CHROM == "3L") %>%
  select(!c(CHROM, POS))
freq_Chr3R <- freq_df %>%
  filter(CHROM == "3R") %>%
  select(!c(CHROM, POS))
freq_ChrX <- freq_df %>%
  filter(CHROM == "X") %>%
  select(!c(CHROM, POS))


############################################################################################################################################


#main function
run_fst4pg <- function(focalpop, frequency_df, sample_size_df, information_df, chromosome){
  
  #building database
  n_freq_DB <- BuildFreqNbG(frequency_df, sample_size_df)
  
  #estimate Hudson Fst
  HFm <- HudsonFst.m(n_freq_DB)
  
  #pick pops
  Two_Pops <- list(First=focalpop, Second=c("WCAfr", "EAfr", "Zam", "OOA.OW", "OOA.NW", "SAfr"))
  
  #Fst profile for two pops
  HFprofile <- HudsonFst.prof(HFm, Contrast=Two_Pops)
  
  #computing reference levels
  Ref_Levels <- rapply(HFprofile, median, classes = "numeric",how='list')
  print(chromosome)
  print(Ref_Levels)

  #summarize profile
  sum_profile <- ProfilingSummary(HFprofile, information_df)
  
  
  
  #saving plots
  png(file=paste0(focalpop, "_WCAfr_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_WCAfr, HFprofile$ZR_WCAfr, Ref=Ref_Levels$ZR_WCAfr, Threshold = abs(3*Ref_Levels$ZR_WCAfr)))
  dev.off()
  png(file=paste0(focalpop, "_EAfr_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_EAfr, HFprofile$ZR_EAfr, Ref=Ref_Levels$ZR_EAfr, Threshold = abs(3*Ref_Levels$ZR_EAfr)))
  dev.off()
  png(file=paste0(focalpop, "_Zam_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_Zam, HFprofile$ZR_Zam, Ref=Ref_Levels$ZR_Zam, Threshold = abs(3*Ref_Levels$ZR_Zam)))
  dev.off()
  png(file=paste0(focalpop, "_OOA.OW_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_OOA.OW, HFprofile$ZR_OOA.OW, Ref=Ref_Levels$ZR_OOA.OW, Threshold = abs(3*Ref_Levels$ZR_OOA.OW)))
  dev.off()
  png(file=paste0(focalpop, "_OOA.NW_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_OOA.NW, HFprofile$ZR_OOA.NW, Ref=Ref_Levels$ZR_OOA.NW, Threshold = abs(3*Ref_Levels$ZR_OOA.NW)))
  dev.off()
  png(file=paste0(focalpop, "_SAfr_", chromosome, "_plot.png"), width=1000, height=350)
  print(HudsonFst.plot(information_df, HFm$ZR_SAfr, HFprofile$ZR_SAfr, Ref=Ref_Levels$ZR_SAfr, Threshold = abs(3*Ref_Levels$ZR_SAfr)))
  dev.off()
  
  #saving csvs
  write.csv(sum_profile$ZR_WCAfr, paste0(focalpop, "_WCAfr_", chromosome, ".csv"), row.names = FALSE)
  write.csv(sum_profile$ZR_EAfr, paste0(focalpop, "_EAfr_", chromosome, ".csv"), row.names = FALSE)
  write.csv(sum_profile$ZR_Zam, paste0(focalpop, "_Zam_", chromosome, ".csv"), row.names = FALSE)
  write.csv(sum_profile$ZR_OOA.OW, paste0(focalpop, "_OOA.OW_", chromosome, ".csv"), row.names = FALSE)
  write.csv(sum_profile$ZR_OOA.NW, paste0(focalpop, "_OOA.NW_", chromosome, ".csv"), row.names = FALSE)
  write.csv(sum_profile$ZR_SAfr, paste0(focalpop, "_SAfr_", chromosome, ".csv"), row.names = FALSE)
  
  
  
  #filtering for regions above threshold
  top_Z_WCAfr <- sum_profile$ZR_WCAfr %>%
    filter(Fst > abs(3*Ref_Levels$ZR_WCAfr))
  
  top_Z_EAfr <- sum_profile$ZR_EAfr %>%
    filter(Fst > abs(3*Ref_Levels$ZR_EAfr))
  
  top_Z_Zam <- sum_profile$ZR_Zam %>%
    filter(Fst > abs(3*Ref_Levels$ZR_Zam))
  
  top_Z_OOA.OW <- sum_profile$ZR_OOA.OW %>%
    filter(Fst > abs(3*Ref_Levels$ZR_OOA.OW))
  
  top_Z_OOA.NW <- sum_profile$ZR_OOA.NW %>%
    filter(Fst > abs(3*Ref_Levels$ZR_OOA.NW))
  
  top_Z_SAfr <- sum_profile$ZR_SAfr %>%
    filter(Fst > abs(3*Ref_Levels$ZR_SAfr))
  
  
  #saving top region csvs
  write.csv(top_Z_WCAfr, paste0("top_regios_", focalpop, "_WCAfr_", chromosome, ".csv"), row.names = FALSE)
  write.csv(top_Z_EAfr, paste0("top_regios_", focalpop, "_EAfr_", chromosome, ".csv"), row.names = FALSE)
  write.csv(top_Z_Zam, paste0("top_regios_", focalpop, "_Zam_", chromosome, ".csv"), row.names = FALSE)
  write.csv(top_Z_OOA.OW, paste0("top_regios_", focalpop, "_OOA.OW_", chromosome, ".csv"), row.names = FALSE)
  write.csv(top_Z_OOA.NW, paste0("top_regios_", focalpop, "_OOA.NW_", chromosome, ".csv"), row.names = FALSE)
  write.csv(top_Z_SAfr, paste0("top_regios_", focalpop, "_SAfr_", chromosome, ".csv"), row.names = FALSE)
  
  
  
}



############################################################################################################################################


run_fst4pg("ZR", freq_Chr2L, n_Chr2L, info_Chr2L, "Chr2L")
run_fst4pg("ZR", freq_Chr2R, n_Chr2R, info_Chr2R, "Chr2R")
run_fst4pg("ZR", freq_Chr3L, n_Chr3L, info_Chr3L, "Chr3L")
run_fst4pg("ZR", freq_Chr3R, n_Chr3R, info_Chr3R, "Chr3R")
run_fst4pg("ZR", freq_ChrX, n_ChrX, info_ChrX, "ChrX")



















