library(tidyverse)
library(fst4pg)
library(gplots)
library(RColorBrewer)


#set directory
setwd("~/Desktop/dros_speciation/Fst/fst4pg")



#opening files
pop_n_maf_df <- read.csv("~/Desktop/dros_speciation/Fst/Hudson_FST/population_allele_frequencies.csv")


#data cleaning
n_df <- pop_n_maf_df %>%
  #select(CHROM, POS, ZS_N, ZH_N, ZR_N, Zam_N, SAfr_N, WCAfr_N, EAfr_N, OOA.OW_N, OOA.NW_N)
  select(CHROM, POS, Zim_N, Zam_N, SAfr_N, WCAfr_N, EAfr_N, OOA.OW_N, OOA.NW_N)

freq_df <- pop_n_maf_df %>%
  #select(CHROM, POS, ZS_maf, ZH_maf, ZR_maf, Zam_maf, SAfr_maf, WCAfr_maf, EAfr_maf, OOA.OW_maf, OOA.NW_maf)
  select(CHROM, POS, Zim_maf, Zam_maf, SAfr_maf, WCAfr_maf, EAfr_maf, OOA.OW_maf, OOA.NW_maf)

names(n_df) <- gsub('_N', '', names(n_df))
names(freq_df) <- gsub('_maf', '', names(freq_df))



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



n_A <- n_df %>%
  filter(CHROM != "X") %>%
  select(!c(CHROM, POS))


freq_A <- freq_df %>%
  filter(CHROM != "X") %>%
  select(!c(CHROM, POS))


############################################################################################################################################

  
#building database
n_freq_DB_X <- BuildFreqNbG(freq_ChrX, n_ChrX)
#estimate Hudson Fst
HFmX <- HudsonFst.m(n_freq_DB_X)
#estimate gw Fst
HFgwX <-  HudsonFst.gw(HFmX)



#building database
n_freq_DB_2L <- BuildFreqNbG(freq_Chr2L, n_Chr2L)
#estimate Hudson Fst
HFm2L <- HudsonFst.m(n_freq_DB_2L)
#estimate gw Fst
HFgw2L <-  HudsonFst.gw(HFm2L)



#building database
n_freq_DB_2R <- BuildFreqNbG(freq_Chr2R, n_Chr2R)
#estimate Hudson Fst
HFm2R <- HudsonFst.m(n_freq_DB_2R)
#estimate gw Fst
HFgw2R <-  HudsonFst.gw(HFm2R)



#building database
n_freq_DB_3L <- BuildFreqNbG(freq_Chr3L, n_Chr3L)
#estimate Hudson Fst
HFm3L <- HudsonFst.m(n_freq_DB_3L)
#estimate gw Fst
HFgw3L <-  HudsonFst.gw(HFm3L)



#building database
n_freq_DB_3R <- BuildFreqNbG(freq_Chr3R, n_Chr3R)
#estimate Hudson Fst
HFm3R <- HudsonFst.m(n_freq_DB_3R)
#estimate gw Fst
HFgw3R <-  HudsonFst.gw(HFm3R)



#building database
n_freq_DB_A <- BuildFreqNbG(freq_A, n_A)
#estimate Hudson Fst
HFmA <- HudsonFst.m(n_freq_DB_A)
#estimate gw Fst
HFgwA <-  HudsonFst.gw(HFmA)



#X and Autosome matrix
gwm <- HFgwA
lower_diag_indices <- lower.tri(HFgwA)
gwm[lower_diag_indices] <- HFgwX[lower_diag_indices]



#heatmaps
make_heatmap <- function(matrix, title){
  
  breaks <- seq(min(matrix, na.rm = TRUE), max(matrix, na.rm = TRUE), length.out = 101)
  colors <- colorRampPalette(c("white","firebrick1"))(100)
  
  heatmap.2(
    matrix,
    trace = "none",  # No trace lines
    dendrogram = "none",  # No dendrogram
    col = colors,  # Color palette
    breaks = breaks,
    Rowv = FALSE,  # No row dendrogram
    Colv = FALSE,  # No column dendrogram
    labRow = rownames(matrix),  # Row labels
    labCol = colnames(matrix),  # Column labels
    density.info = "none",  # No density info
    key = FALSE,  # Show color key
    margins = c(10, 10),  # Adjust margins
    #main = "Mean Fst",  # Main title
    cellnote = round(matrix, 3),  # Display rounded values on the heatmap
    notecol = "black",  # Color of the cell text
    notecex = 0.8,  # Size of the cell text
    colsep = 1:ncol(matrix),  # Column separators
    rowsep = 1:nrow(matrix),  # Row separators
  )
  
  
  title(main = title, line = -2)
}



make_heatmap(gwm, "Genome-Wide FST")


make_heatmap(HFgw2L, "Chr2L FST")
make_heatmap(HFgw2R, "Chr2R FST")
make_heatmap(HFgw3L, "Chr3L FST")
make_heatmap(HFgw3R, "Chr3R FST")
make_heatmap(HFgwX, "ChrX FST")





