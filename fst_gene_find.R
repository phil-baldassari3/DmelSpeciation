#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)


#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_downsampled")

#reading files
Fst_ZS_RAL_ZI_ChrX <- read.csv("Fst_ZS_RAL_ZI/ChrX_ZS_RAL_ZI_Fst.csv")
Fst_ZS_RAL_ZI_Chr2L <- read.csv("Fst_ZS_RAL_ZI/Chr2L_ZS_RAL_ZI_Fst.csv")
Fst_ZS_RAL_ZI_Chr2R <- read.csv("Fst_ZS_RAL_ZI/Chr2R_ZS_RAL_ZI_Fst.csv")
Fst_ZS_RAL_ZI_Chr3L <- read.csv("Fst_ZS_RAL_ZI/Chr3L_ZS_RAL_ZI_Fst.csv")
Fst_ZS_RAL_ZI_Chr3R <- read.csv("Fst_ZS_RAL_ZI/Chr3R_ZS_RAL_ZI_Fst.csv")

Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- read.csv("Fst_ZS_RAL_ZI_FR_SAfr/ChrX_ZS_RAL_ZI_FR_SAfr_Fst.csv")
Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- read.csv("Fst_ZS_RAL_ZI_FR_SAfr/Chr2L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- read.csv("Fst_ZS_RAL_ZI_FR_SAfr/Chr2R_ZS_RAL_ZI_FR_SAfr_Fst.csv")
Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- read.csv("Fst_ZS_RAL_ZI_FR_SAfr/Chr3L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- read.csv("Fst_ZS_RAL_ZI_FR_SAfr/Chr3R_ZS_RAL_ZI_FR_SAfr_Fst.csv")

Fst_Zim_RAL_ZI_ChrX <- read.csv("Fst_Zim_RAL_ZI/ChrX_Zim_RAL_ZI_Fst.csv")
Fst_Zim_RAL_ZI_Chr2L <- read.csv("Fst_Zim_RAL_ZI/Chr2L_Zim_RAL_ZI_Fst.csv")
Fst_Zim_RAL_ZI_Chr2R <- read.csv("Fst_Zim_RAL_ZI/Chr2R_Zim_RAL_ZI_Fst.csv")
Fst_Zim_RAL_ZI_Chr3L <- read.csv("Fst_Zim_RAL_ZI/Chr3L_Zim_RAL_ZI_Fst.csv")
Fst_Zim_RAL_ZI_Chr3R <- read.csv("Fst_Zim_RAL_ZI/Chr3R_Zim_RAL_ZI_Fst.csv")

Fst_ZH_RAL_ZI_ChrX <- read.csv("Fst_ZH_RAL_ZI/ChrX_ZH_RAL_ZI_Fst.csv")
Fst_ZH_RAL_ZI_Chr2L <- read.csv("Fst_ZH_RAL_ZI/Chr2L_ZH_RAL_ZI_Fst.csv")
Fst_ZH_RAL_ZI_Chr2R <- read.csv("Fst_ZH_RAL_ZI/Chr2R_ZH_RAL_ZI_Fst.csv")
Fst_ZH_RAL_ZI_Chr3L <- read.csv("Fst_ZH_RAL_ZI/Chr3L_ZH_RAL_ZI_Fst.csv")
Fst_ZH_RAL_ZI_Chr3R <- read.csv("Fst_ZH_RAL_ZI/Chr3R_ZH_RAL_ZI_Fst.csv")

Fst_ZW_RAL_ZI_ChrX <- read.csv("Fst_ZW_RAL_ZI/ChrX_ZW_RAL_ZI_Fst.csv")
Fst_ZW_RAL_ZI_Chr2L <- read.csv("Fst_ZW_RAL_ZI/Chr2L_ZW_RAL_ZI_Fst.csv")
Fst_ZW_RAL_ZI_Chr2R <- read.csv("Fst_ZW_RAL_ZI/Chr2R_ZW_RAL_ZI_Fst.csv")
Fst_ZW_RAL_ZI_Chr3L <- read.csv("Fst_ZW_RAL_ZI/Chr3L_ZW_RAL_ZI_Fst.csv")
Fst_ZW_RAL_ZI_Chr3R <- read.csv("Fst_ZW_RAL_ZI/Chr3R_ZW_RAL_ZI_Fst.csv")

Fst_ZS_ZH_ZW_ChrX <- read.csv("Fst_ZS_ZH_ZW/ChrX_ZS_ZH_ZW_Fst.csv")
Fst_ZS_ZH_ZW_Chr2L <- read.csv("Fst_ZS_ZH_ZW/Chr2L_ZS_ZH_ZW_Fst.csv")
Fst_ZS_ZH_ZW_Chr2R <- read.csv("Fst_ZS_ZH_ZW/Chr2R_ZS_ZH_ZW_Fst.csv")
Fst_ZS_ZH_ZW_Chr3L <- read.csv("Fst_ZS_ZH_ZW/Chr3L_ZS_ZH_ZW_Fst.csv")
Fst_ZS_ZH_ZW_Chr3R <- read.csv("Fst_ZS_ZH_ZW/Chr3R_ZS_ZH_ZW_Fst.csv")



#percentage subset functions
head_percent <- function(x, percent) {
  head(x, ceiling( nrow(x)*percent/100)) 
}

# last percent of a dataframe
tail_percent <- function(x, percent) {
  tail(x, ceiling( nrow(x)*percent/100)) 
}

#autosome dataframes
Fst_ZS_RAL_ZI_autosome <- bind_rows(Fst_ZS_RAL_ZI_Chr2L, Fst_ZS_RAL_ZI_Chr2R, Fst_ZS_RAL_ZI_Chr3L, Fst_ZS_RAL_ZI_Chr3R)
Fst_ZS_RAL_ZI_FR_SAfr_autosome <- bind_rows(Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, Fst_ZS_RAL_ZI_FR_SAfr_Chr3R)
Fst_Zim_RAL_ZI_autosome <- bind_rows(Fst_Zim_RAL_ZI_Chr2L, Fst_Zim_RAL_ZI_Chr2R, Fst_Zim_RAL_ZI_Chr3L, Fst_Zim_RAL_ZI_Chr3R)
Fst_ZH_RAL_ZI_autosome <- bind_rows(Fst_ZH_RAL_ZI_Chr2L, Fst_ZH_RAL_ZI_Chr2R, Fst_ZH_RAL_ZI_Chr3L, Fst_ZH_RAL_ZI_Chr3R)
Fst_ZW_RAL_ZI_autosome <- bind_rows(Fst_ZW_RAL_ZI_Chr2L, Fst_ZW_RAL_ZI_Chr2R, Fst_ZW_RAL_ZI_Chr3L, Fst_ZW_RAL_ZI_Chr3R)
Fst_ZS_ZH_ZW_autosome <- bind_rows(Fst_ZS_ZH_ZW_Chr2L, Fst_ZS_ZH_ZW_Chr2R, Fst_ZS_ZH_ZW_Chr3L, Fst_ZS_ZH_ZW_Chr3R)

#arrange
arranged_Fst_ZS_RAL_ZI_ChrX <- Fst_ZS_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_Chr2L <- Fst_ZS_RAL_ZI_Chr2L %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_Chr2R <- Fst_ZS_RAL_ZI_Chr2R %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_Chr3L <- Fst_ZS_RAL_ZI_Chr3L %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_Chr3R <- Fst_ZS_RAL_ZI_Chr3R %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_autosome <- Fst_ZS_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))


arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- Fst_ZS_RAL_ZI_FR_SAfr_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- Fst_ZS_RAL_ZI_FR_SAfr_Chr2L %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- Fst_ZS_RAL_ZI_FR_SAfr_Chr2R %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- Fst_ZS_RAL_ZI_FR_SAfr_Chr3L %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- Fst_ZS_RAL_ZI_FR_SAfr_Chr3R %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- Fst_ZS_RAL_ZI_FR_SAfr_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))


arranged_Fst_Zim_RAL_ZI_ChrX <- Fst_Zim_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))

arranged_Fst_Zim_RAL_ZI_Chr2L <- Fst_Zim_RAL_ZI_Chr2L %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))

arranged_Fst_Zim_RAL_ZI_Chr2R <- Fst_Zim_RAL_ZI_Chr2R %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))

arranged_Fst_Zim_RAL_ZI_Chr3L <- Fst_Zim_RAL_ZI_Chr3L %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))

arranged_Fst_Zim_RAL_ZI_Chr3R <- Fst_Zim_RAL_ZI_Chr3R %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))

arranged_Fst_Zim_RAL_ZI_autosome <- Fst_Zim_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_Zim.vs.RAL_ZI))


arranged_Fst_ZH_RAL_ZI_ChrX <- Fst_ZH_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_Chr2L <- Fst_ZH_RAL_ZI_Chr2L %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_Chr2R <- Fst_ZH_RAL_ZI_Chr2R %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_Chr3L <- Fst_ZH_RAL_ZI_Chr3L %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_Chr3R <- Fst_ZH_RAL_ZI_Chr3R %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_autosome <- Fst_ZH_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))


arranged_Fst_ZW_RAL_ZI_ChrX <- Fst_ZW_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_Chr2L <- Fst_ZW_RAL_ZI_Chr2L %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_Chr2R <- Fst_ZW_RAL_ZI_Chr2R %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_Chr3L <- Fst_ZW_RAL_ZI_Chr3L %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_Chr3R <- Fst_ZW_RAL_ZI_Chr3R %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_autosome <- Fst_ZW_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))


arranged_Fst_ZS_ZH_ZW_ChrX <- Fst_ZS_ZH_ZW_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_Chr2L <- Fst_ZS_ZH_ZW_Chr2L %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_Chr2R <- Fst_ZS_ZH_ZW_Chr2R %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_Chr3L <- Fst_ZS_ZH_ZW_Chr3L %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_Chr3R <- Fst_ZS_ZH_ZW_Chr3R %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_autosome <- Fst_ZS_ZH_ZW_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))


#variables for 1% subset files
chunk_Fst_ZS_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_ChrX, 1)
chunk_Fst_ZS_RAL_ZI_Chr2L <- head_percent(arranged_Fst_ZS_RAL_ZI_Chr2L, 1)
chunk_Fst_ZS_RAL_ZI_Chr2R <- head_percent(arranged_Fst_ZS_RAL_ZI_Chr2R, 1)
chunk_Fst_ZS_RAL_ZI_Chr3L <- head_percent(arranged_Fst_ZS_RAL_ZI_Chr3L, 1)
chunk_Fst_ZS_RAL_ZI_Chr3R <- head_percent(arranged_Fst_ZS_RAL_ZI_Chr3R, 1)
chunk_Fst_ZS_RAL_ZI_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_autosome, 1)

chunk_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, 1)
chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, 1)
chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, 1)
chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, 1)
chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, 1)
chunk_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome, 1)

chunk_Fst_Zim_RAL_ZI_ChrX <- head_percent(arranged_Fst_Zim_RAL_ZI_ChrX, 1)
chunk_Fst_Zim_RAL_ZI_Chr2L <- head_percent(arranged_Fst_Zim_RAL_ZI_Chr2L, 1)
chunk_Fst_Zim_RAL_ZI_Chr2R <- head_percent(arranged_Fst_Zim_RAL_ZI_Chr2R, 1)
chunk_Fst_Zim_RAL_ZI_Chr3L <- head_percent(arranged_Fst_Zim_RAL_ZI_Chr3L, 1)
chunk_Fst_Zim_RAL_ZI_Chr3R <- head_percent(arranged_Fst_Zim_RAL_ZI_Chr3R, 1)
chunk_Fst_Zim_RAL_ZI_autosome <- head_percent(arranged_Fst_Zim_RAL_ZI_autosome, 1)

chunk_Fst_ZH_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZH_RAL_ZI_ChrX, 1)
chunk_Fst_ZH_RAL_ZI_Chr2L <- head_percent(arranged_Fst_ZH_RAL_ZI_Chr2L, 1)
chunk_Fst_ZH_RAL_ZI_Chr2R <- head_percent(arranged_Fst_ZH_RAL_ZI_Chr2R, 1)
chunk_Fst_ZH_RAL_ZI_Chr3L <- head_percent(arranged_Fst_ZH_RAL_ZI_Chr3L, 1)
chunk_Fst_ZH_RAL_ZI_Chr3R <- head_percent(arranged_Fst_ZH_RAL_ZI_Chr3R, 1)
chunk_Fst_ZH_RAL_ZI_autosome <- head_percent(arranged_Fst_ZH_RAL_ZI_autosome, 1)

chunk_Fst_ZW_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZW_RAL_ZI_ChrX, 1)
chunk_Fst_ZW_RAL_ZI_Chr2L <- head_percent(arranged_Fst_ZW_RAL_ZI_Chr2L, 1)
chunk_Fst_ZW_RAL_ZI_Chr2R <- head_percent(arranged_Fst_ZW_RAL_ZI_Chr2R, 1)
chunk_Fst_ZW_RAL_ZI_Chr3L <- head_percent(arranged_Fst_ZW_RAL_ZI_Chr3L, 1)
chunk_Fst_ZW_RAL_ZI_Chr3R <- head_percent(arranged_Fst_ZW_RAL_ZI_Chr3R, 1)
chunk_Fst_ZW_RAL_ZI_autosome <- head_percent(arranged_Fst_ZW_RAL_ZI_autosome, 1)

chunk_Fst_ZS_ZH_ZW_ChrX <- head_percent(arranged_Fst_ZS_ZH_ZW_ChrX, 1)
chunk_Fst_ZS_ZH_ZW_Chr2L <- head_percent(arranged_Fst_ZS_ZH_ZW_Chr2L, 1)
chunk_Fst_ZS_ZH_ZW_Chr2R <- head_percent(arranged_Fst_ZS_ZH_ZW_Chr2R, 1)
chunk_Fst_ZS_ZH_ZW_Chr3L <- head_percent(arranged_Fst_ZS_ZH_ZW_Chr3L, 1)
chunk_Fst_ZS_ZH_ZW_Chr3R <- head_percent(arranged_Fst_ZS_ZH_ZW_Chr3R, 1)
chunk_Fst_ZS_ZH_ZW_autosome <- head_percent(arranged_Fst_ZS_ZH_ZW_autosome, 1)



#send to files
write.csv(chunk_Fst_ZS_RAL_ZI_ChrX, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_ChrX.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_Chr2L, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_Chr2L.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_Chr2R, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_Chr2R.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_Chr3L, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_Chr3L.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_Chr3R, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_Chr3R.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_autosome, file = "Fst_ZS_RAL_ZI/top0.01_Fst_ZS_RAL_ZI_autosome.csv")

write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R.csv")
write.csv(chunk_Fst_ZS_RAL_ZI_FR_SAfr_autosome, file = "Fst_ZS_RAL_ZI_FR_SAfr/top0.01_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv")

write.csv(chunk_Fst_Zim_RAL_ZI_ChrX, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_ChrX.csv")
write.csv(chunk_Fst_Zim_RAL_ZI_Chr2L, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_Chr2L.csv")
write.csv(chunk_Fst_Zim_RAL_ZI_Chr2R, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_Chr2R.csv")
write.csv(chunk_Fst_Zim_RAL_ZI_Chr3L, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_Chr3L.csv")
write.csv(chunk_Fst_Zim_RAL_ZI_Chr3R, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_Chr3R.csv")
write.csv(chunk_Fst_Zim_RAL_ZI_autosome, file = "Fst_Zim_RAL_ZI/top0.01_Fst_Zim_RAL_ZI_autosome.csv")

write.csv(chunk_Fst_ZH_RAL_ZI_ChrX, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_ChrX.csv")
write.csv(chunk_Fst_ZH_RAL_ZI_Chr2L, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_Chr2L.csv")
write.csv(chunk_Fst_ZH_RAL_ZI_Chr2R, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_Chr2R.csv")
write.csv(chunk_Fst_ZH_RAL_ZI_Chr3L, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_Chr3L.csv")
write.csv(chunk_Fst_ZH_RAL_ZI_Chr3R, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_Chr3R.csv")
write.csv(chunk_Fst_ZH_RAL_ZI_autosome, file = "Fst_ZH_RAL_ZI/top0.01_Fst_ZH_RAL_ZI_autosome.csv")

write.csv(chunk_Fst_ZW_RAL_ZI_ChrX, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_ChrX.csv")
write.csv(chunk_Fst_ZW_RAL_ZI_Chr2L, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_Chr2L.csv")
write.csv(chunk_Fst_ZW_RAL_ZI_Chr2R, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_Chr2R.csv")
write.csv(chunk_Fst_ZW_RAL_ZI_Chr3L, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_Chr3L.csv")
write.csv(chunk_Fst_ZW_RAL_ZI_Chr3R, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_Chr3R.csv")
write.csv(chunk_Fst_ZW_RAL_ZI_autosome, file = "Fst_ZW_RAL_ZI/top0.01_Fst_ZW_RAL_ZI_autosome.csv")

write.csv(chunk_Fst_ZS_ZH_ZW_ChrX, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_ChrX.csv")
write.csv(chunk_Fst_ZS_ZH_ZW_Chr2L, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_Chr2L.csv")
write.csv(chunk_Fst_ZS_ZH_ZW_Chr2R, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_Chr2R.csv")
write.csv(chunk_Fst_ZS_ZH_ZW_Chr3L, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_Chr3L.csv")
write.csv(chunk_Fst_ZS_ZH_ZW_Chr3R, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_Chr3R.csv")
write.csv(chunk_Fst_ZS_ZH_ZW_autosome, file = "Fst_ZS_ZH_ZW/top0.01_Fst_ZS_ZH_ZW_autosome.csv")









