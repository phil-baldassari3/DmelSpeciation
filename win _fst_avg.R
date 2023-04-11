library(dplyr)
library(tidyr)
library(tidyverse)


#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann/1kbp_windowed_GO_ann")


#reading files
windowed_Fst_ZS_RAL_ZI_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv")
windowed_Fst_ZS_RAL_ZI_A <- read.csv("windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv")

windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv")
windowed_Fst_ZS_RAL_ZI_FR_SAfr_A <- read.csv("windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv")


windowed_Fst_ZH_RAL_ZI_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv")
windowed_Fst_ZH_RAL_ZI_A <- read.csv("windowed_1kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv")


windowed_Fst_ZW_RAL_ZI_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv")
windowed_Fst_ZW_RAL_ZI_A <- read.csv("windowed_1kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv")


windowed_Fst_ZS_ZH_ZW_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv")
windowed_Fst_ZS_ZH_ZW_A <- read.csv("windowed_1kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv")


windowed_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv")
windowed_Fst_FR_vs_RAL_ZI_SAfr_A <- read.csv("windowed_1kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv")


windowed_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv")
windowed_Fst_RAL_vs_FR_ZI_SAfr_A <- read.csv("windowed_1kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv")


windowed_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv")
windowed_Fst_ZI_vs_RAL_FR_SAfr_A <- read.csv("windowed_1kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv")


windowed_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- read.csv("windowed_1kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv")
windowed_Fst_SAfr_vs_RAL_FR_ZI_A <- read.csv("windowed_1kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv")



#averaging
windowed_Fst_ZS_RAL_ZI_ChrX <- windowed_Fst_ZS_RAL_ZI_ChrX %>%
  mutate(Avg_Fst = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_ZS_RAL_ZI_A <- windowed_Fst_ZS_RAL_ZI_A %>%
  mutate(Avg_Fst = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)





windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX %>%
  mutate(Avg_Fst = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_ZS_RAL_ZI_FR_SAfr_A <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_A %>%
  mutate(Avg_Fst = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4) %>%
  relocate(Avg_Fst, .before=FBgn)





windowed_Fst_ZH_RAL_ZI_ChrX <- windowed_Fst_ZH_RAL_ZI_ChrX %>%
  mutate(Avg_Fst = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)


windowed_Fst_ZH_RAL_ZI_A <- windowed_Fst_ZH_RAL_ZI_A %>%
  mutate(Avg_Fst = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)






windowed_Fst_ZW_RAL_ZI_ChrX <- windowed_Fst_ZW_RAL_ZI_ChrX %>%
  mutate(Avg_Fst = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_ZW_RAL_ZI_A <- windowed_Fst_ZW_RAL_ZI_A %>%
  mutate(Avg_Fst = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)





windowed_Fst_ZS_ZH_ZW_ChrX <- windowed_Fst_ZS_ZH_ZW_ChrX %>%
  mutate(Avg_Fst = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_ZS_ZH_ZW_A <- windowed_Fst_ZS_ZH_ZW_A %>%
  mutate(Avg_Fst= (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2) %>%
  relocate(Avg_Fst, .before=FBgn)




windowed_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- windowed_Fst_FR_vs_RAL_ZI_SAfr_ChrX %>%
  mutate(Avg_Fst = (FR.vs.RAL_Fst + FR.vs.ZI_Fst + FR.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_FR_vs_RAL_ZI_SAfr_A <- windowed_Fst_FR_vs_RAL_ZI_SAfr_A %>%
  mutate(Avg_Fst= (FR.vs.RAL_Fst + FR.vs.ZI_Fst + FR.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)



windowed_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- windowed_Fst_RAL_vs_FR_ZI_SAfr_ChrX %>%
  mutate(Avg_Fst = (RAL.vs.FR_Fst + RAL.vs.ZI_Fst + RAL.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_RAL_vs_FR_ZI_SAfr_A <- windowed_Fst_RAL_vs_FR_ZI_SAfr_A %>%
  mutate(Avg_Fst= (RAL.vs.FR_Fst + RAL.vs.ZI_Fst + RAL.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)



windowed_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- windowed_Fst_ZI_vs_RAL_FR_SAfr_ChrX %>%
  mutate(Avg_Fst = (ZI.vs.RAL_Fst + ZI.vs.FR_Fst + ZI.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_ZI_vs_RAL_FR_SAfr_A <- windowed_Fst_ZI_vs_RAL_FR_SAfr_A %>%
  mutate(Avg_Fst= (ZI.vs.RAL_Fst + ZI.vs.FR_Fst + ZI.vs.SAfr_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)



windowed_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- windowed_Fst_SAfr_vs_RAL_FR_ZI_ChrX %>%
  mutate(Avg_Fst = (SAfr.vs.RAL_Fst + SAfr.vs.FR_Fst + SAfr.vs.ZI_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)

windowed_Fst_SAfr_vs_RAL_FR_ZI_A <- windowed_Fst_SAfr_vs_RAL_FR_ZI_A %>%
  mutate(Avg_Fst= (SAfr.vs.RAL_Fst + SAfr.vs.FR_Fst + SAfr.vs.ZI_Fst) / 3) %>%
  relocate(Avg_Fst, .before=FBgn)




#writing to files
write.csv(windowed_Fst_ZS_RAL_ZI_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv")
write.csv(windowed_Fst_ZS_RAL_ZI_A, file = "avg_windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv")

write.csv(windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv")
write.csv(windowed_Fst_ZS_RAL_ZI_FR_SAfr_A, file = "avg_windowed_1kbp_GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv")

write.csv(windowed_Fst_ZH_RAL_ZI_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv")
write.csv(windowed_Fst_ZH_RAL_ZI_A, file = "avg_windowed_1kbp_GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv")

write.csv(windowed_Fst_ZW_RAL_ZI_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv")
write.csv(windowed_Fst_ZW_RAL_ZI_A, file = "avg_windowed_1kbp_GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv")

write.csv(windowed_Fst_ZS_ZH_ZW_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv")
write.csv(windowed_Fst_ZS_ZH_ZW_A, file = "avg_windowed_1kbp_GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv")

write.csv(windowed_Fst_FR_vs_RAL_ZI_SAfr_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv")
write.csv(windowed_Fst_FR_vs_RAL_ZI_SAfr_A, file = "avg_windowed_1kbp_GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv")

write.csv(windowed_Fst_RAL_vs_FR_ZI_SAfr_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv")
write.csv(windowed_Fst_RAL_vs_FR_ZI_SAfr_A, file = "avg_windowed_1kbp_GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv")

write.csv(windowed_Fst_ZI_vs_RAL_FR_SAfr_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv")
write.csv(windowed_Fst_ZI_vs_RAL_FR_SAfr_A, file = "avg_windowed_1kbp_GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv")

write.csv(windowed_Fst_SAfr_vs_RAL_FR_ZI_ChrX, file = "avg_windowed_1kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv")
write.csv(windowed_Fst_SAfr_vs_RAL_FR_ZI_A, file = "avg_windowed_1kbp_GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv")








