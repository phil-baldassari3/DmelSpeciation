#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)


#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_gene_GO_ann")

#reading files
Fst_ZS_RAL_ZI_ChrX <- read.csv("GO_sites_genes_ZS_RAL_ZI_Fst_ChrX.csv")
Fst_ZS_RAL_ZI_autosome <- read.csv("GO_sites_genes_ZS_RAL_ZI_Fst_autosomes.csv")

Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- read.csv("GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_ChrX.csv")
Fst_ZS_RAL_ZI_FR_SAfr_autosome <- read.csv("GO_sites_genes_ZS_RAL_ZI_FR_SAfr_Fst_autosomes.csv")

Fst_ZH_RAL_ZI_ChrX <- read.csv("GO_sites_genes_ZH_RAL_ZI_Fst_ChrX.csv")
Fst_ZH_RAL_ZI_autosome <- read.csv("GO_sites_genes_ZH_RAL_ZI_Fst_autosomes.csv")

Fst_ZW_RAL_ZI_ChrX <- read.csv("GO_sites_genes_ZW_RAL_ZI_Fst_ChrX.csv")
Fst_ZW_RAL_ZI_autosome <- read.csv("GO_sites_genes_ZW_RAL_ZI_Fst_autosomes.csv")

Fst_ZS_ZH_ZW_ChrX <- read.csv("GO_sites_genes_ZS_ZH_ZW_Fst_ChrX.csv")
Fst_ZS_ZH_ZW_autosome <- read.csv("GO_sites_genes_ZS_ZH_ZW_Fst_autosomes.csv")

Fst_FR_vs_RAL_ZI_SAfr_ChrX <- read.csv("GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_ChrX.csv")
Fst_FR_vs_RAL_ZI_SAfr_autosome <- read.csv("GO_sites_genes_FR_vs_RAL_ZI_SAfr_Fst_autosomes.csv")

Fst_RAL_vs_FR_ZI_SAfr_ChrX <- read.csv("GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_ChrX.csv")
Fst_RAL_vs_FR_ZI_SAfr_autosome <- read.csv("GO_sites_genes_RAL_vs_FR_ZI_SAfr_Fst_autosomes.csv")

Fst_SAfr_vs_RAL_FR_ZI_ChrX <- read.csv("GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_ChrX.csv")
Fst_SAfr_vs_RAL_FR_ZI_autosome <- read.csv("GO_sites_genes_SAfr_vs_RAL_FR_ZI_Fst_autosomes.csv")

Fst_ZI_vs_RAL_FR_SAfr_ChrX <- read.csv("GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_ChrX.csv")
Fst_ZI_vs_RAL_FR_SAfr_autosome <- read.csv("GO_sites_genes_ZI_vs_RAL_FR_SAfr_Fst_autosomes.csv")



#percentage subset functions
head_percent <- function(x, percent) {
  head(x, ceiling( nrow(x)*percent/100)) 
}

# last percent of a dataframe
tail_percent <- function(x, percent) {
  tail(x, ceiling( nrow(x)*percent/100)) 
}





#arrange
arranged_Fst_ZS_RAL_ZI_ChrX <- Fst_ZS_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))

arranged_Fst_ZS_RAL_ZI_autosome <- Fst_ZS_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI))




arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- Fst_ZS_RAL_ZI_FR_SAfr_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))

arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- Fst_ZS_RAL_ZI_FR_SAfr_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr))




arranged_Fst_ZH_RAL_ZI_ChrX <- Fst_ZH_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))

arranged_Fst_ZH_RAL_ZI_autosome <- Fst_ZH_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZH.vs.RAL_ZI))



arranged_Fst_ZW_RAL_ZI_ChrX <- Fst_ZW_RAL_ZI_ChrX %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))

arranged_Fst_ZW_RAL_ZI_autosome <- Fst_ZW_RAL_ZI_autosome %>%
  arrange(desc(Avg_Fst_ZW.vs.RAL_ZI))



arranged_Fst_ZS_ZH_ZW_ChrX <- Fst_ZS_ZH_ZW_ChrX %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))

arranged_Fst_ZS_ZH_ZW_autosome <- Fst_ZS_ZH_ZW_autosome %>%
  arrange(desc(Avg_Fst_ZS.vs.ZH_ZW))



arranged_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- Fst_FR_vs_RAL_ZI_SAfr_ChrX %>%
  arrange(desc(Avg_Fst_FR.vs.RAL_ZI_SAfr))

arranged_Fst_FR_vs_RAL_ZI_SAfr_autosome <- Fst_FR_vs_RAL_ZI_SAfr_autosome %>%
  arrange(desc(Avg_Fst_FR.vs.RAL_ZI_SAfr))



arranged_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- Fst_RAL_vs_FR_ZI_SAfr_ChrX %>%
  arrange(desc(Avg_Fst_RAL.vs.FR_ZI_SAfr))

arranged_Fst_RAL_vs_FR_ZI_SAfr_autosome <- Fst_RAL_vs_FR_ZI_SAfr_autosome %>%
  arrange(desc(Avg_Fst_RAL.vs.FR_ZI_SAfr))



arranged_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- Fst_SAfr_vs_RAL_FR_ZI_ChrX %>%
  arrange(desc(Avg_Fst_SAfr.vs.RAL_FR_ZI))

arranged_Fst_SAfr_vs_RAL_FR_ZI_autosome <- Fst_SAfr_vs_RAL_FR_ZI_autosome %>%
  arrange(desc(Avg_Fst_SAfr.vs.RAL_FR_ZI))



arranged_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- Fst_ZI_vs_RAL_FR_SAfr_ChrX %>%
  arrange(desc(Avg_Fst_ZI.vs.RAL_FR_SAfr))

arranged_Fst_ZI_vs_RAL_FR_SAfr_autosome <- Fst_ZI_vs_RAL_FR_SAfr_autosome %>%
  arrange(desc(Avg_Fst_ZI.vs.RAL_FR_SAfr))


#variables for 1% subset files
chunk5_Fst_ZS_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_ChrX, 5)
chunk5_Fst_ZS_RAL_ZI_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_autosome, 5)

chunk5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, 5)
chunk5_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome, 5)

chunk5_Fst_ZH_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZH_RAL_ZI_ChrX, 5)
chunk5_Fst_ZH_RAL_ZI_autosome <- head_percent(arranged_Fst_ZH_RAL_ZI_autosome, 5)

chunk5_Fst_ZW_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZW_RAL_ZI_ChrX, 5)
chunk5_Fst_ZW_RAL_ZI_autosome <- head_percent(arranged_Fst_ZW_RAL_ZI_autosome, 5)

chunk5_Fst_ZS_ZH_ZW_ChrX <- head_percent(arranged_Fst_ZS_ZH_ZW_ChrX, 5)
chunk5_Fst_ZS_ZH_ZW_autosome <- head_percent(arranged_Fst_ZS_ZH_ZW_autosome, 5)

chunk5_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_ChrX, 5)
chunk5_Fst_FR_vs_RAL_ZI_SAfr_autosome <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_autosome, 5)

chunk5_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_ChrX, 5)
chunk5_Fst_RAL_vs_FR_ZI_SAfr_autosome <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_autosome, 5)

chunk5_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_ChrX, 5)
chunk5_Fst_SAfr_vs_RAL_FR_ZI_autosome <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_autosome, 5)

chunk5_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_ChrX, 5)
chunk5_Fst_ZI_vs_RAL_FR_SAfr_autosome <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_autosome, 5)







chunk1_Fst_ZS_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_ChrX, 1)
chunk1_Fst_ZS_RAL_ZI_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_autosome, 1)

chunk1_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, 1)
chunk1_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome, 1)

chunk1_Fst_ZH_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZH_RAL_ZI_ChrX, 1)
chunk1_Fst_ZH_RAL_ZI_autosome <- head_percent(arranged_Fst_ZH_RAL_ZI_autosome, 1)

chunk1_Fst_ZW_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZW_RAL_ZI_ChrX, 1)
chunk1_Fst_ZW_RAL_ZI_autosome <- head_percent(arranged_Fst_ZW_RAL_ZI_autosome, 1)

chunk1_Fst_ZS_ZH_ZW_ChrX <- head_percent(arranged_Fst_ZS_ZH_ZW_ChrX, 1)
chunk1_Fst_ZS_ZH_ZW_autosome <- head_percent(arranged_Fst_ZS_ZH_ZW_autosome, 1)

chunk1_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_ChrX, 1)
chunk1_Fst_FR_vs_RAL_ZI_SAfr_autosome <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_autosome, 1)

chunk1_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_ChrX, 1)
chunk1_Fst_RAL_vs_FR_ZI_SAfr_autosome <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_autosome, 1)

chunk1_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_ChrX, 1)
chunk1_Fst_SAfr_vs_RAL_FR_ZI_autosome <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_autosome, 1)

chunk1_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_ChrX, 1)
chunk1_Fst_ZI_vs_RAL_FR_SAfr_autosome <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_autosome, 1)






chunk0.5_Fst_ZS_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_ChrX, 0.5)
chunk0.5_Fst_ZS_RAL_ZI_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_autosome, 0.5)

chunk0.5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, 0.5)
chunk0.5_Fst_ZS_RAL_ZI_FR_SAfr_autosome <- head_percent(arranged_Fst_ZS_RAL_ZI_FR_SAfr_autosome, 0.5)

chunk0.5_Fst_ZH_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZH_RAL_ZI_ChrX, 0.5)
chunk0.5_Fst_ZH_RAL_ZI_autosome <- head_percent(arranged_Fst_ZH_RAL_ZI_autosome, 0.5)

chunk0.5_Fst_ZW_RAL_ZI_ChrX <- head_percent(arranged_Fst_ZW_RAL_ZI_ChrX, 0.5)
chunk0.5_Fst_ZW_RAL_ZI_autosome <- head_percent(arranged_Fst_ZW_RAL_ZI_autosome, 0.5)

chunk0.5_Fst_ZS_ZH_ZW_ChrX <- head_percent(arranged_Fst_ZS_ZH_ZW_ChrX, 0.5)
chunk0.5_Fst_ZS_ZH_ZW_autosome <- head_percent(arranged_Fst_ZS_ZH_ZW_autosome, 0.5)

chunk0.5_Fst_FR_vs_RAL_ZI_SAfr_ChrX <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_ChrX, 0.5)
chunk0.5_Fst_FR_vs_RAL_ZI_SAfr_autosome <- head_percent(arranged_Fst_FR_vs_RAL_ZI_SAfr_autosome, 0.5)

chunk0.5_Fst_RAL_vs_FR_ZI_SAfr_ChrX <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_ChrX, 0.5)
chunk0.5_Fst_RAL_vs_FR_ZI_SAfr_autosome <- head_percent(arranged_Fst_RAL_vs_FR_ZI_SAfr_autosome, 0.5)

chunk0.5_Fst_SAfr_vs_RAL_FR_ZI_ChrX <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_ChrX, 0.5)
chunk0.5_Fst_SAfr_vs_RAL_FR_ZI_autosome <- head_percent(arranged_Fst_SAfr_vs_RAL_FR_ZI_autosome, 0.5)

chunk0.5_Fst_ZI_vs_RAL_FR_SAfr_ChrX <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_ChrX, 0.5)
chunk0.5_Fst_ZI_vs_RAL_FR_SAfr_autosome <- head_percent(arranged_Fst_ZI_vs_RAL_FR_SAfr_autosome, 0.5)



#send to files
write.csv(chunk5_Fst_ZS_RAL_ZI_ChrX, file = "top5_Fst_ZS_RAL_ZI_ChrX.csv")
write.csv(chunk5_Fst_ZS_RAL_ZI_autosome, file = "top5_Fst_ZS_RAL_ZI_autosome.csv")

write.csv(chunk5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, file = "top5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv")
write.csv(chunk5_Fst_ZS_RAL_ZI_FR_SAfr_autosome, file = "top5_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv")

write.csv(chunk5_Fst_ZH_RAL_ZI_ChrX, file = "top5_Fst_ZH_RAL_ZI_ChrX.csv")
write.csv(chunk5_Fst_ZH_RAL_ZI_autosome, file = "top5_Fst_ZH_RAL_ZI_autosome.csv")

write.csv(chunk5_Fst_ZW_RAL_ZI_ChrX, file = "top5_Fst_ZW_RAL_ZI_ChrX.csv")
write.csv(chunk5_Fst_ZW_RAL_ZI_autosome, file = "top5_Fst_ZW_RAL_ZI_autosome.csv")

write.csv(chunk5_Fst_ZS_ZH_ZW_ChrX, file = "top5_Fst_ZS_ZH_ZW_ChrX.csv")
write.csv(chunk5_Fst_ZS_ZH_ZW_autosome, file = "top5_Fst_ZS_ZH_ZW_autosome.csv")

write.csv(chunk5_Fst_FR_vs_RAL_ZI_SAfr_ChrX, file = "top5_Fst_FR_vs_RAL_ZI_SAfr_ChrX.csv")
write.csv(chunk5_Fst_FR_vs_RAL_ZI_SAfr_autosome, file = "top5_Fst_FR_vs_RAL_ZI_SAfr_autosome.csv")

write.csv(chunk5_Fst_RAL_vs_FR_ZI_SAfr_ChrX, file = "top5_Fst_RAL_vs_FR_ZI_SAfr_ChrX.csv")
write.csv(chunk5_Fst_RAL_vs_FR_ZI_SAfr_autosome, file = "top5_Fst_RAL_vs_FR_ZI_SAfr_autosome.csv")

write.csv(chunk5_Fst_SAfr_vs_RAL_FR_ZI_ChrX, file = "top5_Fst_SAfr_vs_RAL_FR_ZI_ChrX.csv")
write.csv(chunk5_Fst_SAfr_vs_RAL_FR_ZI_autosome, file = "top5_Fst_SAfr_vs_RAL_FR_ZI_autosome.csv")

write.csv(chunk5_Fst_ZI_vs_RAL_FR_SAfr_ChrX, file = "top5_Fst_ZI_vs_RAL_FR_SAfr_ChrX.csv")
write.csv(chunk5_Fst_ZI_vs_RAL_FR_SAfr_autosome, file = "top5_Fst_ZI_vs_RAL_FR_SAfr_autosome.csv")



write.csv(chunk1_Fst_ZS_RAL_ZI_ChrX, file = "top1_Fst_ZS_RAL_ZI_ChrX.csv")
write.csv(chunk1_Fst_ZS_RAL_ZI_autosome, file = "top1_Fst_ZS_RAL_ZI_autosome.csv")

write.csv(chunk1_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, file = "top1_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv")
write.csv(chunk1_Fst_ZS_RAL_ZI_FR_SAfr_autosome, file = "top1_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv")

write.csv(chunk1_Fst_ZH_RAL_ZI_ChrX, file = "top1_Fst_ZH_RAL_ZI_ChrX.csv")
write.csv(chunk1_Fst_ZH_RAL_ZI_autosome, file = "top1_Fst_ZH_RAL_ZI_autosome.csv")

write.csv(chunk1_Fst_ZW_RAL_ZI_ChrX, file = "top1_Fst_ZW_RAL_ZI_ChrX.csv")
write.csv(chunk1_Fst_ZW_RAL_ZI_autosome, file = "top1_Fst_ZW_RAL_ZI_autosome.csv")

write.csv(chunk1_Fst_ZS_ZH_ZW_ChrX, file = "top1_Fst_ZS_ZH_ZW_ChrX.csv")
write.csv(chunk1_Fst_ZS_ZH_ZW_autosome, file = "top1_Fst_ZS_ZH_ZW_autosome.csv")

write.csv(chunk1_Fst_FR_vs_RAL_ZI_SAfr_ChrX, file = "top1_Fst_FR_vs_RAL_ZI_SAfr_ChrX.csv")
write.csv(chunk1_Fst_FR_vs_RAL_ZI_SAfr_autosome, file = "top1_Fst_FR_vs_RAL_ZI_SAfr_autosome.csv")

write.csv(chunk1_Fst_RAL_vs_FR_ZI_SAfr_ChrX, file = "top1_Fst_RAL_vs_FR_ZI_SAfr_ChrX.csv")
write.csv(chunk1_Fst_RAL_vs_FR_ZI_SAfr_autosome, file = "top1_Fst_RAL_vs_FR_ZI_SAfr_autosome.csv")

write.csv(chunk1_Fst_SAfr_vs_RAL_FR_ZI_ChrX, file = "top1_Fst_SAfr_vs_RAL_FR_ZI_ChrX.csv")
write.csv(chunk1_Fst_SAfr_vs_RAL_FR_ZI_autosome, file = "top1_Fst_SAfr_vs_RAL_FR_ZI_autosome.csv")

write.csv(chunk1_Fst_ZI_vs_RAL_FR_SAfr_ChrX, file = "top1_Fst_ZI_vs_RAL_FR_SAfr_ChrX.csv")
write.csv(chunk1_Fst_ZI_vs_RAL_FR_SAfr_autosome, file = "top1_Fst_ZI_vs_RAL_FR_SAfr_autosome.csv")



write.csv(chunk0.5_Fst_ZS_RAL_ZI_ChrX, file = "top0.5_Fst_ZS_RAL_ZI_ChrX.csv")
write.csv(chunk0.5_Fst_ZS_RAL_ZI_autosome, file = "top0.5_Fst_ZS_RAL_ZI_autosome.csv")

write.csv(chunk0.5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, file = "top0.5_Fst_ZS_RAL_ZI_FR_SAfr_ChrX.csv")
write.csv(chunk0.5_Fst_ZS_RAL_ZI_FR_SAfr_autosome, file = "top0.5_Fst_ZS_RAL_ZI_FR_SAfr_autosome.csv")

write.csv(chunk0.5_Fst_ZH_RAL_ZI_ChrX, file = "top0.5_Fst_ZH_RAL_ZI_ChrX.csv")
write.csv(chunk0.5_Fst_ZH_RAL_ZI_autosome, file = "top0.5_Fst_ZH_RAL_ZI_autosome.csv")

write.csv(chunk0.5_Fst_ZW_RAL_ZI_ChrX, file = "top0.5_Fst_ZW_RAL_ZI_ChrX.csv")
write.csv(chunk0.5_Fst_ZW_RAL_ZI_autosome, file = "top0.5_Fst_ZW_RAL_ZI_autosome.csv")

write.csv(chunk0.5_Fst_ZS_ZH_ZW_ChrX, file = "top0.5_Fst_ZS_ZH_ZW_ChrX.csv")
write.csv(chunk0.5_Fst_ZS_ZH_ZW_autosome, file = "top0.5_Fst_ZS_ZH_ZW_autosome.csv")

write.csv(chunk0.5_Fst_FR_vs_RAL_ZI_SAfr_ChrX, file = "top0.5_Fst_FR_vs_RAL_ZI_SAfr_ChrX.csv")
write.csv(chunk0.5_Fst_FR_vs_RAL_ZI_SAfr_autosome, file = "top0.5_Fst_FR_vs_RAL_ZI_SAfr_autosome.csv")

write.csv(chunk0.5_Fst_RAL_vs_FR_ZI_SAfr_ChrX, file = "top0.5_Fst_RAL_vs_FR_ZI_SAfr_ChrX.csv")
write.csv(chunk0.5_Fst_RAL_vs_FR_ZI_SAfr_autosome, file = "top0.5_Fst_RAL_vs_FR_ZI_SAfr_autosome.csv")

write.csv(chunk0.5_Fst_SAfr_vs_RAL_FR_ZI_ChrX, file = "top0.5_Fst_SAfr_vs_RAL_FR_ZI_ChrX.csv")
write.csv(chunk0.5_Fst_SAfr_vs_RAL_FR_ZI_autosome, file = "top0.5_Fst_SAfr_vs_RAL_FR_ZI_autosome.csv")

write.csv(chunk0.5_Fst_ZI_vs_RAL_FR_SAfr_ChrX, file = "top0.5_Fst_ZI_vs_RAL_FR_SAfr_ChrX.csv")
write.csv(chunk0.5_Fst_ZI_vs_RAL_FR_SAfr_autosome, file = "top0.5_Fst_ZI_vs_RAL_FR_SAfr_autosome.csv")




