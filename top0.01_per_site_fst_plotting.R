#loading packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(grid)


#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_downsampled/sliding_window_fst")

#reading files
windowed_Fst_ZS_RAL_ZI_ChrX <- read.csv("windowed_10kbp_ChrX_ZS_RAL_ZI_Fst.csv")
windowed_Fst_ZS_RAL_ZI_Chr2L <- read.csv("windowed_10kbp_Chr2L_ZS_RAL_ZI_Fst.csv")
windowed_Fst_ZS_RAL_ZI_Chr2R <- read.csv("windowed_10kbp_Chr2R_ZS_RAL_ZI_Fst.csv")
windowed_Fst_ZS_RAL_ZI_Chr3L <- read.csv("windowed_10kbp_Chr3L_ZS_RAL_ZI_Fst.csv")
windowed_Fst_ZS_RAL_ZI_Chr3R <- read.csv("windowed_10kbp_Chr3R_ZS_RAL_ZI_Fst.csv")

windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- read.csv("windowed_10kbp_ChrX_ZS_RAL_ZI_FR_SAfr_Fst.csv")
windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- read.csv("windowed_10kbp_Chr2L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- read.csv("windowed_10kbp_Chr2R_ZS_RAL_ZI_FR_SAfr_Fst.csv")
windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- read.csv("windowed_10kbp_Chr3L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- read.csv("windowed_10kbp_Chr3R_ZS_RAL_ZI_FR_SAfr_Fst.csv")

windowed_Fst_Zim_RAL_ZI_ChrX <- read.csv("windowed_10kbp_ChrX_Zim_RAL_ZI_Fst.csv")
windowed_Fst_Zim_RAL_ZI_Chr2L <- read.csv("windowed_10kbp_Chr2L_Zim_RAL_ZI_Fst.csv")
windowed_Fst_Zim_RAL_ZI_Chr2R <- read.csv("windowed_10kbp_Chr2R_Zim_RAL_ZI_Fst.csv")
windowed_Fst_Zim_RAL_ZI_Chr3L <- read.csv("windowed_10kbp_Chr3L_Zim_RAL_ZI_Fst.csv")
windowed_Fst_Zim_RAL_ZI_Chr3R <- read.csv("windowed_10kbp_Chr3R_Zim_RAL_ZI_Fst.csv")

windowed_Fst_ZH_RAL_ZI_ChrX <- read.csv("windowed_10kbp_ChrX_ZH_RAL_ZI_Fst.csv")

windowed_Fst_ZH_RAL_ZI_Chr2R <- read.csv("windowed_10kbp_Chr2R_ZH_RAL_ZI_Fst.csv")
windowed_Fst_ZH_RAL_ZI_Chr3L <- read.csv("windowed_10kbp_Chr3L_ZH_RAL_ZI_Fst.csv")
windowed_Fst_ZH_RAL_ZI_Chr3R <- read.csv("windowed_10kbp_Chr3R_ZH_RAL_ZI_Fst.csv")

windowed_Fst_ZW_RAL_ZI_ChrX <- read.csv("windowed_10kbp_ChrX_ZW_RAL_ZI_Fst.csv")
windowed_Fst_ZW_RAL_ZI_Chr2L <- read.csv("windowed_10kbp_Chr2L_ZW_RAL_ZI_Fst.csv")
windowed_Fst_ZW_RAL_ZI_Chr2R <- read.csv("windowed_10kbp_Chr2R_ZW_RAL_ZI_Fst.csv")
windowed_Fst_ZW_RAL_ZI_Chr3L <- read.csv("windowed_10kbp_Chr3L_ZW_RAL_ZI_Fst.csv")
windowed_Fst_ZW_RAL_ZI_Chr3R <- read.csv("windowed_10kbp_Chr3R_ZW_RAL_ZI_Fst.csv")

windowed_Fst_ZS_ZH_ZW_ChrX <- read.csv("windowed_10kbp_ChrX_ZS_ZH_ZW_Fst.csv")
windowed_Fst_ZS_ZH_ZW_Chr2L <- read.csv("windowed_10kbp_Chr2L_ZS_ZH_ZW_Fst.csv")
windowed_Fst_ZS_ZH_ZW_Chr2R <- read.csv("windowed_10kbp_Chr2R_ZS_ZH_ZW_Fst.csv")
windowed_Fst_ZS_ZH_ZW_Chr3L <- read.csv("windowed_10kbp_Chr3L_ZS_ZH_ZW_Fst.csv")
windowed_Fst_ZS_ZH_ZW_Chr3R <- read.csv("windowed_10kbp_Chr3R_ZS_ZH_ZW_Fst.csv")


#averaging
windowed_Fst_ZS_RAL_ZI_ChrX <- windowed_Fst_ZS_RAL_ZI_ChrX %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)

windowed_Fst_ZS_RAL_ZI_Chr2L <- windowed_Fst_ZS_RAL_ZI_Chr2L %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)

windowed_Fst_ZS_RAL_ZI_Chr2R <- windowed_Fst_ZS_RAL_ZI_Chr2R %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)

windowed_Fst_ZS_RAL_ZI_Chr3L <- windowed_Fst_ZS_RAL_ZI_Chr3L %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)

windowed_Fst_ZS_RAL_ZI_Chr3R <- windowed_Fst_ZS_RAL_ZI_Chr3R %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst) / 2)




windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)

windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)

windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)

windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)

windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R %>%
  mutate(Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr = (ZS.vs.RAL_Fst + ZS.vs.ZI_Fst + ZS.vs.FR_Fst + ZS.vs.SAfr_Fst) / 4)




windowed_Fst_Zim_RAL_ZI_ChrX <- windowed_Fst_Zim_RAL_ZI_ChrX %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)

windowed_Fst_Zim_RAL_ZI_Chr2L <- windowed_Fst_Zim_RAL_ZI_Chr2L %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)

windowed_Fst_Zim_RAL_ZI_Chr2R <- windowed_Fst_Zim_RAL_ZI_Chr2R %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)

windowed_Fst_Zim_RAL_ZI_Chr3L <- windowed_Fst_Zim_RAL_ZI_Chr3L %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)

windowed_Fst_Zim_RAL_ZI_Chr3R <- windowed_Fst_Zim_RAL_ZI_Chr3R %>%
  mutate(Avg_Fst_Zim.vs.RAL_ZI = (Zim.vs.RAL_Fst + Zim.vs.ZI_Fst) / 2)




windowed_Fst_ZH_RAL_ZI_ChrX <- windowed_Fst_ZH_RAL_ZI_ChrX %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)


windowed_Fst_ZH_RAL_ZI_Chr2R <- windowed_Fst_ZH_RAL_ZI_Chr2R %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)

windowed_Fst_ZH_RAL_ZI_Chr3L <- windowed_Fst_ZH_RAL_ZI_Chr3L %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)

windowed_Fst_ZH_RAL_ZI_Chr3R <- windowed_Fst_ZH_RAL_ZI_Chr3R %>%
  mutate(Avg_Fst_ZH.vs.RAL_ZI = (ZH.vs.RAL_Fst + ZH.vs.ZI_Fst) / 2)




windowed_Fst_ZW_RAL_ZI_ChrX <- windowed_Fst_ZW_RAL_ZI_ChrX %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)

windowed_Fst_ZW_RAL_ZI_Chr2L <- windowed_Fst_ZW_RAL_ZI_Chr2L %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)

windowed_Fst_ZW_RAL_ZI_Chr2R <- windowed_Fst_ZW_RAL_ZI_Chr2R %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)

windowed_Fst_ZW_RAL_ZI_Chr3L <- windowed_Fst_ZW_RAL_ZI_Chr3L %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)

windowed_Fst_ZW_RAL_ZI_Chr3R <- windowed_Fst_ZW_RAL_ZI_Chr3R %>%
  mutate(Avg_Fst_ZW.vs.RAL_ZI = (ZW.vs.RAL_Fst + ZW.vs.ZI_Fst) / 2)




windowed_Fst_ZS_ZH_ZW_ChrX <- windowed_Fst_ZS_ZH_ZW_ChrX %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW = (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)

windowed_Fst_ZS_ZH_ZW_Chr2L <- windowed_Fst_ZS_ZH_ZW_Chr2L %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW= (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)

windowed_Fst_ZS_ZH_ZW_Chr2R <- windowed_Fst_ZS_ZH_ZW_Chr2R %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW= (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)

windowed_Fst_ZS_ZH_ZW_Chr3L <- windowed_Fst_ZS_ZH_ZW_Chr3L %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW= (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)

windowed_Fst_ZS_ZH_ZW_Chr3R <- windowed_Fst_ZS_ZH_ZW_Chr3R %>%
  mutate(Avg_Fst_ZS.vs.ZH_ZW= (ZS.vs.ZH_Fst + ZS.vs.ZW_Fst) / 2)




#set directory
setwd("/Users/philipbaldassari/Desktop/zim-cos_downsampled/per_site_top0.01_fst")

#reading files
per_site_Fst_ZS_RAL_ZI_ChrX <- read.csv("per_site_top0.01_ChrX_ZS_RAL_ZI_Fst.csv")
per_site_Fst_ZS_RAL_ZI_Chr2L <- read.csv("per_site_top0.01_Chr2L_ZS_RAL_ZI_Fst.csv")
per_site_Fst_ZS_RAL_ZI_Chr2R <- read.csv("per_site_top0.01_Chr2R_ZS_RAL_ZI_Fst.csv")
per_site_Fst_ZS_RAL_ZI_Chr3L <- read.csv("per_site_top0.01_Chr3L_ZS_RAL_ZI_Fst.csv")
per_site_Fst_ZS_RAL_ZI_Chr3R <- read.csv("per_site_top0.01_Chr3R_ZS_RAL_ZI_Fst.csv")

per_site_Fst_ZS_RAL_ZI_FR_SAfr_ChrX <- read.csv("per_site_top0.01_ChrX_ZS_RAL_ZI_FR_SAfr_Fst.csv")
per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L <- read.csv("per_site_top0.01_Chr2L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R <- read.csv("per_site_top0.01_Chr2R_ZS_RAL_ZI_FR_SAfr_Fst.csv")
per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L <- read.csv("per_site_top0.01_Chr3L_ZS_RAL_ZI_FR_SAfr_Fst.csv")
per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R <- read.csv("per_site_top0.01_Chr3R_ZS_RAL_ZI_FR_SAfr_Fst.csv")

per_site_Fst_Zim_RAL_ZI_ChrX <- read.csv("per_site_top0.01_ChrX_Zim_RAL_ZI_Fst.csv")
per_site_Fst_Zim_RAL_ZI_Chr2L <- read.csv("per_site_top0.01_Chr2L_Zim_RAL_ZI_Fst.csv")
per_site_Fst_Zim_RAL_ZI_Chr2R <- read.csv("per_site_top0.01_Chr2R_Zim_RAL_ZI_Fst.csv")
per_site_Fst_Zim_RAL_ZI_Chr3L <- read.csv("per_site_top0.01_Chr3L_Zim_RAL_ZI_Fst.csv")
per_site_Fst_Zim_RAL_ZI_Chr3R <- read.csv("per_site_top0.01_Chr3R_Zim_RAL_ZI_Fst.csv")

per_site_Fst_ZH_RAL_ZI_ChrX <- read.csv("per_site_top0.01_ChrX_ZH_RAL_ZI_Fst.csv")

per_site_Fst_ZH_RAL_ZI_Chr2R <- read.csv("per_site_top0.01_Chr2R_ZH_RAL_ZI_Fst.csv")
per_site_Fst_ZH_RAL_ZI_Chr3L <- read.csv("per_site_top0.01_Chr3L_ZH_RAL_ZI_Fst.csv")
per_site_Fst_ZH_RAL_ZI_Chr3R <- read.csv("per_site_top0.01_Chr3R_ZH_RAL_ZI_Fst.csv")

per_site_Fst_ZW_RAL_ZI_ChrX <- read.csv("per_site_top0.01_ChrX_ZW_RAL_ZI_Fst.csv")
per_site_Fst_ZW_RAL_ZI_Chr2L <- read.csv("per_site_top0.01_Chr2L_ZW_RAL_ZI_Fst.csv")
per_site_Fst_ZW_RAL_ZI_Chr2R <- read.csv("per_site_top0.01_Chr2R_ZW_RAL_ZI_Fst.csv")
per_site_Fst_ZW_RAL_ZI_Chr3L <- read.csv("per_site_top0.01_Chr3L_ZW_RAL_ZI_Fst.csv")
per_site_Fst_ZW_RAL_ZI_Chr3R <- read.csv("per_site_top0.01_Chr3R_ZW_RAL_ZI_Fst.csv")

per_site_Fst_ZS_ZH_ZW_ChrX <- read.csv("per_site_top0.01_ChrX_ZS_ZH_ZW_Fst.csv")
per_site_Fst_ZS_ZH_ZW_Chr2L <- read.csv("per_site_top0.01_Chr2L_ZS_ZH_ZW_Fst.csv")
per_site_Fst_ZS_ZH_ZW_Chr2R <- read.csv("per_site_top0.01_Chr2R_ZS_ZH_ZW_Fst.csv")
per_site_Fst_ZS_ZH_ZW_Chr3L <- read.csv("per_site_top0.01_Chr3L_ZS_ZH_ZW_Fst.csv")
per_site_Fst_ZS_ZH_ZW_Chr3R <- read.csv("per_site_top0.01_Chr3R_ZS_ZH_ZW_Fst.csv")





##plotting
#ChrX
png(filename = "plots/ChrX_per_site_Fst_pairwise_comparisons.png", width = 1800, height = 900)
ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI, color = "ZS vs. RAL & ZI"), shape=15, size=1) + 
  geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr, color = "ZS vs. RAL, ZI, FR, & SAfr"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_Zim_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI, color = "Zim vs. RAL & ZI"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_ZH_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI, color = "ZH vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZW_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI, color = "ZW vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZS_ZH_ZW_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW, color = "ZS vs. ZH & ZW"), shape=15, size=1) +
  xlab(" ") + ylab("Mean Fst") + ggtitle("ChrX Fst (Top 1% Fst Sites)") + labs(colour=' ')
dev.off()


p1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("ChrX Fst (Top 1% Fst Sites)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(0.31,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(0.31,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(0.31,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=2, size=0.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(0.31,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(0.31,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
p6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=2, size=0.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(0.31,1) + theme_linedraw()
png(filename = "plots/stack_ChrX_windowed_Fst_pairwise_comparisons.png", width = 2000, height = 1000)
grid.newpage()
grid.draw(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4), ggplotGrob(p5), ggplotGrob(p6),size = "last"))
dev.off()


g1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_ChrX, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI), linetype="solid", size = 1.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("ChrX Fst (10kbp sliding window and top 1% SNPs)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(-0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
g2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_FR_SAfr_ChrX, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), linetype="solid", size = 1.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
g3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25)  +geom_line(data=windowed_Fst_Zim_RAL_ZI_ChrX, aes(x=window_start, y = Avg_Fst_Zim.vs.RAL_ZI), linetype="solid", size = 1.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
g4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=2, size=0.25)  +geom_line(data=windowed_Fst_ZH_RAL_ZI_ChrX, aes(x=window_start, y = Avg_Fst_ZH.vs.RAL_ZI), linetype="solid", size = 1.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
g5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_ChrX, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25) +geom_line(data=windowed_Fst_ZW_RAL_ZI_ChrX, aes(x=window_start, y = Avg_Fst_ZW.vs.RAL_ZI), linetype="solid", size = 1.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
g6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_ChrX, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=2, size=0.25) +geom_line(data=windowed_Fst_ZS_ZH_ZW_ChrX, aes(x=window_start, y = Avg_Fst_ZS.vs.ZH_ZW), linetype="solid", size = 1.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(-0.1,1) + theme_linedraw() + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
png(filename = "plots/stack_ChrX_per_site_Fst_pairwise_comparisons.png", width = 2000, height = 1500)
grid.newpage()
grid.draw(rbind(ggplotGrob(g1), ggplotGrob(g2), ggplotGrob(g3), ggplotGrob(g4), ggplotGrob(g5), ggplotGrob(g6),size = "last"))
dev.off()



#Chr2L
png(filename = "plots/Chr2L_per_site_Fst_pairwise_comparisons.png", width = 1800, height = 900)
ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI, color = "ZS vs. RAL & ZI"), shape=15, size=1) + 
  geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr, color = "ZS vs. RAL, ZI, FR, & SAfr"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI, color = "Zim vs. RAL & ZI"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI, color = "ZW vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW, color = "ZS vs. ZH & ZW"), shape=15, size=1) +
  xlab(" ") + ylab("Mean Fst") + ggtitle("Chr2L Fst (Top 1% Fst Sites)") + labs(colour=' ')
dev.off()


q1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr2L Fst (Top 1% Fst Sites)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(0,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
q2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
q3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
q4 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
q5 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(0,1) + theme_linedraw()
png(filename = "plots/stack_Chr2L_windowed_Fst_pairwise_comparisons.png", width = 2000, height = 1000)
grid.newpage()
grid.draw(rbind(ggplotGrob(q1), ggplotGrob(q2), ggplotGrob(q3), ggplotGrob(q4), ggplotGrob(q5), size = "last"))
dev.off()


w1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_Chr2L, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI), linetype="solid", size = 1.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr2L Fst (10kbp sliding window and top 1% SNPs)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(-0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
w2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2L, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), linetype="solid", size = 1.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
w3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25)  +geom_line(data=windowed_Fst_Zim_RAL_ZI_Chr2L, aes(x=window_start, y = Avg_Fst_Zim.vs.RAL_ZI), linetype="solid", size = 1.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
w4 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25) +geom_line(data=windowed_Fst_ZW_RAL_ZI_Chr2L, aes(x=window_start, y = Avg_Fst_ZW.vs.RAL_ZI), linetype="solid", size = 1.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
w5 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5) +geom_line(data=windowed_Fst_ZS_ZH_ZW_Chr2L, aes(x=window_start, y = Avg_Fst_ZS.vs.ZH_ZW), linetype="solid", size = 1.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(-0.1,1) + theme_linedraw() + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
png(filename = "plots/stack_Chr2L_per_site_Fst_pairwise_comparisons.png", width = 2000, height = 1500)
grid.newpage()
grid.draw(rbind(ggplotGrob(w1), ggplotGrob(w2), ggplotGrob(w3), ggplotGrob(w4), ggplotGrob(w5), size = "last"))
dev.off()



#Chr2R
png(filename = "plots/Chr2R_per_site_Fst_pairwise_comparisons.png", width = 1800, height = 900)
ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI, color = "ZS vs. RAL & ZI"), shape=15, size=1) + 
  geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr, color = "ZS vs. RAL, ZI, FR, & SAfr"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI, color = "Zim vs. RAL & ZI"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI, color = "ZH vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI, color = "ZW vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW, color = "ZS vs. ZH & ZW"), shape=15, size=1) +
  xlab(" ") + ylab("Mean Fst") + ggtitle("Chr2R Fst (Top 1% Fst Sites)") + labs(colour=' ')
dev.off()


t1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr2R Fst (Top 1% Fst Sites)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(0.2,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
t2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(0.2,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
t3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(0.2,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
t4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=2, size=0.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(0.2,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
t5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(0.2,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
t6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=2, size=0.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(0.2,1) + theme_linedraw()
png(filename = "plots/stack_Chr2R_windowed_Fst_pairwise_comparisons.png", width = 2000, height = 1000)
grid.newpage()
grid.draw(rbind(ggplotGrob(t1), ggplotGrob(t2), ggplotGrob(t3), ggplotGrob(t4), ggplotGrob(t5), ggplotGrob(t6),size = "last"))
dev.off()


y1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_Chr2R, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI), linetype="solid", size = 1.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr2R Fst (10kbp sliding window and top 1% SNPs)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(-0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
y2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr2R, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), linetype="solid", size = 1.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
y3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=2, size=0.25)  +geom_line(data=windowed_Fst_Zim_RAL_ZI_Chr2R, aes(x=window_start, y = Avg_Fst_Zim.vs.RAL_ZI), linetype="solid", size = 1.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
y4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=2, size=0.25)  +geom_line(data=windowed_Fst_ZH_RAL_ZI_Chr2R, aes(x=window_start, y = Avg_Fst_ZH.vs.RAL_ZI), linetype="solid", size = 1.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
y5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=2, size=0.25) +geom_line(data=windowed_Fst_ZW_RAL_ZI_Chr2R, aes(x=window_start, y = Avg_Fst_ZW.vs.RAL_ZI), linetype="solid", size = 1.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
y6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr2R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=2, size=0.25) +geom_line(data=windowed_Fst_ZS_ZH_ZW_Chr2R, aes(x=window_start, y = Avg_Fst_ZS.vs.ZH_ZW), linetype="solid", size = 1.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(-0.1,1) + theme_linedraw() + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
png(filename = "plots/stack_Chr2R_per_site_Fst_pairwise_comparisons.png", width = 2000, height = 1500)
grid.newpage()
grid.draw(rbind(ggplotGrob(y1), ggplotGrob(y2), ggplotGrob(y3), ggplotGrob(y4), ggplotGrob(y5), ggplotGrob(y6),size = "last"))
dev.off()



#Chr3L
png(filename = "plots/Chr3L_per_site_Fst_pairwise_comparisons.png", width = 1800, height = 900)
ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI, color = "ZS vs. RAL & ZI"), shape=15, size=1) + 
  geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr, color = "ZS vs. RAL, ZI, FR, & SAfr"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI, color = "Zim vs. RAL & ZI"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI, color = "ZH vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI, color = "ZW vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW, color = "ZS vs. ZH & ZW"), shape=15, size=1) +
  xlab(" ") + ylab("Mean Fst") + ggtitle("Chr3L Fst (Top 1% Fst Sites)") + labs(colour=' ')
dev.off()


o1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr3L Fst (Top 1% Fst Sites)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
o2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
o3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=19, size=1.5, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
o4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=19, size=1.5, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
o5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=19, size=1.5, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
o6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(0.1,1) + theme_linedraw()
png(filename = "plots/stack_Chr3L_windowed_Fst_pairwise_comparisons.png", width = 2000, height = 1000)
grid.newpage()
grid.draw(rbind(ggplotGrob(o1), ggplotGrob(o2), ggplotGrob(o3), ggplotGrob(o4), ggplotGrob(o5), ggplotGrob(o6),size = "last"))
dev.off()


l1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_Chr3L, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI), linetype="solid", size = 1.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr3L Fst (10kbp sliding window and top 1% SNPs)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(-0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
l2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3L, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), linetype="solid", size = 1.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
l3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=19, size=1.5)  +geom_line(data=windowed_Fst_Zim_RAL_ZI_Chr3L, aes(x=window_start, y = Avg_Fst_Zim.vs.RAL_ZI), linetype="solid", size = 1.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
l4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=19, size=1.5)  +geom_line(data=windowed_Fst_ZH_RAL_ZI_Chr3L, aes(x=window_start, y = Avg_Fst_ZH.vs.RAL_ZI), linetype="solid", size = 1.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
l5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=19, size=1.5) +geom_line(data=windowed_Fst_ZW_RAL_ZI_Chr3L, aes(x=window_start, y = Avg_Fst_ZW.vs.RAL_ZI), linetype="solid", size = 1.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
l6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3L, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5) +geom_line(data=windowed_Fst_ZS_ZH_ZW_Chr3L, aes(x=window_start, y = Avg_Fst_ZS.vs.ZH_ZW), linetype="solid", size = 1.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(-0.1,1) + theme_linedraw() + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
png(filename = "plots/stack_Chr3L_per_site_Fst_pairwise_comparisons.png", width = 2000, height = 1500)
grid.newpage()
grid.draw(rbind(ggplotGrob(l1), ggplotGrob(l2), ggplotGrob(l3), ggplotGrob(l4), ggplotGrob(l5), ggplotGrob(l6),size = "last"))
dev.off()



#Chr3R
png(filename = "plots/Chr3R_per_site_Fst_pairwise_comparisons.png", width = 1800, height = 900)
ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI, color = "ZS vs. RAL & ZI"), shape=15, size=1) + 
  geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr, color = "ZS vs. RAL, ZI, FR, & SAfr"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI, color = "Zim vs. RAL & ZI"), shape=15, size=1)  +
  geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI, color = "ZH vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI, color = "ZW vs. RAL & ZI"), shape=15, size=1) +
  geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW, color = "ZS vs. ZH & ZW"), shape=15, size=1) +
  xlab(" ") + ylab("Mean Fst") + ggtitle("Chr3R Fst (Top 1% Fst Sites)") + labs(colour=' ')
dev.off()


n1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr3R Fst (Top 1% Fst Sites)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(0,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
n2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
n3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=19, size=1, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
n4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=19, size=1.5, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
n5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=19, size=1.5, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(0,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank())
n6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(0,1) + theme_linedraw()
png(filename = "plots/stack_Chr3R_windowed_Fst_pairwise_comparisons.png", width = 2000, height = 1000)
grid.newpage()
grid.draw(rbind(ggplotGrob(n1), ggplotGrob(n2), ggplotGrob(n3), ggplotGrob(n4), ggplotGrob(n5), ggplotGrob(n6),size = "last"))
dev.off()


m1 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_Chr3R, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI), linetype="solid", size = 1.25, color = "azure4") + xlab(" ") + ylab("Mean Fst") + ggtitle(paste("Chr3R Fst (10kbp sliding window and top 1% SNPs)","\n\n", "ZS vs. RAL & ZI", sep=" "))  + ylim(-0.1,1) + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
m2 <- ggplot() + geom_point(data=per_site_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), shape=2, size=0.25) + geom_line(data=windowed_Fst_ZS_RAL_ZI_FR_SAfr_Chr3R, aes(x=window_start, y = Avg_Fst_ZS.vs.RAL_ZI_FR_SAfr), linetype="solid", size = 1.25, color = "blue1") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. RAL, ZI, FR, & SAfr") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
m3 <- ggplot() + geom_point(data=per_site_Fst_Zim_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_Zim.vs.RAL_ZI), shape=19, size=1)  +geom_line(data=windowed_Fst_Zim_RAL_ZI_Chr3R, aes(x=window_start, y = Avg_Fst_Zim.vs.RAL_ZI), linetype="solid", size = 1.25, color = "darkred") + xlab(" ") + ylab("Mean Fst") + ggtitle("Zim vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
m4 <- ggplot() + geom_point(data=per_site_Fst_ZH_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZH.vs.RAL_ZI), shape=19, size=1.5)  +geom_line(data=windowed_Fst_ZH_RAL_ZI_Chr3R, aes(x=window_start, y = Avg_Fst_ZH.vs.RAL_ZI), linetype="solid", size = 1.25, color = "chartreuse") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZH vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
m5 <- ggplot() + geom_point(data=per_site_Fst_ZW_RAL_ZI_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZW.vs.RAL_ZI), shape=19, size=1.5) +geom_line(data=windowed_Fst_ZW_RAL_ZI_Chr3R, aes(x=window_start, y = Avg_Fst_ZW.vs.RAL_ZI), linetype="solid", size = 1.25, color = "coral") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZW vs. RAL & ZI") + ylim(-0.1,1) + theme_linedraw() + theme_linedraw()+ theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
m6 <- ggplot() + geom_point(data=per_site_Fst_ZS_ZH_ZW_Chr3R, aes(x=Site_r6, y = Avg_Fst_ZS.vs.ZH_ZW), shape=19, size=1.5) +geom_line(data=windowed_Fst_ZS_ZH_ZW_Chr3R, aes(x=window_start, y = Avg_Fst_ZS.vs.ZH_ZW), linetype="solid", size = 1.25, color = "darkgoldenrod3") + xlab(" ") + ylab("Mean Fst") + ggtitle("ZS vs. ZH & ZW") + ylim(-0.1,1) + theme_linedraw() + geom_hline(yintercept=0, linetype="dotted", color="black", size=1.5)
png(filename = "plots/stack_Chr3R_per_site_Fst_pairwise_comparisons.png", width = 2000, height = 1500)
grid.newpage()
grid.draw(rbind(ggplotGrob(m1), ggplotGrob(m2), ggplotGrob(m3), ggplotGrob(m4), ggplotGrob(m5), ggplotGrob(m6),size = "last"))
dev.off()










