library(dplyr)
library(BioCircos)

setwd("~/Desktop/dros_speciation/manuscript/figures/making_figures/figure4")


######################################################################################################################################################



#opening gene data
ZSZH_specific_genes <- read.csv("ZS_ZH_specific_genes.tsv", sep = "\t")
ZSZH_specific_NS_genes <- read.csv("ZS_ZH_specific_NS_genes.tsv", sep = "\t")
ZSZH_specific_NS_top10_genes <- read.csv("ZS_ZH_specific_top10_NS_genes.tsv", sep = "\t")

ZS_99pct_mean <- read.csv("ZS_mean_99pct.tsv", sep = "\t")
ZH_99pct_mean <- read.csv("ZH_mean_99pct.tsv", sep = "\t")
ZR_99pct_mean <- read.csv("ZR_mean_99pct.tsv", sep = "\t")


#opening island data
ZS_islands_mean <- read.csv("ZS_mean_island_coordinates4circos.csv")
ZH_islands_mean <- read.csv("ZH_mean_island_coordinates4circos.csv")
ZR_islands_mean <- read.csv("ZR_mean_island_coordinates4circos.csv")
Cos_islands_mean <- read.csv("Cos_mean_island_coordinates4circos.csv")


#opening Fst data
win_mean_Fst <- read.csv("windowed_mean_5kbp_FST_averages.csv")



#genome
dmel_r6 = list("2L" = 23513712,
               "2R" = 25286936,
               "3L" = 28110227,
               "3R" = 32079331,
               "X" = 23542271)





#ZS and ZH specific gene tracks
tracks <- BioCircosSNPTrack("ZSZH_specific_genes", ZSZH_specific_genes$Chrom, ZSZH_specific_genes$Coordinate, rep(0.25,length(ZSZH_specific_genes$Chrom)),
                                     color = "black", size = 2.5, maxRadius = 1.25, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZSZH_specific_NS_genes", ZSZH_specific_NS_genes$Chrom, ZSZH_specific_NS_genes$Coordinate, rep(0.25,length(ZSZH_specific_NS_genes$Chrom)),
                            color = "orange", size = 3, maxRadius = 1.25, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZSZH_specific_top10_NS_genes", ZSZH_specific_NS_top10_genes$Chrom, ZSZH_specific_NS_top10_genes$Coordinate, rep(0.25,length(ZSZH_specific_NS_top10_genes$Chrom)),
                            color = "red", size = 4, maxRadius = 1.25, minRadius = 1, range = c(0,1))


#99th pct genes
tracks <- tracks + BioCircosSNPTrack("ZS_99", ZS_99pct_mean$Chrom, ZS_99pct_mean$Coordinate, rep(0.45,length(ZS_99pct_mean$Chrom)),
                                     color = "green", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZH_99", ZH_99pct_mean$Chrom, ZH_99pct_mean$Coordinate, rep(0.675,length(ZH_99pct_mean$Chrom)),
                                     color = "blue", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZR_99", ZR_99pct_mean$Chrom, ZR_99pct_mean$Coordinate, rep(0.9,length(ZR_99pct_mean$Chrom)),
                                     color = "purple", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))


#islands
tracks = tracks + BioCircosArcTrack('ZS_islands_mean', ZS_islands_mean$Chrom, ZS_islands_mean$Start, ZS_islands_mean$Stop,
                                    color = "green", maxRadius = 0.994, minRadius = 0.965)

tracks = tracks + BioCircosArcTrack('ZH_islands_mean', ZH_islands_mean$Chrom, ZH_islands_mean$Start, ZH_islands_mean$Stop,
                                    color = "blue", maxRadius = 0.96, minRadius = 0.93)

tracks = tracks + BioCircosArcTrack('ZR_islands_mean', ZR_islands_mean$Chrom, ZR_islands_mean$Start, ZR_islands_mean$Stop,
                                    color = "purple", maxRadius = 0.925, minRadius = 0.895)

tracks = tracks + BioCircosArcTrack('ZS_islands_mean', Cos_islands_mean$Chrom, Cos_islands_mean$Start, Cos_islands_mean$Stop,
                                    color = "#964B00", maxRadius = 0.89, minRadius = 0.866)


#Fst tracks
tracks <- tracks + BioCircosLineTrack("ZS_win_mean_fst", win_mean_Fst$Chrom, win_mean_Fst$window_start, win_mean_Fst$Win_FST_ZS,
                             color = "green", width = 0.5, maxRadius = 0.86, minRadius = 0.72, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("ZH_win_mean_fst", win_mean_Fst$Chrom, win_mean_Fst$window_start, win_mean_Fst$Win_FST_ZH,
                                      color = "blue", width = 0.5, maxRadius = 0.72, minRadius = 0.58, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("ZR_win_mean_fst", win_mean_Fst$Chrom, win_mean_Fst$window_start, win_mean_Fst$Win_FST_ZR,
                                      color = "purple", width = 0.5, maxRadius = 0.58, minRadius = 0.44, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("Cos_win_mean_fst", win_mean_Fst$Chrom, win_mean_Fst$window_start, win_mean_Fst$Win_FST_Cos,
                                      color = "#964B00", width = 0.5, maxRadius = 0.44, minRadius = 0.30, range = c(0,1))







#background tracks
tracks <- tracks + BioCircosBackgroundTrack("islands", fillColors = "snow",
                                            borderColors = "#000000", maxRadius = 1, minRadius = 0.86, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZS", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.86, minRadius = 0.72, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZH", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.72, minRadius = 0.58, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZR", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.58, minRadius = 0.44, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("Cos", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.44, minRadius = 0.30, borderSize = 3)






BioCircos(genome = dmel_r6,
          tracks,
          genomeLabelTextSize = "14pt",
          genomeTicksScale = 1000000,
          genomeTicksLen = 4,
          genomeTicksTextSize = 0,
          genomeLabelDy = 20,
          displayGenomeBorder = TRUE,
          genomeBorderSize = 3,
          chrPad = 0.04,
          genomeFillColor = c("white", "white", "white", "white", "white"))





######################################################################################################################################################



#opening gene data
ZSZH_specific_genes <- read.csv("ZS_ZH_specific_genes.tsv", sep = "\t")
ZSZH_specific_NS_genes <- read.csv("ZS_ZH_specific_NS_genes.tsv", sep = "\t")
ZSZH_specific_NS_top10_genes <- read.csv("ZS_ZH_specific_top10_NS_genes.tsv", sep = "\t")

ZS_99pct_max <- read.csv("ZS_max_99pct.tsv", sep = "\t")
ZH_99pct_max <- read.csv("ZH_max_99pct.tsv", sep = "\t")
ZR_99pct_max <- read.csv("ZR_max_99pct.tsv", sep = "\t")


#opening island data
ZS_islands_max <- read.csv("ZS_max_island_coordinates4circos.csv")
ZH_islands_max <- read.csv("ZH_max_island_coordinates4circos.csv")
ZR_islands_max <- read.csv("ZR_max_island_coordinates4circos.csv")
Cos_islands_max <- read.csv("Cos_max_island_coordinates4circos.csv")


#opening Fst data
win_max_Fst <- read.csv("windowed_max_5kbp_FST_averages.csv")
win_max_Fst <- win_max_Fst %>%
  add_row(Chrom = "X", 
          window_start = 23542271, 
          window_end = 23542272, 
          Win_FST_Zim = 1.0, 
          Win_FST_ZS = 1.0, 
          Win_FST_ZH = 1.0, 
          Win_FST_ZC = 1.0, 
          Win_FST_ZR = 1.0, 
          Win_FST_Cos = 1.0)



#genome
dmel_r6 = list("2L" = 23513712,
               "2R" = 25286936,
               "3L" = 28110227,
               "3R" = 32079331,
               "X" = 23542271)





#ZS and ZH specific gene tracks
tracks <- BioCircosSNPTrack("ZSZH_specific_genes", ZSZH_specific_genes$Chrom, ZSZH_specific_genes$Coordinate, rep(0.25,length(ZSZH_specific_genes$Chrom)),
                            color = "black", size = 2.5, maxRadius = 1.25, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZSZH_specific_NS_genes", ZSZH_specific_NS_genes$Chrom, ZSZH_specific_NS_genes$Coordinate, rep(0.25,length(ZSZH_specific_NS_genes$Chrom)),
                                     color = "orange", size = 3, maxRadius = 1.25, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZSZH_specific_top10_NS_genes", ZSZH_specific_NS_top10_genes$Chrom, ZSZH_specific_NS_top10_genes$Coordinate, rep(0.25,length(ZSZH_specific_NS_top10_genes$Chrom)),
                                     color = "red", size = 4, maxRadius = 1.25, minRadius = 1, range = c(0,1))


#99th pct genes
tracks <- tracks + BioCircosSNPTrack("ZS_99", ZS_99pct_max$Chrom, ZS_99pct_max$Coordinate, rep(0.45,length(ZS_99pct_max$Chrom)),
                                     color = "green", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZH_99", ZH_99pct_max$Chrom, ZH_99pct_max$Coordinate, rep(0.675,length(ZH_99pct_max$Chrom)),
                                     color = "blue", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZR_99", ZR_99pct_max$Chrom, ZR_99pct_max$Coordinate, rep(0.9,length(ZR_99pct_max$Chrom)),
                                     color = "purple", size = 4, maxRadius = 1.215, minRadius = 1, range = c(0,1))


#islands
tracks = tracks + BioCircosArcTrack('ZS_islands_max', ZS_islands_max$Chrom, ZS_islands_max$Start, ZS_islands_max$Stop,
                                    color = "green", maxRadius = 0.994, minRadius = 0.965)

tracks = tracks + BioCircosArcTrack('ZH_islands_max', ZH_islands_max$Chrom, ZH_islands_max$Start, ZH_islands_max$Stop,
                                    color = "blue", maxRadius = 0.96, minRadius = 0.93)

tracks = tracks + BioCircosArcTrack('ZR_islands_max', ZR_islands_max$Chrom, ZR_islands_max$Start, ZR_islands_max$Stop,
                                    color = "purple", maxRadius = 0.925, minRadius = 0.895)

tracks = tracks + BioCircosArcTrack('ZS_islands_max', Cos_islands_max$Chrom, Cos_islands_max$Start, Cos_islands_max$Stop,
                                    color = "#964B00", maxRadius = 0.89, minRadius = 0.866)


#Fst tracks
tracks <- tracks + BioCircosLineTrack("ZS_win_max_fst", win_max_Fst$Chrom, win_max_Fst$window_start, win_max_Fst$Win_FST_ZS,
                                      color = "green", width = 0.5, maxRadius = 0.86, minRadius = 0.72, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("ZH_win_max_fst", win_max_Fst$Chrom, win_max_Fst$window_start, win_max_Fst$Win_FST_ZH,
                                      color = "blue", width = 0.5, maxRadius = 0.72, minRadius = 0.58, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("ZR_win_max_fst", win_max_Fst$Chrom, win_max_Fst$window_start, win_max_Fst$Win_FST_ZR,
                                      color = "purple", width = 0.5, maxRadius = 0.58, minRadius = 0.44, range = c(0,1))

tracks <- tracks + BioCircosLineTrack("Cos_win_max_fst", win_max_Fst$Chrom, win_max_Fst$window_start, win_max_Fst$Win_FST_Cos,
                                      color = "#964B00", width = 0.5, maxRadius = 0.44, minRadius = 0.30, range = c(0,1))



#fixing the FST = 1 scaling problem
tracks <- tracks + BioCircosLineTrack("Zim.vs.SAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.86, minRadius = 0.72, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.WCAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.72, minRadius = 0.58, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.EAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.58, minRadius = 0.44, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.EAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.44, minRadius = 0.30, range = c(0,1))





#background tracks
tracks <- tracks + BioCircosBackgroundTrack("islands", fillColors = "snow",
                                            borderColors = "#000000", maxRadius = 1, minRadius = 0.86, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZS", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.86, minRadius = 0.72, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZH", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.72, minRadius = 0.58, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("ZR", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.58, minRadius = 0.44, borderSize = 3)

tracks <- tracks + BioCircosBackgroundTrack("Cos", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.44, minRadius = 0.30, borderSize = 3)






BioCircos(genome = dmel_r6,
          tracks,
          genomeLabelTextSize = "14pt",
          genomeTicksScale = 1000000,
          genomeTicksLen = 4,
          genomeTicksTextSize = 0,
          genomeLabelDy = 20,
          displayGenomeBorder = TRUE,
          genomeBorderSize = 3,
          chrPad = 0.04,
          genomeFillColor = c("white", "white", "white", "white", "white"))


















