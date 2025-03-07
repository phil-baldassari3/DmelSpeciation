library(dplyr)
library(BioCircos)

setwd("~/Desktop/dros_speciation/manuscript/figures/making_figures/figure3")


######################################################################################################################################################



#opening fst4pg data
Zim_Zam <- read.csv("Zim_Zam.csv")
Zim_SAfr <- read.csv("Zim_SAfr.csv")
Zim_WCAfr <- read.csv("Zim_WCAfr.csv")
Zim_EAfr <- read.csv("Zim_EAfr.csv")
Zim_OOA.OW <- read.csv("Zim_OOA.OW.csv")
Zim_OOA.NW <- read.csv("Zim_OOA.NW.csv")

#opening Gene data
Zim_Zam_genes <- read.csv("gene_coordinates/Zim_Zam_GO_genes.tsv", sep = "\t")
Zim_SAfr_genes <- read.csv("gene_coordinates/Zim_SAfr_GO_genes.tsv", sep = "\t")
Zim_WCAfr_genes <- read.csv("gene_coordinates/Zim_WCAfr_GO_genes.tsv", sep = "\t")
Zim_EAfr_genes <- read.csv("gene_coordinates/Zim_EAfr_GO_genes.tsv", sep = "\t")



#genome
dmel_r6 = list("2L" = 23513712,
               "2R" = 25286936,
               "3L" = 28110227,
               "3R" = 32079331,
               "X" = 23542271)

#line tracks
tracks <- BioCircosLineTrack("Zim.vs.Zam", Zim_Zam$Segment, Zim_Zam$Coordinate, Zim_Zam$Fst,
                             color = "red", width = 2, maxRadius = 1, minRadius = 0.875, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.SAfr", Zim_SAfr$Segment, Zim_SAfr$Coordinate, Zim_SAfr$Fst,
                                      color = "orange", width = 2, maxRadius = 0.875, minRadius = 0.75, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.WCAfr", Zim_WCAfr$Segment, Zim_WCAfr$Coordinate, Zim_WCAfr$Fst,
                                      color = "#cdcd00", width = 2, maxRadius = 0.75, minRadius = 0.625, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.EAfr", Zim_EAfr$Segment, Zim_EAfr$Coordinate, Zim_EAfr$Fst,
                                      color = "green", width = 2, maxRadius = 0.625, minRadius = 0.5, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.OOA.OW", Zim_OOA.OW$Segment, Zim_OOA.OW$Coordinate, Zim_OOA.OW$Fst,
                                      color = "blue", width = 2, maxRadius = 0.5, minRadius = 0.375, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.OOA.NW", Zim_OOA.NW$Segment, Zim_OOA.NW$Coordinate, Zim_OOA.NW$Fst,
                                      color = "purple", width = 2, maxRadius = 0.375, minRadius = 0.25, range = c(0,1))



#fixing the FST = 1 scaling problem
tracks <- tracks + BioCircosLineTrack("Zim.vs.Zam_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 1, minRadius = 0.875, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.SAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.875, minRadius = 0.75, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.WCAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.75, minRadius = 0.625, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.EAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.625, minRadius = 0.5, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.OOA.OW_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.5, minRadius = 0.375, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("Zim.vs.OOA.NW_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.375, minRadius = 0.25, range = c(0,1))



#background tracks
tracks <- tracks + BioCircosBackgroundTrack("Zam", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 1, minRadius = 0.875,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("SAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.875, minRadius = 0.75,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("WCAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.75, minRadius = 0.625,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("EAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.625, minRadius = 0.5,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("OOA.OW", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.5, minRadius = 0.375,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("OOA.NW", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.375, minRadius = 0.25,
                                            borderSize = 3)


#genes
tracks <- tracks + BioCircosSNPTrack("Zim.vs.Zam_genes", Zim_Zam_genes$Chrom, Zim_Zam_genes$Coordinate, rep(0.45,length(Zim_Zam_genes$Chrom)),
                             color = "red", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("Zim.vs.SAfr_genes", Zim_SAfr_genes$Chrom, Zim_SAfr_genes$Coordinate, rep(0.6,length(Zim_SAfr_genes$Chrom)),
                                      color = "orange", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("Zim.vs.WCAfr_genes", Zim_WCAfr_genes$Chrom, Zim_WCAfr_genes$Coordinate, rep(0.75,length(Zim_WCAfr_genes$Chrom)),
                                      color = "#cdcd00", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("Zim.vs.EAfr_genes", Zim_EAfr_genes$Chrom, Zim_EAfr_genes$Coordinate, rep(0.90,length(Zim_EAfr_genes$Chrom)),
                                      color = "green", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))


tracks <- tracks + BioCircosTextTrack('Title', 'Zim vs. Cos', x = -0.1, y = -1.3)


BioCircos(genome = dmel_r6,
          tracks,
          genomeLabelTextSize = "14pt",
          genomeTicksScale = 1000000,
          genomeTicksLen = 4,
          genomeTicksTextSize = 0,
          displayGenomeBorder = TRUE,
          genomeBorderSize = 3,
          chrPad = 0.04,
          genomeFillColor = c("white", "white", "white", "white", "white"))







######################################################################################################################################################






#opening fst4pg data
ZS_Zam <- read.csv("ZS_Zam.csv")
ZS_SAfr <- read.csv("ZS_SAfr.csv")
ZS_WCAfr <- read.csv("ZS_WCAfr.csv")
ZS_EAfr <- read.csv("ZS_EAfr.csv")
ZS_OOA.OW <- read.csv("ZS_OOA.OW.csv")
ZS_OOA.NW <- read.csv("ZS_OOA.NW.csv")

#opening Gene data
ZS_Zam_genes <- read.csv("gene_coordinates/ZS_Zam_GO_genes.tsv", sep = "\t")
ZS_SAfr_genes <- read.csv("gene_coordinates/ZS_SAfr_GO_genes.tsv", sep = "\t")
ZS_WCAfr_genes <- read.csv("gene_coordinates/ZS_WCAfr_GO_genes.tsv", sep = "\t")
ZS_EAfr_genes <- read.csv("gene_coordinates/ZS_EAfr_GO_genes.tsv", sep = "\t")
ZS_OOA.OW_genes <- read.csv("gene_coordinates/ZS_OOA-OW_GO_genes.tsv", sep = "\t")



#genome
dmel_r6 = list("2L" = 23513712,
               "2R" = 25286936,
               "3L" = 28110227,
               "3R" = 32079331,
               "X" = 23542271)

#line tracks
tracks <- BioCircosLineTrack("ZS.vs.Zam", ZS_Zam$Segment, ZS_Zam$Coordinate, ZS_Zam$Fst,
                             color = "red", width = 2, maxRadius = 1, minRadius = 0.875, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.SAfr", ZS_SAfr$Segment, ZS_SAfr$Coordinate, ZS_SAfr$Fst,
                                      color = "orange", width = 2, maxRadius = 0.875, minRadius = 0.75, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.WCAfr", ZS_WCAfr$Segment, ZS_WCAfr$Coordinate, ZS_WCAfr$Fst,
                                      color = "#cdcd00", width = 2, maxRadius = 0.75, minRadius = 0.625, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.EAfr", ZS_EAfr$Segment, ZS_EAfr$Coordinate, ZS_EAfr$Fst,
                                      color = "green", width = 2, maxRadius = 0.625, minRadius = 0.5, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.OOA.OW", ZS_OOA.OW$Segment, ZS_OOA.OW$Coordinate, ZS_OOA.OW$Fst,
                                      color = "blue", width = 2, maxRadius = 0.5, minRadius = 0.375, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.OOA.NW", ZS_OOA.NW$Segment, ZS_OOA.NW$Coordinate, ZS_OOA.NW$Fst,
                                      color = "purple", width = 2, maxRadius = 0.375, minRadius = 0.25, range = c(0,1))



#fixing the FST = 1 scaling problem
tracks <- tracks + BioCircosLineTrack("ZS.vs.Zam_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 1, minRadius = 0.875, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.SAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.875, minRadius = 0.75, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.WCAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.75, minRadius = 0.625, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.EAfr_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.625, minRadius = 0.5, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.OOA.OW_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.5, minRadius = 0.375, range = c(0,1))
tracks <- tracks + BioCircosLineTrack("ZS.vs.OOA.NW_fst1_correction", c("X", "X"), c(23542271, 23542272), c(0,1),
                                      color = "#000000", width = 2, maxRadius = 0.375, minRadius = 0.25, range = c(0,1))



#background tracks
tracks <- tracks + BioCircosBackgroundTrack("Zam", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 1, minRadius = 0.875,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("SAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.875, minRadius = 0.75,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("WCAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.75, minRadius = 0.625,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("EAfr", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.625, minRadius = 0.5,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("OOA.OW", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.5, minRadius = 0.375,
                                            borderSize = 3)
tracks <- tracks + BioCircosBackgroundTrack("OOA.NW", fillColors = "aliceblue",
                                            borderColors = "#000000", maxRadius = 0.375, minRadius = 0.25,
                                            borderSize = 3)


#genes
tracks <- tracks + BioCircosSNPTrack("ZS.vs.Zam_genes", ZS_Zam_genes$Chrom, ZS_Zam_genes$Coordinate, rep(0.45,length(ZS_Zam_genes$Chrom)),
                                     color = "red", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZS.vs.SAfr_genes", ZS_SAfr_genes$Chrom, ZS_SAfr_genes$Coordinate, rep(0.5625,length(ZS_SAfr_genes$Chrom)),
                                     color = "orange", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZS.vs.WCAfr_genes", ZS_WCAfr_genes$Chrom, ZS_WCAfr_genes$Coordinate, rep(0.675,length(ZS_WCAfr_genes$Chrom)),
                                     color = "#cdcd00", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZS.vs.EAfr_genes", ZS_EAfr_genes$Chrom, ZS_EAfr_genes$Coordinate, rep(0.7875,length(ZS_EAfr_genes$Chrom)),
                                     color = "green", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))

tracks <- tracks + BioCircosSNPTrack("ZS.vs.OOA.OW_genes", ZS_OOA.OW_genes$Chrom, ZS_OOA.OW_genes$Coordinate, rep(0.90,length(ZS_OOA.OW_genes$Chrom)),
                                     color = "blue", size = 2.5, maxRadius = 1.215, minRadius = 1, range = c(0,1))


tracks <- tracks + BioCircosTextTrack('Title', 'ZS vs. Cos', x = -0.1, y = -1.3)


BioCircos(genome = dmel_r6,
          tracks,
          genomeLabelTextSize = "14pt",
          genomeTicksScale = 1000000,
          genomeTicksLen = 4,
          genomeTicksTextSize = 0,
          displayGenomeBorder = TRUE,
          genomeBorderSize = 3,
          chrPad = 0.04,
          genomeFillColor = c("white", "white", "white", "white", "white"))






