rm(list = ls())
setwd("\\spatial\\EMTAB13530")

library(SpaCET)
library(cols4all)

# QC----
visiumPath <- "\\spatial\\EMTAB13530\\P10_T1\\"
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)
# calculate the QC metrics
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj)
# plot the QC metrics
SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "QualityControl", 
  spatialFeatures=c("UMI","Gene"),
  imageBg = F
) 

#deconvolve ST data----
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "LUAD", coreNo = 48)
# show the ST deconvolution results
SpaCET_obj@results$deconvolution$propMat[1:13,1:6]
# show the spatial distribution of malignant cells.
c4a_series()
c4a_palettes()

p1 <- SpaCET.visualize.spatialFeature(
  SpaCET_obj, 
  spatialType = "CellFraction", 
  spatialFeatures=c("Malignant"),
  imageBg = T) +
  scale_color_continuous_c4a_seq("hcl.blue_red2") +
  labs(color = "Fraction") +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12))
#theme(legend.position = "none",
#title = element_text(size = 12))
p1
ggsave(filename = "/P10_T1/P10Malignant.pdf", p1, height = 4, width = 4.6, dpi = 600)


# calculate the cell-cell colocalization.
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
#SpaCET.visualize.colocalization(SpaCET_obj)
dev.off()
# # calculate the L-R network score across ST spots.
# SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 48)
# SpaCET.visualize.spatialFeature(
#   SpaCET_obj,
#   spatialType = "LRNetworkScore",
#   spatialFeatures=c("Network_Score","Network_Score_pv")
# )

#Identify the Tumor-Stroma Interface----
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj, MalignantCutoff = 0.5)
p2 <- SpaCET.visualize.spatialFeature(SpaCET_obj, spatialType = "Interface", 
                                      spatialFeature = "Interface",
                                      imageBg = T) +
  theme(legend.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        title = element_text(size = 12))
theme(legend.position = "none",
      title = element_text(size = 12))
p2
ggsave(filename = "/P10_T1/P10Interface.pdf", p2, height = 4, width = 5, dpi = 600)

####LocationCompare_dotplot####

library(semla)
library(tibble)
library(ggplot2)
library(patchwork)
library(scico)
library(tidyr)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggpubr)

locations <- t(SpaCET_obj@results[["CCI"]][["interface"]])
barcodes <- SpaCET_obj@input[["metaData"]][["barcode"]]

str(locations)
head(locations)

loc_df <- data.frame(Barcodes = barcodes,
                     Location = locations)
rownames(loc_df) <- loc_df$Barcodes
colnames(loc_df)[2] <- "Location"
table(loc_df$Location)
save(loc_df, file = "\\P10_T1/P10_T1_LocInfo.RData")

#semla----
samples <- "PT10_T1/filtered_feature_bc_matrix.h5"
imgs <- "PT10_T1/spatial/tissue_hires_image.png"
spotfiles <- "PT10_T1/spatial/tissue_positions_list.csv"
json <- "PT10_T1/spatial/scalefactors_json.json"

infoTable <- tibble(samples, imgs, spotfiles, json)

p10_se <- ReadVisiumData(infoTable = infoTable)

##add loacation to meta.data

identical(rownames(loc_df), rownames(p10_se@meta.data))
p10_se$Location <- loc_df$Location

table(P10_se$Location)

genes <- c("MAOB","AKAP12","ALDOA","LDHA","SLC47A1","SLCO1B3")

library(RColorBrewer)
Idents(P10_se) <- P10_se$Location

dp <- DotPlot(P10_se, features = genes)$data
dp$pct.exp.clip <- pmax(dp$pct.exp, 0.1)


p_gene<-ggplot(dp, aes(features.plot, id)
) +
  geom_point(aes(size = pct.exp.clip, color = avg.exp.scaled)) +
  scale_size(range = c(2, 9)) +
  scale_colour_gradientn(colors=rev(RColorBrewer::brewer.pal(11, "RdYlBu"))) +
  #theme_pubr() +
  #theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()  +
  theme(axis.text.x = element_text(angle = 90,
                                   color = "black",
                                   size = 14
  ),
  axis.text.y = element_text(color = "black",
                             size = 16),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12)) +RotatedAxis()+
  labs(size="Percent Expressed", x=NULL, y=NULL)
p_gene

ggsave(
  filename = paste0("/P10_T1/P10_LocationCompare.pdf"),
  plot = p_gene,
  height = 4,
  width = 7,
  dpi = 600,
  device = cairo_pdf
)