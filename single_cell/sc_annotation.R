rm(list=ls())
setwd("\\ScRNA\\QC\\predata")

library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(patchwork)
library(viridis)
library(DoubletFinder)

load("/tn tLB GSE131907.RData")
sce.all <- subset(sce.all, Sample %in% c("LUNG_N01", "LUNG_N08","LUNG_N30","LUNG_N34",
                                         "LUNG_T08","LUNG_T30","LUNG_T34","BRONCHO_58"))

table(sce.all$Sample)

sample_list <- SplitObject(sce.all, split.by = "Sample")
samples <- names(sample_list)
samples
## QC + DoubletFinder
processed_list <- vector("list", length(samples))
names(processed_list) <- samples

for (i in seq_along(samples)) {
  cat("Processing sample:", samples[i], "\n")
  processed_list[[i]] <- make_seurat_object_and_doublet_removal(
    currentSample = sample_list[[i]],
    project_name  = samples[i],
    out_dir_root  = "\\QC\\predata"  
  )
}


seu_merged_filtered <- processed_list[[1]]

## merge
if (length(processed_list) > 1) {
  for (i in 2:length(processed_list)) {
    seu_merged_filtered <- merge(seu_merged_filtered, processed_list[[i]])
  }
}


seurat_qc_plots <- function(colon, sample_name){
  # Make some basic qc plots
  pdf(paste0("./seurat_nFeature_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nFeature_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_nCount_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("nCount_RNA"), ncol = 1, pt.size = 0.2))
  dev.off()
  
  pdf(paste0("./seurat_pMT_plots_", sample_name, ".pdf"), width = 40, height = 15)
  print(VlnPlot(colon, features = c("percent.mt"), ncol = 1, pt.size = 0.2))
  dev.off()
}

## QC
seurat_qc_plots(seu_merged_filtered, "all_samples_postDoublet")


########harmony
seu_merged_filtered <- NormalizeData(seu_merged_filtered)
seu_merged_filtered <- FindVariableFeatures(seu_merged_filtered)
seu_merged_filtered <- ScaleData(seu_merged_filtered, vars.to.regress = c("percent.mt"))
seu_merged_filtered <- RunPCA(seu_merged_filtered, verbose=F)
###Integrated with harmony
scRNA_harmony <- IntegrateLayers(object = seu_merged_filtered, 
                                 method = HarmonyIntegration, 
                                 orig.reduction = "pca", 
                                 new.reduction = "harmony",
                                 verbose = FALSE)
Layers(scRNA_harmony)
table(scRNA_harmony$Sample)
remove(seu_merged_filtered)
####JoinLayers
scRNA_harmony[["RNA"]] <- JoinLayers(scRNA_harmony[["RNA"]])


####reduction Umap
ElbowPlot(scRNA_harmony, ndims = 50)
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:20)

scRNA_harmony <- FindClusters(scRNA_harmony, resolution = seq(from = 0.1, to = 1.0, by = 0.1))
scRNA_harmony <- RunUMAP(scRNA_harmony, dims = 1:20, reduction = "harmony")

library(clustree)
clustree(scRNA_harmony)
DimPlot(scRNA_harmony , reduction = "umap",label = T) 
#####annotation####

#####Azimuth
library(AzimuthAPI)
scRNA_harmony<-CloudAzimuth(scRNA_harmony)
obj=scRNA_harmony
DimPlot(obj, group.by="azimuth_medium", label=T, label.size=5)
colnames(obj@meta.data)
library(dplyr)
for(i in levels(obj$RNA_snn_res.0.5)){
  print(paste("Cluster", i))
  obj@meta.data %>%
    filter(RNA_snn_res.0.5 == i)%>%
    group_by(azimuth_label)%>%
    summarize(count = n())%>%
    arrange(desc(count))%>%
    head(n =5)%>%
    print()
}


library(SingleR)
HumanPrimaryCellAtlasData=celldex::HumanPrimaryCellAtlasData()
save(HumanPrimaryCellAtlasData,file = "HumanPrimaryCellAtlasData.rdata")

load("HumanPrimaryCellAtlasData.rdata")
refdata <- HumanPrimaryCellAtlasData



testdata1=obj@assays[["RNA"]]@layers[["data"]]

testdata=LayerData(obj, assay="RNA", layer='data')

clusters <- obj@meta.data$RNA_snn_res.0.5

table(refdata$label.main)
###singler
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.main, 
                    method = "cluster", clusters = clusters, 
                    assay.type.test = "logcounts", assay.type.ref = "logcounts")
celltype_singleR = data.frame(ClusterID=rownames(cellpred), celltype=cellpred$labels)

###############dim20,resolution0.5
obj$group_manual <- NA

# 1) T / NK  cluster：0, 1, 4, 11, 19
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("0", "1", "4", "11", "19")] <- "T/NK cell"

# 2) B & Plasma：8（B cell），15（Plasma cell）
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("8", "15")] <- "B cell"

# 3) Myeloid：2, 3, 5, 6, 13, 20
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("2",  "5", "6", "13", "20")] <- "Myeloid cell"

# 4) Mast cell：7
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("7")] <- "Mast cell"

# 5) Fibroblast：9
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("9")] <- "Fibroblast"

# 6) Endothelial：12
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("12")] <- "Endothelial cell"

# 7) Epithelial ：10(AT2), 14(Club/Goblet/Basal), 17(AT1), 18(Ciliated)
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("10","3", "14", "17", "18")] <- "Epithelial cell"

# 8) Cycling：16
obj$group_manual[obj$RNA_snn_res.0.5 %in% c("16")] <- "Cycling cell"


table(obj$RNA_snn_res.0.5, obj$group_manual)
DimPlot(obj, group.by="group_manual", label=T, label.size=5)
DimPlot(obj, group.by="RNA_snn_res.0.5", label=T, label.size=5)


save(obj,file = "GSE131907_res0.5annotion.RData")

#####fig check####
rm(list=ls())

library(Seurat)
library(tidyverse)
library(dplyr)
library(patchwork)
library(harmony)
library(cowplot)
library(scRNAtoolVis)
library(devtools)
library(reshape2)
library(ggplot2)
library(tidydr)
library(RColorBrewer)

setwd("ScRNA\\annotion")
load("GSE131907_res0.5annotion.RData")
DimPlot(obj)
obj@assays$RNA@layers$counts

rownames(obj)[1:5]

####Umap#####
p1=DimPlot(obj, reduction = "umap", group.by = "group_manual",
           label = F, repel = F, shuffle = T ,pt.size = 0.01,label.size = 10) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),# panel.grid
        axis.title = element_text(face = 3,hjust = 0.03),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 16)
        #,plot.margin = margin(20, 50, 20, 50)
  ) +#axis.title
  #NoLegend()+
  coord_cartesian(clip = "off")+
  scale_color_manual(values = alpha(my36colors,0.3)) +  
  # NoLegend() + 
  scale_fill_manual(values = alpha(my36colors,0.3))
p1
ggsave(filename = "umap.pdf", width = 7, height = 5 )

p2=DimPlot(obj, reduction = "umap", group.by = "RNA_snn_res.0.5",
           label = T, repel = F, shuffle = T ,pt.size = 0.01,label.size = 6) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),# panel.grid
        axis.title = element_text(face = 3,hjust = 0.03),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 16)
        #,plot.margin = margin(20, 50, 20, 50)
  ) +#axis.title
  #NoLegend()+
  coord_cartesian(clip = "off")+
  scale_color_manual(values = alpha(my36colors,0.3)) +  
  # NoLegend() + 
  scale_fill_manual(values = alpha(my36colors,0.3))
p2
ggsave(filename = "umap cluster.pdf", width = 6, height = 5 )

######DotPlot###
genes_to_check = c(
  "EPCAM", "KRT19", "KRT18", "CDH1", ##epi
  "DCN","THY1", "COL1A1", "COL1A2", #Fibro
  "CDH5", "VWF","CLDN5","RAMP2",#end
  "CD3D", "CD3E", "TRAC",##T 
  "FGFBP2", #nk
  
  "LYZ","MARCO", "CD68", "CD163",##meloid
  "MS4A2","CPA3","TPSAB1", "TPSB2", ##Mast
  "IGHG3",  'MS4A1','BANK1', 'CD79A',#B cell
  
  "MKI67", "TOP2A","MCM5", "PCNA"#cycl
)
library(dplyr)

##extract cluster and celltype
cluster_labels <- obj@meta.data %>%
  dplyr::select(RNA_snn_res.0.5, group_manual) %>%
  dplyr::distinct() %>%                                   
  dplyr::arrange(as.numeric(as.character(RNA_snn_res.0.5))) %>%  
  dplyr::mutate(cluster_annot = paste0(RNA_snn_res.0.5,
                                       "\n",
                                       group_manual))


obj$cluster_annot <- paste0(obj$RNA_snn_res.0.5,
                            "\n",
                            obj$group_manual)


obj$cluster_annot <- factor(obj$cluster_annot,
                            levels = cluster_labels$cluster_annot)

##cluster_annot as group.by
p3 <- DotPlot(obj,
             features = unique(genes_to_check),
             assay = "RNA",
             group.by = "cluster_annot") +
  coord_flip() +
  scale_colour_gradientn(colors = brewer.pal(9, "YlGnBu")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black")
  ) +
  RotatedAxis()

p3

ggsave(plot=p3, filename="check cluster.pdf",width=9, height=9)

####violinPlot_gene####
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(RColorBrewer)

seu <- obj
DefaultAssay(seu) <- "RNA"

subtype_levels <- c(
  "B cell",
  "Cycling cell",
  "Endothelial cell",
  "Epithelial cell",
  "Fibroblast",
  "Mast cell",
  "Myeloid cell",
  "T/NK cell"
)

genes_to_plot <- c(
  "EDAR","GABRA2",  "GSTA3" ,"MAOB", "AKAP12", "ALDOA","F2","LDHA",  #"SLCO1B3",
  "SLC47A1","TRPA1" 
)

genes_exist <- genes_to_plot[genes_to_plot %in% rownames(seu)]

expr_df <- FetchData(
  object = seu,
  vars = c(genes_exist, "group_manual", "Sample_Origin"),
  slot = "data"
)

expr_df <- expr_df %>%
  mutate(
    Sample_Origin = ifelse(trimws(Sample_Origin) == "nLung", "nLung", "tLung")
  )

plot_df <- expr_df %>%
  rownames_to_column("cell") %>%
  pivot_longer(
    cols = all_of(genes_exist),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  mutate(
    subtype = group_manual,
    lung_group = Sample_Origin
  ) %>%
  filter(
    !is.na(subtype),
    !is.na(lung_group),
    subtype %in% subtype_levels,
    lung_group %in% c("nLung", "tLung")
  ) %>%
  mutate(
    subtype = factor(subtype, levels = subtype_levels),
    gene = factor(gene, levels = genes_exist),
    lung_group = factor(lung_group, levels = c("nLung", "tLung"))
  )

# 两组颜色
group_cols <- c(
  "nLung" = "#C44E52",
  "tLung" = "#4C72B0"
)

p1 <- ggplot(
  plot_df,
  aes(x = subtype, y = expression, fill = lung_group, color = lung_group)
) +
  geom_violin(
    position = position_dodge(width = 0.7),
    scale = "width",
    width = 0.6,
    trim = FALSE,
    linewidth = 0.25
  ) +
  stat_summary(
    aes(group = lung_group),
    fun = median,
    geom = "point",
    position = position_dodge(width = 0.7),
    shape = 21,
    size = 1.8,
    stroke = 0.3,
    fill = "white",
    color = "black"
  ) +
  facet_grid(
    rows = vars(gene),
    scales = "free_y"
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.08)),
    n.breaks = 3
  ) +
  scale_fill_manual(values = group_cols) +
  scale_color_manual(values = group_cols) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y.right = element_text(size = 16, color = "black", angle = 0),
    strip.text.y.left = element_text(size = 16, color = "black", angle = 0),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(fill = NA, color = "black", linewidth = 0),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.ticks.length.y.left = unit(0.05, "cm"),
    axis.text.x = element_text(
      angle = 35,
      color = "black",
      hjust = 1,
      vjust = 1,
      size = 18
    ),
    axis.ticks.x = element_line(),
    plot.margin = margin(0.5, 0, 0, 0, "cm")
  ) +
  labs(x = "", y = "")

p1

ggsave("violinPlot_gene.pdf",
       plot = p1,
       width = 10,
       height = 10,
       #units = "in",
       dpi = 600 )