rm(list=ls())
setwd("\\ScRNA\\epi")

#####epi_annotion######
load("/GSE131907_res0.5annotion.RData")
#extact Epithelial cell
table(Idents(obj))
Idents(obj) <- "group_manual"
Epi_sce <- subset(obj, idents = c("Epithelial cell"))
table(Epi_sce@meta.data$group_manual)
rownames(Epi_sce)
Epi_sce = CreateSeuratObject(counts = GetAssayData(Epi_sce, assay="RNA", layer='counts'),  
                             meta.data = Epi_sce@meta.data)  # keep meta.data
Epi_sce@assays$RNA@layers$counts
#Normalize the data, identify highly variable genes, perform normalization, and conduct PCA dimensionality reduction
Epi_sce <- NormalizeData(Epi_sce) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)

Epi_sce
ElbowPlot(Epi_sce, ndims=50, reduction="pca")
Epi_sce <- FindNeighbors(Epi_sce, reduction = "pca", dims = 1:20)#reduction="harmony"
Epi_sce = FindClusters(Epi_sce,resolution = 0.2)
table(Epi_sce@meta.data$seurat_clusters)
Epi_sce <- RunUMAP(Epi_sce, reduction = "pca", dims = 1:20)#reduction="harmony"

# Inherit the UMAP coordinates of the original object
Epi_sce@reductions$umap <- obj@reductions$umap

#EPI_umap without annotation
Epi_sce$group
plot1 =DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE) 
plot1
plot2 = DimPlot(Epi_sce, reduction = "umap", group.by='orig.ident',raster=FALSE) 
plot2
plot3 = DimPlot(Epi_sce, reduction = "umap",split.by = "group",label = T,raster=FALSE) 
plot3
plot4 = DimPlot(Epi_sce, reduction = "umap",group.by = "group",shuffle = T,raster=FALSE) 
dev.off()

#save
Epi_sce@meta.data$Origin_group = ifelse(Epi_sce@meta.data$Sample_Origin == "nLung", 
                                        "nLung", 
                                        "tLung")
save(Epi_sce,file = "Epi_sce.Rdata")
load("Epi_sce.Rdata")

#FindAllMarkers#
markers <- FindAllMarkers(object = Epi_sce, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25,
                          min.pct = 0.25)
all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) 
top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

table(Epi_sce$seurat_clusters)

library(Seurat)

Idents(Epi_sce) <- "seurat_clusters"

#subtype
old.clusters <- as.character(Idents(Epi_sce))

epi_subtype4 <- old.clusters 

epi_subtype4[old.clusters %in% c("6")] <- "AT1" 
epi_subtype4[old.clusters %in% c("1","11")] <- "AT2"
epi_subtype4[old.clusters %in% c("8")] <- "Club"      # Club / secretory
epi_subtype4[old.clusters %in% c("7")]     <- "Ciliated"            
epi_subtype4[old.clusters %in% c("9","0","2","3","4","5","10")] <- "Malignant"


Epi_sce$epi_subtype4 <- factor(
  epi_subtype4,
  levels = c("AT1","AT2","Club", 
             "Ciliated", "Malignant")
)

# to epi_subtype4
Idents(Epi_sce) <- "epi_subtype4"

markers.list <- c(
  
  'AGER', 'RTKN2', 'CLIC5',

  "HHIP","CA2",'ETV5' , 'WIF1', 
 
  'DCDC2B','PIFO',"FOXJ1", 'CCDC78',

  "SCGB1A1", "SCGB3A2","SCGB3A1",  "CYP2F1"
  
  )

markers.use <- unique(markers.list)

p_dot_subtype <- DotPlot(
  Epi_sce,
  features = markers.use
) + RotatedAxis()

p_dot_subtype

DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE) 

DimPlot(Epi_sce, pt.size = 1, label = TRUE, label.size = 5, reduction = "umap",shuffle = T)
DimPlot(Epi_sce, pt.size = 1, group.by= "seurat_clusters",label = TRUE, label.size = 5, reduction = "umap",shuffle = T)
DimPlot(Epi_sce, pt.size = 1, group.by= "Sample_Origin",label = TRUE, label.size = 5, reduction = "umap",shuffle = T)

?DimPlot
Idents(Epi_sce) = "MetaCluster"

save(Epi_sce,file = "Epi_annoted.Rdata")


####Umap
custom_colors2=c(
  "#E84D49",
  "#F9E0AB",
  "#6583A8",  
  "#AA83B4", 

  "#AE8867", "#78966D", "#86AFBE", "#A8A588", "#7D8B85", "#9C4A31", "#644839")

p1=DimPlot(Epi_sce, reduction = "umap", group.by = "epi_subtype4",
           label = F, repel = F, shuffle = T ,pt.size = 0.01,label.size = 10) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 3,hjust = 0.03),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
  ) +

  coord_cartesian(clip = "off")+
  scale_color_manual(values = alpha(custom_colors2,0.3)) +  

  scale_fill_manual(values = alpha(custom_colors2,0.3))
p1
ggsave(filename = "epiumap.pdf", width = 6, height = 5 )



p2=DimPlot(Epi_sce, reduction = "umap", group.by = "seurat_clusters",
           label = T, repel = F, shuffle = T ,pt.size = 0.01,label.size = 5) + 
  theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.2, "inches"),type = "closed")) +
  theme(panel.grid = element_blank(),
        axis.title = element_text(face = 3,hjust = 0.03),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)
  ) 
  coord_cartesian(clip = "off")+
  scale_color_manual(values = alpha(custom_colors2,0.3)) +  

  scale_fill_manual(values = alpha(custom_colors2,0.3))
p2
ggsave(filename = "epiumap cluster.pdf", width = 6, height = 5 )

####epicheck cluster
library(dplyr)
cluster_labels <- Epi_sce@meta.data %>%
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
ggsave(plot=p3, filename="epicheck cluster.pdf",width=9, height=9)

#####CNV######
load("Epi_annoted.Rdata")

DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE,label.size = 8)
DimPlot(Epi_sce, reduction = "umap",label = T,split.by = "group",raster=FALSE)
DimPlot(Epi_sce, reduction = "umap",label = T,raster=FALSE,group.by = "Origin_group")


table(Epi_sce$seurat_clusters,Epi_sce$Origin_group)
table(Epi_sce$Sample)

#test
Epi_sce <- Epi_sce[, sample(1:ncol(Epi_sce),round(ncol(Epi_sce)/5))]
table(Epi_sce$seurat_clusters,Epi_sce$Origin_group)

#
matrix_counts <- as.matrix(GetAssayData(Epi_sce, assay="RNA", layer='counts'))

Epi_sce$newcelltype <- paste0(Epi_sce$Origin_group,"_",Epi_sce$seurat_clusters)
unique(Epi_sce$newcelltype)

library(readr)
# https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
gencode <- read_tsv("hg38_gencode_v27.txt", col_names = c("gene", "chr", "start", "end"))
gencode <- gencode[!duplicated(gencode$gene),]
common_genes <- intersect(gencode$gene, rownames(matrix_counts))

table(Epi_sce$newcelltype)
sort(unique(Epi_sce$newcelltype))
#inferfercnvObject
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = matrix_counts[common_genes,],#counts
                                    annotations_file="./celltype.label1128_2.txt",#celltype
                                    delim="\t",
                                    gene_order_file="./hg38_gencode_v27.txt",#ref
                                    ref_group_names=c("nLung_1"  ,"nLung_11","nLung_8" ),#reference
                                    chr_exclude=c("chrY", "chrM"))



# infercnv analysis
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, 
                             out_dir="infercnv1128_2/", 
                             no_prelim_plot = T,
                             cluster_by_groups=T,
                             denoise=TRUE,
                             HMM=F,
                             min_cells_per_gene = 10,
                             num_threads=10,
                             # num_threads = parallel::detectCores(),  
                             write_expr_matrix = T

)

save(infercnv_obj,file = "infercnv.Rdata")
#infercnv.png
#infercnv.observations.txt
#infercnv.references.txt
#run.final.infercnv_obj


#infercnv.png
library(RColorBrewer)
infercnv::plot_cnv(infercnv_obj, 
                   output_filename = "inferCNV_heatmap1128_2",
                   output_format = "pdf",
                   custom_color_pal =  color.palette(c("#2067AE","white","#B11F2B"))) #改颜色


#####dotplot######
table(Epi_sce$epi_subtype4)

genes <- c("EDAR", "GSTA3","MAOB", "AKAP12", "ALDOA",  "F2","LDHA", "SLC47A1")

library(RColorBrewer)
Idents(Epi_sce) <- Epi_sce$epi_subtype4

dp <- DotPlot(Epi_sce, features = genes)$data
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
  filename = paste0("/epi_subtype4Compare.pdf"),
  plot = p_gene,
  height = 4,
  width = 7,
  dpi = 600,
  device = cairo_pdf
)
