rm(list=ls())
setwd("/ML")

library(dplyr)
depens<-c('GSEABase', 'GSVA', 'cancerclass', 'mixOmics', 'sparrow', 'sva' , 'ComplexHeatmap' )
for(i in 1:length(depens)){
  depen<-depens[i]
  if (!requireNamespace(depen, quietly = TRUE))  BiocManager::install(depen,update = FALSE)
}
install.packages("usethis")
library(devtools)
install_github("binderh/CoxBoost")
if (!requireNamespace("CoxBoost", quietly = TRUE))
  devtools::install_github("binderh/CoxBoost")

if (!requireNamespace("fastAdaboost", quietly = TRUE))
  devtools::install_github("souravc83/fastAdaboost")

if (!requireNamespace("Mime", quietly = TRUE))
  devtools::install_github("l-magnificence/Mime")
library(Mime1)
library(data.table)

#####ML######
#load data
TCGA <- fread("TCGA_OS_exp.csv")
GSE30219 <- fread("GSE30219_OS_exp.csv")
GSE31210 <- fread("GSE31210_OS_exp.csv")
GSE72094 <- fread("GSE72094_OS_exp.csv")
gene<-read.csv("Intersection_Genes.csv")
gene <- gene$Gene

col<-c("ID","OS.time","OS",gene)
TCGA <- TCGA[,..col]  
GSE30219 <- GSE30219[,..col]
GSE31210 <- GSE31210[,..col]
GSE72094 <- GSE72094[,..col]

list_train_vali_Data <- list(TCGA=TCGA, 
                             GSE30219=GSE30219, GSE31210=GSE31210, 
                             GSE72094=GSE72094)
rm(TCGA, GSE30219, GSE31210,GSE72094)

# check data
sum(list_train_vali_Data$TCGA$OS.time < 30)  #11
sum(list_train_vali_Data$GSE30219$OS.time < 30)#2
sum(list_train_vali_Data$GSE31210$OS.time < 30)
sum(list_train_vali_Data$GSE72094$OS.time < 30)#12

list_train_vali_Data$TCGA <- subset(list_train_vali_Data$TCGA, OS.time > 30)
list_train_vali_Data$GSE30219 <- subset(list_train_vali_Data$GSE30219, OS.time > 30)
list_train_vali_Data$GSE31210 <- subset(list_train_vali_Data$GSE31210, OS.time > 30)
list_train_vali_Data$GSE72094 <- subset(list_train_vali_Data$GSE72094, OS.time > 30)


rm(res)
res <- ML.Dev.Prog.Sig(train_data = list_train_vali_Data$TCGA,
                       list_train_vali_Data = list_train_vali_Data,
                       unicox.filter.for.candi = T,
                       unicox_p_cutoff = 0.05,
                       candidate_genes = gene,
                       mode = 'all',nodesize =5,seed = 1314 )
save(res, file = "ml.RData")

########drawheatmap######
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)

df <- res$Cindex.res

rm(res)
tail(df)
table(df$ID)

df_wide <- df %>%
  mutate(Cindex = as.numeric(Cindex)) %>%        # ensure numeric
  pivot_wider(names_from = ID, values_from = Cindex) # wide format

df_wide <- as.data.frame(df_wide)


rownames(df_wide) <- df_wide$Model
df_wide <- df_wide[, -1]  


head(rownames(df_wide))  

Cindex_mat <- as.matrix(df_wide)

# Compute each model's average C-index
df_wide$avg_Cindex <- rowMeans(Cindex_mat
                               , na.rm = TRUE)
library(dplyr)
df_sorted <- df_wide %>% arrange(desc(avg_Cindex))

avg_Cindex <- round(rowMeans(df_sorted, na.rm = TRUE), 3) 
Cindex_mat <- as.matrix(df_sorted)

#Define color palettes

CohortCol <- c("TCGA" = "#4195C1", 
               "GSE30219" = "#e64b35",
               "GSE31210" = "#f39c12",
               "GSE72094" = "#2ecc71")

SimpleHeatmap <- function(Cindex_mat = NULL, 
                          avg_Cindex = NULL, 
                          CohortCol = NULL, 
                          barCol = NULL,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"), 
                          cellwidth = 1, cellheight = 0.5, 
                          cluster_columns, cluster_rows){
  col_ha = columnAnnotation("Cohort" = colnames(Cindex_mat),
                            col = list("Cohort" = CohortCol),
                            show_annotation_name = F)
  
  row_ha = rowAnnotation(bar = anno_barplot(avg_Cindex, bar_width = 0.8,
                                            border = FALSE,
                                            gp = gpar(fill = barCol, col = NA),
                                            add_numbers = TRUE,  
                                            labels = function(avg_Cindex) sprintf("%.3f", avg_Cindex),number_format = "%.3f",  # 直接指定格式  # 通过labels参数格式化
                                            
                                            numbers_offset = unit(-10, "mm"),
                                            axis_param = list("labels_rot" = 0),
                                            numbers_gp = gpar(fontsize = 9, col = "white"),
                                            width = unit(3, "cm")),
                         show_annotation_name = F)
  
  Heatmap(as.matrix(Cindex_mat), name = "C-index",
          right_annotation = row_ha, 
          top_annotation = col_ha,
          col = colorRamp2(c(0.4, 0.6, 0.95), col),
          
          rect_gp = gpar(col = "black", lwd = 1), 
          cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
          show_column_names = FALSE, 
          show_row_names = TRUE,
          row_names_side = "left",
          width = unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
          height = unit(cellheight * nrow(Cindex_mat), "cm"),
          column_split = factor(colnames(Cindex_mat),
                                levels = colnames(Cindex_mat)), 
          column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, col) { 
            grid.text(label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
                      x, y, gp = gpar(fontsize = 10))
          }
  )
}


p<-SimpleHeatmap(
  Cindex_mat = Cindex_mat,
  avg_Cindex = avg_Cindex,
  CohortCol = c("TCGA" = "#4195C1", 
                "GSE30219" = "#e64b35",
                "GSE31210" = "#f39c12",
                "GSE72094" = "#2ecc71"),
  barCol = rep("#e64b35",length(avg_Cindex)),
  cluster_columns = FALSE, 
  cluster_rows = FALSE
)

tiff("groupCindex_Heatmap.tiff", width = 12, height = 30, units = "in", res = 600, compression = "lzw")
draw(p, merge_legend = TRUE)
dev.off()
library(ggplot2)
library(ggpubr)     
library(Cairo) 

Cairo::CairoPDF(
  "groupCindex_Heatmap.pdf",
  width = 12,         
  height = 30,         
  dpi = 600, 
  units = "in",
  family = "Arial"   
)
draw(p)               
dev.off()