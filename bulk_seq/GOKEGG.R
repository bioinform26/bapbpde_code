rm(list=ls())
setwd("\\KEGGGO")

library(tidyverse)
library(dplyr) 
library(clusterProfiler) 
library(DOSE)
library(ggplot2) 
library(RColorBrewer)
library(org.Hs.eg.db)

DEG <- as.data.frame(read.csv("Intersection_Genes.csv"))

DEG$SYMBOL<-DEG$Gene
genelist <- bitr(DEG$Gene, fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(DEG,genelist,by="SYMBOL")

####GO####
ego <- enrichGO(gene = DEG$ENTREZID,
                OrgDb = org.Hs.eg.db, 
                ont = "all",
                pAdjustMethod = "BH",
                minGSSize = 1,
                pvalueCutoff =0.05, 
                qvalueCutoff =0.05,
                readable = TRUE)

ego_res <- ego@result
ego_res <- ego_res[ego_res$p.adjust < 0.05, ]
write.csv(ego_res,file = "gores.csv")



####KEGG####
kk <- enrichKEGG(gene         = DEG$ENTREZID,
                 organism     = 'hsa',
                 pvalueCutoff = 0.05,
                 qvalueCutoff =0.05,
                 pAdjustMethod= "BH",
                 minGSSize= 5,
                 maxGSSize= 500)
kk<- setReadable(kk,
                 OrgDb= org.Hs.eg.db,
                 keyType= "ENTREZID")

#Rich Factor：
KEGG_diff2<- mutate(kk,
                    RichFactor= Count / as.numeric(sub("/\\d+", "", BgRatio)))

#Fold Enrichment：
KEGG_diff2<- mutate(KEGG_diff2, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
KEGG_diff2@result$RichFactor[1:6]
KEGG_diff2@result$FoldEnrichment[1:6]

#save output：
KEGG_result<- KEGG_diff2@result

KEGG_result <- KEGG_result[KEGG_result$p.adjust < 0.05, ]
write.csv(KEGG_result, file = c('KEGGres.csv'))

