rm(list=ls())
setwd("\\2deg")

library(data.table)
library(ggrepel)
library(tidyverse)
library(ggsignif) 
library(limma)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(beepr)
library(gplots)
library(parallel)
library(future.apply)
library(pheatmap)
setDTthreads(threads =detectCores())

data<-fread("TCGA_bind_GTEx.txt",header = T,check.names = F,data.table = F)
rownames(data)=data$V1
data=data[,2:ncol(data)]

group=data.frame(ID=colnames(data),
                 group=ifelse(sapply(strsplit(colnames(data),"\\-"), "[", 1)=="GTEx",
                              "Normal",
                              ifelse(sapply(strsplit(sapply(strsplit(colnames(data),"\\-"), "[", 4),""),"[",1)=="0",
                                     "Tumor",
                                     "Normal")
                 )
)
rownames(group)=group$ID
design<-data.frame("con"=ifelse(group$group=="Normal",1,0),
                   "treat"=ifelse(group$group=="Tumor",1,0))
#limma differential analysis#
df.fit <- lmFit(data,design)
df.matrix<- makeContrasts(treat - con ,levels=design)
fit<- contrasts.fit(df.fit,df.matrix)

fit2 <- eBayes(fit)

allDEG1 = topTable(fit2,coef=1,n=Inf,adjust="BH") 
allDEG1 = na.omit(allDEG1)

limma_diff = allDEG1[(allDEG1$adj.P.Val < 0.05 & 
                        (allDEG1$logFC>=1 | allDEG1$logFC<=(-1))),]
limma_diff = limma_diff[order(limma_diff$logFC),]
all_limma=allDEG1
all_limma=all_limma[order(all_limma$adj.P.Val),]
fwrite(limma_diff,"TCGA_GTEx_limma_diff.txt",quote = F,sep = "\t",row.names = T)
fwrite(all_limma,"TCGA_GTEx_all_limma.txt",quote = F,sep = "\t",row.names = T)
limma_diff_exp=data[rownames(limma_diff),]
fwrite(limma_diff_exp,"TCGA_GTEx_limma_diff_exp.txt",quote = F,sep = "\t",row.names = T)

#volcano plot
all_limma$logp<--log10(all_limma$adj.P.Val) 
all_limma$Group<-ifelse((all_limma$logFC>=1 & all_limma$adj.P.Val<0.05),"Up",ifelse((all_limma$logFC<=-1 & all_limma$adj.P.Val<0.05),"Down","Not"))
plot<-ggscatter(all_limma,
                x="logFC",
                y="logp",
                color="Group", 
                palette=c("#2f5688","#BBBBBB","#CC0000"),
                size=0.25, 
                xlab="log2FoldChange", 
                ylab="-log10(Adjust P-value)",)+
  theme_base()+ 
  geom_hline(yintercept=-log10(0.05),linetype="dashed")+ 
  geom_vline(xintercept=c(-1,1),linetype="dashed")
pdf("volcano plot.pdf")
print(plot)
dev.off()


#heatmap
diffgroup0=data.frame(ID=colnames(limma_diff_exp),
                      group=ifelse(sapply(strsplit(colnames(limma_diff_exp),"\\-"), "[", 1)=="GTEx",
                                   "Normal",
                                   ifelse(sapply(strsplit(sapply(strsplit(colnames(limma_diff_exp),"\\-"), "[", 4),""),"[",1)=="0",
                                          "Tumor",
                                          "Normal")
                      )
)
rownames(diffgroup0)<-diffgroup0$ID
diffTumor<-limma_diff_exp[,rownames(diffgroup0[diffgroup0$group=="Tumor",])]
diffNormal<-limma_diff_exp[,rownames(diffgroup0[diffgroup0$group=="Normal",])]
diffexp<-cbind(diffNormal,diffTumor)
diffgroup=data.frame(ID=colnames(diffexp),
                     group=ifelse(sapply(strsplit(colnames(diffexp),"\\-"), "[", 1)=="GTEx",
                                  "Normal",
                                  ifelse(sapply(strsplit(sapply(strsplit(colnames(diffexp),"\\-"), "[", 4),""),"[",1)=="0",
                                         "Tumor",
                                         "Normal")
                     )
)
diffgroupmap<-data.frame(group=diffgroup$group)
rownames(diffgroupmap)<-diffgroup$ID
rownames(diffgroup)=diffgroup$ID
diffexp<-diffexp[1:50,]
pheatmaps<-pheatmap(diffexp, color=colorRampPalette(c("#3300CC","#3399FF","white","#FF3333","#CC0000"))(50),border_color=NA,
                    labels_row=NULL,show_rownames=T,show_colnames=F,annotation_col =diffgroupmap,scale = "none",cluster_cols=F)
pdf("heatmap.pdf",width = 15,height = 15)
print(pheatmaps)
dev.off()
