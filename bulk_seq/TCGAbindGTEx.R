rm(list=ls())
setwd("\\predada")

library(sva)
library(FactoMineR)
library(factoextra)
library(limma)
library(data.table)
library(parallel)
library(cowplot)
library(WGCNA)

randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(10,randomColor())
setDTthreads(threads =detectCores())
######TCGA####
TCGAdata=fread("TCGA.txt",header = T,check.names = F,data.table = F)
rownames(TCGAdata)=TCGAdata$V1
TCGAdata=TCGAdata[,2:ncol(TCGAdata)]
TCGAdata=log2(TCGAdata+0.001)

group=data.frame(ID=colnames(TCGAdata),
                 group=ifelse(sapply(strsplit(sapply(strsplit(colnames(TCGAdata),"\\-"), "[", 4),""),"[",1)=="0","Tumor","Normal"))
rownames(group)=group$ID

dat0=as.data.frame(t(TCGAdata))
dat.pca0 <- PCA(dat0, graph = FALSE)
pca_plot0 <- fviz_pca_ind(dat.pca0,geom.ind = "point",
                          col.ind = group$group,
                          palette = randomColors[1:2],
                          addEllipses = TRUE, 
                          legend.title = "Groups")
ggsave(plot = pca_plot0,filename ="PCA_TCGA.pdf")

TCGAdatahclust=as.matrix(t(TCGAdata))
sampleTree = hclust(dist(TCGAdatahclust), method = "average")
pdf(file = "beforeclustering_TCGA.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
dev.off()

hlevel=400

#Draw the dividing line

pdf(file = "dividing line_TCGA.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
abline(h = hlevel, col = "red")#
dev.off()


clust = cutreeStatic(sampleTree, cutHeight = hlevel, minSize = 10)

table(clust)
keepSamples = (clust==1)
TCGAdata_in=TCGAdata[,keepSamples]
TCGAdata_in=normalizeBetweenArrays(TCGAdata_in)
group_in=group[keepSamples,]
dat_in=as.data.frame(t(TCGAdata_in))
dat.pca_in <- PCA(dat_in, graph = FALSE)
pca_plot_in <- fviz_pca_ind(dat.pca_in,geom.ind = "point",
                            col.ind = group_in$group,
                            palette = randomColors[3:4],
                            addEllipses = TRUE, 
                            legend.title = "Groups")
ggsave(plot = pca_plot_in,filename ="afterPCA_TCGA.pdf")
fwrite(TCGAdata_in,"TCGA_out.txt",quote = F,sep = "\t",row.names = T)

#GTEx####
GTEx<-fread("GTEx.txt",,header = T,check.names = F,data.table = F)
rownames(GTEx)=GTEx$V1
GTEx<-GTEx[,2:ncol(GTEx)]
#PCA
dat1=as.data.frame(t(GTEx))
dat.pca1 <- PCA(dat1, graph = FALSE)
pca_plot1 <- fviz_pca_ind(dat.pca1,geom.ind = "point",
                          addEllipses = TRUE)
ggsave(plot = pca_plot1,filename ="beforePCA_GTEx.pdf")
#Sample clustering
GTExhclust<-as.matrix(t(GTEx))
sampleTree_GTEx = hclust(dist(GTExhclust), method = "average")
pdf(file = "Sample clustering_GTEx.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree_GTEx, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
dev.off()

hlevel_GTEx=800
#Draw the dividing line

pdf(file = "samplecut_GTEx.pdf", width = 10, height = 10)
par(cex = 0.5)
par(mar = c(0,4,2,0))
plot(sampleTree_GTEx, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
abline(h = hlevel_GTEx, col = "red")
dev.off()

clust_GTEx = cutreeStatic(sampleTree_GTEx, cutHeight = hlevel_GTEx, minSize = 10)

table(clust_GTEx)

keepSamples_GTEx = (clust_GTEx==1)

GTEx_in=GTEx[,keepSamples_GTEx]

GTEx_in=normalizeBetweenArrays(GTEx_in)

dat_in_GTEx=as.data.frame(t(GTEx_in))
dat.pca_in_GTEx <- PCA(dat_in_GTEx, graph = FALSE)
pca_plot_in_GTEx <- fviz_pca_ind(dat.pca_in_GTEx,geom.ind = "point",
                                 addEllipses = TRUE)
ggsave(plot = pca_plot_in_GTEx,filename ="PCA_GTEx_in.pdf")
GTEx_in<-as.data.frame(GTEx_in)

fwrite(GTEx_in,"GTEx_out.txt",quote = F,sep = "\t",row.names = T)

#######TCGA_bind_GTEx#######
rm(list=ls())
setwd("E:\\1predata")

library(sva)
library(FactoMineR)
library(factoextra)
library(limma)
library(data.table)
library(parallel)
randomColor <- function() {
  paste0("#",paste0(sample(c(0:9, letters[1:6]), 6, replace = TRUE),collapse = ""))
}
randomColors <- replicate(10,randomColor())
setDTthreads(threads =detectCores())

TCGA_data=fread("TCGA_out.txt",header = T,data.table = F,check.names = F)
rownames(TCGA_data)=TCGA_data$V1
TCGA_data=TCGA_data[,2:ncol(TCGA_data)]

GTEx_data=fread("GTEx_out.txt",header = T,check.names = F,data.table = F)
rownames(GTEx_data)=GTEx_data$V1
GTEx_data=GTEx_data[,2:ncol(GTEx_data)]

#data merging
rownames(TCGA_data)=make.names(rownames(TCGA_data))
rownames(GTEx_data)=make.names(rownames(GTEx_data))
sameid=intersect(rownames(TCGA_data),rownames(GTEx_data))
TCGA_data=TCGA_data[sameid,]
GTEx_data=GTEx_data[sameid,]
bind_data=cbind(TCGA_data,GTEx_data)
data2=apply(bind_data,2,as.numeric)
rownames(data2)=rownames(bind_data)

#batch removal
batch_TCGA=data.frame(ID=colnames(TCGA_data),batch="TCGA",group=ifelse(sapply(strsplit(sapply(strsplit(colnames(TCGA_data),"\\-"), "[", 4),""),"[",1)=="0","Tumor","Normal"))
batch_GTEx=data.frame(ID=colnames(GTEx_data),batch="GTEx",group="Normal")
batch=rbind(batch_TCGA,batch_GTEx)
batchs=batch$batch
groups=batch$group


pdf(file = "beforebatchremove_TCGAGTEx.pdf")
boxplot(data2,outline=F, notch=F,las=2,xaxt="n")
dev.off()

dat_batchs=as.data.frame(t(data2))
dat_batchs <- PCA(dat_batchs, graph = FALSE)
pca_batchs <- fviz_pca_ind(dat_batchs,geom.ind = "point",
                           col.ind = batchs,
                           palette = randomColors[1:2],
                           addEllipses = TRUE, 
                           legend.title = "Batchs")
ggsave(plot = pca_batchs,filename ="beforebatchremove_TCGAGTEx_PCA.pdf")

dat_groups=as.data.frame(t(data2))
dat_groups <- PCA(dat_groups, graph = FALSE)
pca_groups <- fviz_pca_ind(dat_groups,geom.ind = "point",
                           col.ind = groups,
                           palette = randomColors[1:2],
                           addEllipses = TRUE, 
                           legend.title = "Groups")
ggsave(plot = pca_groups,filename ="beforebatchremove_TCGAGTEx_groupPCA.pdf")


mod=model.matrix(~group,data = batch)
bind_combat=ComBat(dat=data2,batch = batchs,mod=mod,par.prior = T)
bind_combat<-normalizeBetweenArrays(bind_combat)

#check after batch removal

pdf(file = "afterbatchremove_TCGAGTEx.pdf")
boxplot(bind_combat,outline=F, notch=F,las=2,xaxt="n")
dev.off()

dat1=as.data.frame(t(bind_combat))
dat_batchs_combat <- PCA(dat1, graph = FALSE)
pca_batchs_combat <- fviz_pca_ind(dat_batchs_combat,geom.ind = "point",
                                  col.ind = batchs,
                                  palette = randomColors[5:6],
                                  addEllipses = TRUE, 
                                  legend.title = "Batchs")
ggsave(plot =pca_batchs_combat,filename ="afterbatchremove_TCGAGTEx_PCA.pdf")
dat_groups_combat <- PCA(dat1, graph = FALSE)
pca_groups_combat <- fviz_pca_ind(dat_groups_combat,geom.ind = "point",
                                  col.ind = batch$group,
                                  palette = randomColors[7:8],
                                  addEllipses = TRUE, 
                                  legend.title = "Groups")
ggsave(plot = pca_groups_combat,filename ="afterbatchremove_TCGAGTEx_groupPCA.pdf")

bind_combat<-as.data.frame(bind_combat)
fwrite(bind_combat,"TCGA_bind_GTEx.txt",quote = F,sep = "\t",row.names = T)
