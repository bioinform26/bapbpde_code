rm(list = ls())
setwd("/WGCNA")

library(WGCNA)
library(data.table)

WGCNA0=fread("TCGA_bind_GTEx.txt",header = T,data.table = F,check.names = F)
rownames(WGCNA0)<-WGCNA0[,1]
WGCNA0<-WGCNA0[,-1]

WGCNA0 <- t(WGCNA0)
dim(WGCNA0)

# top 50% variable gene
vars_res <- apply(WGCNA0, 2, var)
per_res <- quantile(vars_res, probs = seq(0, 1, 0.25)) 
per_res
upperGene <- WGCNA0[, which(vars_res > per_res[3])] 
dim(upperGene)

WGCNA1 <- data.matrix(upperGene)
nGenes <- ncol(WGCNA1)
nSamples <- nrow(WGCNA1)
head(WGCNA1[1:6,1:6])

####clinical information####
library(tidyverse)
##GTEX
clinical=data.table::fread("GTEX_phenotype.gz")
table(clinical$`_primary_site`)
s=intersect(rownames(WGCNA1),clinical$Sample)
clinical=clinical[match(s,clinical$Sample),]
CTExclin<-clinical
####LUAD
clinical=data.table::fread("TCGA-LUAD.clinical.tsv.gz")
s=intersect(rownames(WGCNA1),clinical$sample)
clinical=clinical[match(s,clinical$sample),]
clinical <-clinical[,c("sample","gender.demographic")]
TCGAclin<-clinical

batch_TCGA <- data.frame(
  ID = TCGAclin$sample,
  Gender = TCGAclin$gender.demographic,
  batch = "TCGA",
  group = ifelse(
    substr(sapply(strsplit(TCGAclin$sample, "-"), `[`, 4), 1, 2) %in% 
      c("01", "02", "03", "04", "05", "06", "07", "08", "09"),
    "Tumor",
    ifelse(
      substr(sapply(strsplit(TCGAclin$sample, "-"), `[`, 4), 1, 2) %in% 
        c("10", "11", "12", "13", "14", "15", "16", "17", "18", "19"),
      "Normal",
      "Unknown"
    )
  )
)
table(batch_TCGA$group)  
batch_GTEx <- data.frame(
  ID = CTExclin$Sample,  
  Gender = CTExclin$`_gender`,
  batch = "GTEx",
  group = "Normal"
)

batch <- rbind(batch_TCGA, batch_GTEx)
table(batch$group) 
write.csv(batch,"TCGA_GTEx_bind_group.csv")
rm(batch_TCGA,batch_GTEx,clinical,TCGAclin,CTExclin,s)

batch<-fread("TCGA_GTEx_bind_group.csv",header = T,data.table = F,check.names = F)
datTraits <- batch[,-1]

colnames(batch)[2]
head(batch[[2]])
head(rownames(WGCNA1))

# check sample 
identical(as.character(datTraits$ID), rownames(WGCNA1))

# check ID 
anyDuplicated(datTraits$ID)
sum(is.na(datTraits$ID))
sum(!rownames(WGCNA1) %in% datTraits$ID)
stopifnot(
  identical(as.character(datTraits$ID), rownames(WGCNA1)),
  !anyDuplicated(datTraits$ID)
)

rownames(datTraits) <- datTraits$ID

library(stringr)
datTraits=as.data.frame(datTraits)

###check
gsg <- goodSamplesGenes(WGCNA1, verbose = 3);
gsg[["allOK"]] 

sampleTree <- hclust(dist(WGCNA1), method = "average")

traitColors <- numbers2colors(as.numeric(factor(datTraits$group)), 
                              colors = rainbow(length(table(datTraits$group))), 
                              signed = FALSE)

plotDendroAndColors(sampleTree, 
                    traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

powers <- c(c(1:10), 
            seq(from = 12, 
                to = 30,
                by = 2))

###power=sft$powerEstimate####
sft <- pickSoftThreshold(as.matrix(WGCNA1),
                         powerVector = powers, 
                         networkType = "signed")

power = sft$fitIndices$Power[min(which(sft$fitIndices$SFT.R.sq > 0.8))]


par(mfrow = c(1,2))  
cex1 = 0.8  

# scale independence and mean connectivity

plot(sft$fitIndices[,1],  
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],  
     xlab = "Soft Threshold (power)",  
     ylab = "Scale Free Topology Model Fit,signed R^2",  
     type = "n",  
     main = paste("Scale independence")) +  
  text(sft$fitIndices[,1],  
       -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels = powers, 
       cex = cex1,
       col = "steelblue") +
  abline(h = 0.80,  
         col = "red")

plot(sft$fitIndices[,1],  
     sft$fitIndices[,5],  
     xlab = "Soft Threshold (power)",  
     ylab = "Mean Connectivity",  
     type="n",  
     main = paste("Mean connectivity")) +  
  text(sft$fitIndices[,1],  
       sft$fitIndices[,5],
       labels = powers, 
       cex = cex1, 
       col = "steelblue")

enableWGCNAThreads()  

if(T){ 
  net = blockwiseModules(
    as.matrix(WGCNA1),  
    power = 2,  #sft$powerEstimate,  
    maxBlockSize = nGenes, 
    networkType = "signed", 
    #networkType = "unsigned",
    TOMType = "unsigned",  
    minModuleSize = 25,  
    deepSplit = 3,           
    reassignThreshold = 0,  
    mergeCutHeight = 0.20,  
    numericLabels = TRUE,  
    pamRespectsDendro = FALSE,  
    saveTOMs = F,  
    verbose = 3  
  )
  table(net$colors)  
}
if(T){  
  moduleColors=labels2colors(net$colors)
  table(moduleColors)  
    plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

library(WGCNA) 
library(forcats)  
datTraits$group <- factor(datTraits$group)  
levels(datTraits$group) 
if(T){ 
  nGenes = ncol(WGCNA1) 
  nSamples = nrow(WGCNA1) 
  design <- model.matrix(~0+datTraits$group)  
  colnames(design)= levels(datTraits$group) 
  MES0 <- moduleEigengenes(WGCNA1,moduleColors)$eigengenes  
  MEs = orderMEs(MES0) 
  moduleTraitCor <- cor(MEs,design,use = "p") 
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  
  dim(textMatrix)=dim(moduleTraitCor)  
  par(mar=c(6, 8.5, 3, 3))  
  labeledHeatmap(Matrix = moduleTraitCor,  
                 xLabels = colnames(design), 
                 yLabels = names(MEs),  
                 ySymbols = names(MEs),  
                 colorLabels = FALSE,  
                 colors = blueWhiteRed(50),  
                 textMatrix = textMatrix,  
                 setStdMargins = FALSE,  
                 cex.text = 0.5,  
                 zlim = c(-1,1),  
                 main = paste("Module-trait relationships"))  
}

module="turquoise"  
probes = colnames(WGCNA1)  
inModule = (moduleColors==module)  
modProbes=probes[inModule]  
write.csv(modProbes,"power2turquoise_gene.csv")

#####power12############

rm(net)
if(T){ 
  net = blockwiseModules(
    as.matrix(WGCNA1), 
    power = 12,  
    maxBlockSize = nGenes,
    networkType = "signed", 
    #networkType = "unsigned",
    TOMType = "unsigned",  
    minModuleSize = 25,  
    deepSplit = 3,           
    reassignThreshold = 0,  
    mergeCutHeight = 0.20,  
    numericLabels = TRUE,  
    pamRespectsDendro = FALSE,  
    saveTOMs = F,  
    verbose = 3  
  )
  table(net$colors)  
}
if(T){ 
  moduleColors=labels2colors(net$colors)
  table(moduleColors) 
  
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

library(WGCNA)  
library(forcats)  
datTraits$group <- factor(datTraits$group)  
levels(datTraits$group)  
if(T){ 
  nGenes = ncol(WGCNA1)  
  nSamples = nrow(WGCNA1)  
  design <- model.matrix(~0+datTraits$group)  
  colnames(design)= levels(datTraits$group)  
  MES0 <- moduleEigengenes(WGCNA1,moduleColors)$eigengenes  
  MEs = orderMEs(MES0)  
  moduleTraitCor <- cor(MEs,design,use = "p")  
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  
  dim(textMatrix)=dim(moduleTraitCor) 
  par(mar=c(6, 8.5, 3, 3)) 
  labeledHeatmap(Matrix = moduleTraitCor,  
                 xLabels = colnames(design),  
                 yLabels = names(MEs),  
                 ySymbols = names(MEs),  
                 colorLabels = FALSE,  
                 colors = blueWhiteRed(50),  
                 textMatrix = textMatrix,  
                 setStdMargins = FALSE,  
                 cex.text = 0.5,  
                 zlim = c(-1,1),  
                 main = paste("Module-trait relationships12"))  
}

module="turquoise"  
probes = colnames(WGCNA1)  
inModule = (moduleColors==module)  
modProbes=probes[inModule]  
write.csv(modProbes,"power12turquoise_gene.csv")

#####power2############
rm(net)
if(T){ 
  net = blockwiseModules(
    as.matrix(WGCNA1),  
    power = 14,
    maxBlockSize = nGenes, 
    networkType = "signed", 
    #networkType = "unsigned",
    TOMType = "unsigned",  
    minModuleSize = 25,  
    deepSplit = 3,           
    reassignThreshold = 0,  
    mergeCutHeight = 0.20,  
    numericLabels = TRUE,  
    pamRespectsDendro = FALSE,  
    saveTOMs = F,  
    verbose = 3 
  )
  table(net$colors)  
}
if(T){  
  moduleColors=labels2colors(net$colors)
  table(moduleColors)  
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
}

library(WGCNA)  
library(forcats)  
datTraits$group <- factor(datTraits$group)  
levels(datTraits$group)  
if(T){ 
  nGenes = ncol(WGCNA1)  
  nSamples = nrow(WGCNA1)  
  design <- model.matrix(~0+datTraits$group) 
  colnames(design)= levels(datTraits$group)  
  MES0 <- moduleEigengenes(WGCNA1,moduleColors)$eigengenes  
  MEs = orderMEs(MES0)  
  moduleTraitCor <- cor(MEs,design,use = "p")  
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor,nSamples)  
  textMatrix = paste(signif(moduleTraitCor,2),"\n(",
                     signif(moduleTraitPvalue,1),")",sep = "")  
  dim(textMatrix)=dim(moduleTraitCor)  
  par(mar=c(6, 8.5, 3, 3))  
  labeledHeatmap(Matrix = moduleTraitCor,  
                 xLabels = colnames(design),  
                 yLabels = names(MEs),  
                 ySymbols = names(MEs),  
                 colorLabels = FALSE,  
                 colors = blueWhiteRed(50),  
                 textMatrix = textMatrix,  
                 setStdMargins = FALSE,  
                 cex.text = 0.5,  
                 zlim = c(-1,1),  
                 main = paste("Module-trait relationships14"))  
}

module="turquoise"  
probes = colnames(WGCNA1)  
inModule = (moduleColors==module)  
modProbes=probes[inModule]  
write.csv(modProbes,"power14turquoise_gene.csv")


