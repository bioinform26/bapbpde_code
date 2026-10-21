rm(list = ls())
setwd("\\ROC")

library(tidyverse)
library(dplyr)
library(pROC)
library(data.table)
data<-fread("TCGA.csv",header = T,check.names = F,data.table = F)
rownames(data)=data$V1
data=data[,2:ncol(data)]

#group####
group=data.frame(ID=colnames(data),
                 group=ifelse(sapply(strsplit(sapply(strsplit(colnames(data),"\\-"), "[", 4),""),"[",1)=="0",
                              "Tumor",
                              "Normal")
)

#########"SLC47A1","AKAP12","LDHA","F2","GABRA2"  ROC#######
genes = c("SLC47A1","AKAP12","LDHA","F2","GABRA2")
fdata<- cbind(group, t(data[genes,]))
head(data)

#---- ROC ----
roc_formula <-paste("group ~",
                    paste(genes, 
                          collapse ="+"))%>% as.formula
roc <- roc(roc_formula, 
           data = fdata,
           percent = F,
           ci= TRUE,
           smoooth = F,
           levels=c("Normal","Tumor"))

### ROC_curves
library(ggplot2)
library(ggsci)
library(ggprism)

roc1<-roc[c("SLC47A1","AKAP12","LDHA","F2","GABRA2")]
names(roc1) <- sapply(names(roc1), function(gene) {
  sprintf("%s (AUC = %.3f)", gene, roc[[gene]]$auc)
})


p <- ggroc(roc1, legacy.axes = TRUE, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linewidth = 0.6) +
  ggsci::scale_color_lancet() +
  theme_prism(border = TRUE) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13, face = "italic", color = "gray30"),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 12, face = "plain"),
    legend.position = c(0.65, 0.17),
    legend.background = element_rect(fill = "NA", color = NA),
    legend.key.width = unit(1.2, "line"),
    legend.key.height = unit(0.6, "line"),
    legend.text = element_text(size = 15),
    text = element_text(family = "sans")
  ) 

print(p)
ggsave("5ROC_curves.pdf", plot = p, width = 5, height = 5)

##########"ALDOA","EDAR","TRPA1","MAOB","SLCO1B3","GSTA3"##########
genes = c("ALDOA","EDAR","TRPA1","MAOB","SLCO1B3","GSTA3")
fdata<- cbind(group, t(data[genes,]))
head(data)

#---- ROC ----
roc_formula <-paste("group ~",
                    paste(genes, 
                          collapse ="+"))%>% as.formula
roc <- roc(roc_formula, 
           data = fdata,
           percent = F,
           ci= TRUE,
           smoooth = F,
           levels=c("Normal","Tumor"))

### ROC_curves
library(ggplot2)
library(ggsci)
library(ggprism)

roc1<-roc[c("ALDOA","EDAR","TRPA1","MAOB","SLCO1B3","GSTA3")]
names(roc1) <- sapply(names(roc1), function(gene) {
  sprintf("%s (AUC = %.3f)", gene, roc[[gene]]$auc)
})


p <- ggroc(roc1, legacy.axes = TRUE, size = 1) +
  geom_abline(slope = 1, intercept = 0, color = "gray60", linewidth = 0.6) +
  ggsci::scale_color_lancet() +
  theme_prism(border = TRUE) +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13, face = "italic", color = "gray30"),
    axis.title = element_text(size = 14, face = "plain"),
    axis.text = element_text(size = 12, face = "plain"),
    legend.position = c(0.65, 0.17),
    legend.background = element_rect(fill = "NA", color = NA),
    legend.key.width = unit(1.2, "line"),
    legend.key.height = unit(0.6, "line"),
    legend.text = element_text(size = 15),
    text = element_text(family = "sans")
  ) 

print(p)
ggsave("6ROC_curves.pdf", plot = p, width = 5, height = 5)