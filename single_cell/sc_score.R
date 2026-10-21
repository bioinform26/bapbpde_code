rm(list=ls())
setwd("\\ScRNA\\score")

load("/GSE131907_res0.5.RData")

library(Seurat)
library(gplots)
library(ggplot2) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(GSVA) 
library(GSEABase)
library(data.table)
score1 <- list(c("SLC47A1", "AKAP12", "LDHA", "F2", "ALDOA", "TRPA1","SLCO1B3",
                 "GABRA2", "EDAR", "MAOB", "GSTA3"))

score <- AddModuleScore(obj,
                        features = score1,
                        ctrl = 100,
                        name = "score")
save(score,file = "GSE131907_res0.5_score.RData")

library(dplyr)
library(ggplot2)
mydata<- FetchData(score,vars = c("umap_1","umap_2","score1","group_manual","Sample_Origin"))

## group
mydata <- mydata %>% 
  mutate(
    Origin_group = ifelse(Sample_Origin == "nLung", 
                          "nLung", 
                          "tLung")
  )
mydata$Origin_group <- factor(mydata$Origin_group,
                              levels = c("nLung", "tLung"))

range(mydata$score1)

####UMAP visualization of target gene scores####
library(ggplot2)

a <- ggplot(mydata, aes(x = umap_1, y = umap_2, colour = score1)) +
  geom_point(size = 0.5) +

  scale_color_gradientn(
    values  = seq(0, 1, 0.2),
    colours = c('#333366', "#6666FF", '#CC3333', '#FFCC33'),
    name    = "Score"
  ) +

  facet_wrap(~ Origin_group, nrow = 1) +

  coord_cartesian(clip = "off") +
  theme(
    panel.grid       = element_blank(),         
    panel.background = element_blank(),         
    panel.border     = element_blank(),         
    axis.title       = element_text(face = 3, hjust = 0.03),
    axis.text        = element_blank(),         
    axis.ticks       = element_blank(),
    
    legend.text      = element_text(size = 14),
    legend.title     = element_text(size = 16),
    
    strip.background = element_blank(),         
    strip.text       = element_text(size = 16, face = "bold")
  )

a
####Score distribution across cell subtypes and Sample_Origin#####
library(ggplot2)
library(Seurat)

score$Sample_Origin = ifelse(score$Sample_Origin == "nLung", 
                             "nLung", 
                             "tLung")
table(score$Sample_Origin)

library(Seurat)
library(ggplot2)
library(dplyr)
library(rstatix)   

df <- score@meta.data %>%
  dplyr::select(group_manual, score1, Sample_Origin)

# group
df$group_manual <- factor(df$group_manual)
df$Sample_Origin <- factor(df$Sample_Origin)

p <- ggplot(df, aes(x = group_manual, y = score1, fill = Sample_Origin
)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA, width = 0.7) +
  #  geom_jitter(
  #    aes(color = Sample_Origin),
  #    size = 0,
  #    width = 0.2,
  #    alpha = 0.5
  #  ) +
  scale_fill_brewer(palette = "Set1") +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 16) +
  theme(
    axis.text = element_text(color = "black"),          
    axis.title = element_text(color = "black"),         
    plot.title = element_text(color = "black", hjust=0.5), 
    legend.text = element_text(color = "black"),        
    legend.title = element_text(color = "black"),       
    
    panel.grid = element_blank(),                      
    panel.border = element_blank(),                     
    axis.line = element_line(color = "black"),          
    
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  #xlab("Cell Subtype (group_manual)") +
  ylab("Score1")
#ggtitle("Score distribution across cell subtypes and Sample_Origin")
##Wilcoxon 
stat.test <- df %>%
  group_by(group_manual) %>%
  wilcox_test(score1 ~ Sample_Origin) %>%           
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")
## p.signif
stat.test <- stat.test %>%
  mutate(
    p.signif = case_when(
      p.adj < 0.001 ~ "***",
      p.adj < 0.01  ~ "**",
      p.adj < 0.05  ~ "*",
      TRUE ~ "ns"
    )
  )
y.pos <- df %>%
  group_by(group_manual) %>%
  summarise(y.position = max(score1, na.rm = TRUE) * 1.05)
stat.test <- stat.test %>%
  left_join(y.pos, by = "group_manual")
p_sig <- p +
  geom_text(
    data = stat.test,
    aes(x = group_manual, y = y.position, label = p.signif),
    inherit.aes = FALSE,      
    color = "black",
    size = 4,
    vjust = 0
  )
print(p_sig)
#########featureplot######
marker = c("EDAR","GABRA2",  "GSTA3" ,"MAOB", "AKAP12", "ALDOA","F2","LDHA",  "SLCO1B3",
           "SLC47A1","TRPA1" )
p3<-FeaturePlot(
  obj,
  features   = marker,
  ncol       = 4,
  pt.size    = 0.1,
  order      = TRUE                           
)
ggsave("featurePlot.pdf",
       plot = p3,
       width = 12,
       height = 7,
       #units = "in",
       dpi = 600 )

dev.off()