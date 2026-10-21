rm(list=ls())
setwd("\\venn")

library(data.table)
library(ggvenn)
library(ggplot2)
library(eulerr)
library(VennDiagram)

#### venn DEG WGCNA####
file1 <- "\\turquoise_gene.csv"
file2 <- "\\TCGA_GTEx_limma_diff.txt"

WGCNAall <- fread(file1,header = T)
limmall <- fread(file2)

DEG <- unique(limmall[["V1"]])
WGCNA <- unique(WGCNAall[["gene"]])

gene_list <- list(DEG = DEG,
                 WGCNA = WGCNA)
venn.diagram(x=gene_list,
             scaled = F,
             alpha= 0.5, #
             lwd=1,lty=1,col=c('pink','skyblue'),
             label.col ='black',
             cex = 2,
             fontface = "bold",
             fill=c('pink','skyblue'),
             category.names = c("DEG","WGCNA") ,
             cat.dist = 0.02,
             cat.pos = -180, 
             cat.cex = 2, 
             cat.fontface = "bold",  
             cat.col='black' ,  
             cat.default.pos = "outer",  
             output=TRUE,
             filename='venn DEG WGCNA.tiff',
             imagetype="tiff",  
             resolution = 600,  
             compression = "lzw"
)

#### venn bapBPde luad####
file1 <- "\\bapbpde.csv"
file2 <- "\\LUAD.txt"

bapbpdeall <- fread(file1,header = T)
LUADall <- fread(file2)

bapbpde <- unique(bapbpdeall[["Gene"]])
LUAD <- unique(LUADall[["x"]])

gene_list <-list(bapbpde = bapbpde,
                 LUAD = LUAD)
venn.diagram(x=gene_list,
             scaled = F,
             alpha= 0.5, #
             lwd=1,lty=1,col=c('pink','skyblue'),
             label.col ='black',
             cex = 2,
             fontface = "bold",
             fill=c('pink','skyblue'),
             category.names = c("bapbpde","LUAD") ,
             cat.dist = 0.02,
             cat.pos = -180, 
             cat.cex = 2, 
             cat.fontface = "bold",  
             cat.col='black' ,  
             cat.default.pos = "outer",  
             output=TRUE,
             filename='venn bapBPde luad.tiff',
             imagetype="tiff",  
             resolution = 600,  
             compression = "lzw"
)

####Venn power####
read_genes <- function(path) {
  x <- trimws(read.csv(path, stringsAsFactors=FALSE)$x)
  unique(x[!is.na(x) & nzchar(x)])
}
gene_list <- list(
  power12=read_genes('power12turquoise_genes.csv'),
  power14=read_genes('power14turquoise_genes.csv'),
  power2=read_genes('power2turquoise_gene.csv')
)
all_genes <- sort(unique(unlist(gene_list)))
membership <- sapply(gene_list, function(s) all_genes %in% s)
patterns <- apply(membership, 1, paste0, collapse='')
print(sapply(gene_list, length))
print(table(patterns))
write.csv(data.frame(Gene=all_genes, membership), file.path(out, 'gene_membership.csv'), row.names=FALSE)
v <- venn.diagram(x=gene_list, filename=NULL, scaled=TRUE,
                  alpha=0.5, lwd=1, lty=1, col=c('pink','skyblue','#90C99A'),
                  fill=c('pink','skyblue','#90C99A'), label.col='black',
                  cex=2, fontface='bold', fontfamily='sans',
                  category.names=c('power12','power14','power2'),
                  cat.dist=c(0.055,0.055,0.045), cat.pos=c(-25,25,180),
                  cat.cex=1.6, cat.fontface='bold', cat.fontfamily='sans',
                  cat.col='black', cat.default.pos='outer', margin=0.14,
                  disable.logging=TRUE)
tiff(file.path(out,'venn_power.tiff'), width=8, height=8, units='in', res=600, compression='lzw')
grid.draw(v)
dev.off()
png(file.path(out,'venn_power.png'), width=1600, height=1600, res=200)
grid.draw(v)
dev.off()