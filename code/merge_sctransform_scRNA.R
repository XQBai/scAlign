## Author: Sue Grimes <sgrimes@stanford.edu> 
## Rerun Seurat pipeline for merged object
library(dplyr)
library(Seurat)
library(viridis)
library(cowplot)
library(future)

dims.reduce = 20
cluster.res = 1

seu_merged <- loadRObj('seurat_merged.Robj')

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=12)

# SCtransform, PCA/clustering, UMAP
aggr <- seu_merged %>% SCTransform(return.only.var.genes=FALSE) %>% RunPCA(algorithm=4, approx=F) 
aggr <- aggr %>% FindNeighbors(dims=1:dims.reduce) %>% FindClusters(resolution=cluster.res) 
aggr <- aggr %>% BuildClusterTree(reorder=T, reorder.numeric=T, dims=1:dims.reduce) %>% RunUMAP(dims=1:dims.reduce) 
saveRDS(aggr, file="seurat_aggr.rds")

p1 <- DimPlot(aggr, reduction="umap", label=TRUE)
p2 <- DimPlot(aggr, reduction="umap", group.by="condition", label=TRUE)
pdf("umap.pdf")
plot(p1)
plot(p2)
dev.off()
gc()

cellnumbers <- table(Idents(aggr), aggr@meta.data$condition)
write.csv(cellnumbers, file="cellnumbers.csv")

# Normalize and scale RNA counts in preparation for find markers step
DefaultAssay(aggr) <- "RNA"

aggr <- NormalizeData(aggr, assay="RNA")
aggr <- FindVariableFeatures(aggr, assay="RNA", selection.method="vst", nfeatures=2000)

# Using default of just variable genes for scaling since on extremely larged merged object using all genes takes huge amount
#   of resources (esp memory) to do this, and can fail after 16+ hrs with memory errors, vs about 30 mins 
#   execution time with just variable genes.  Can use a custom set of features if needed (eg combine 
#   variable genes with a set of standard markers needed for heatmap plot, in case any of those are not in variable gene list)  
aggr <- ScaleData(aggr, assay="RNA", vars.to.regress="nCount_RNA")
  
saveRDS(aggr, file="seurat_aggr.rds")
