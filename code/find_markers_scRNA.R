# File: find_markers.R
# Author: Sue Grimes
# Desc: Find markers on merged Seurat objects
#
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Parse input parameters                                                      #                                                 
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
suppressPackageStartupMessages(require(optparse))
# manual: http://cran.r-project.org/web/packages/optparse/optparse.pdf

option_list = list(
  make_option(c("-r", "--seurat"), action="store", default=NA, type='character',
              help="Seurat R object")
)
opt = parse_args(OptionParser(option_list=option_list))

if ( is.null(opt$seurat) ) {
  print(paste0("Required parameter: seurat is missing"))
  quit(status=10)
} else {
  seu_fn = opt$seurat
}

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Load packages, and custom helper functions                                  #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(dplyr)
library(Seurat)
library(viridis)
library(cowplot)
library(future)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers= 24) 

#loads merged object and finds markers
seu <- readRDS(seu_fn)

# add rescaling of data if needed (eg if merge/SCT step only used variable genes)
#seu <- ScaleData(seu, assay="RNA", features=rownames(seu), vars.to.regress="nCount_RNA") 
#saveRDS(seu, file=paste0('mrg_', opt$cell_lineage, '.sct.rds'))

markers_RNA <- FindAllMarkers(seu, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers_RNA, file = "markers_RNA.csv")

top20 <- markers_RNA %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
write.csv(top20, file = "top20_markers_RNA.csv")

top5 <- markers_RNA %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

pdf("top5_cluster.pdf")
DoHeatmap(seu, features = top5$gene) + scale_fill_viridis(option = "D") + theme(axis.text.y = element_text(size = 6))
dev.off()

cell_markers<- c("EPCAM", "TFF3", "ALB", "DCN", "ACTA2", "VWF", "CD14", "CD68", "CD3D", "CD8A","NKG7","CD79A","CD19")
epi_markers <- c("EPCAM","CLDN4","KRT7","KRT17","KRT18","KRT19","MUC1","MUC6","MUC5AC","PGC","TFF1","TFF2","TFF3")
immune_markers <- c("CD3D", "CD4","IL7R", "FOXP3", "IL2RA","CD8A", "CD8B","GNLY", "NKG7","TRDC", "SELL", "CCR7", "CCR6", "KLRB1", "CD19", "MS4A1", "SDC1")

pdf("cell_markers.pdf")
DoHeatmap(seu, features=cell_markers) + scale_fill_viridis(option="D")
dev.off()

pdf("epi_markers.pdf")
DoHeatmap(seu, features=epi_markers) + scale_fill_viridis(option="D")
dev.off()

pdf("immune_markers.pdf")
DoHeatmap(seu, features=immune_markers) + scale_fill_viridis(option="D")
dev.off()

