library(ape)
library(Seurat)
library(dplyr)
library(RColorBrewer)

seurat_scDNA <- readRDS('seurat_scDNA.rds')

celltype <- seurat_scDNA$celltype
seurat_scDNA_filter <- subset(seurat_scDNA, cells= union(which(celltype == 'normal'), 
                                                         which(celltype == 'G0G1')))

tmp_matrix <- t(as.matrix(seurat_scDNA_filter@assays$RNA@counts))
tmp_matrix[tmp_matrix == 0] <- NA

tmp_matrix <- tmp_matrix[colnames(seurat_scDNA_filter), ]
tmp_matrix <- as.data.frame(tmp_matrix)

tmp_matrix$ident <- seurat_scDNA_filter@meta.data$subclone

tmp_matrix <- tmp_matrix %>% group_by(ident) %>% summarise_all(list(median), na.rm=T)
ident <- tmp_matrix$ident
tmp_matrix$ident <- NULL
tmp_matrix <- as.data.frame(t(tmp_matrix))
tmp_matrix <- na.omit(tmp_matrix)

# Compute the distance for subclones
tmp_matrix <- t(round(tmp_matrix))
rownames(tmp_matrix) <- ident
tmp_matrix <- tmp_matrix[sort(as.character(ident)), ]

dist <- stats::dist(tmp_matrix)

# Construct hierarchical clusters for subclones
hc <- hclust(dist)

mycol = brewer.pal(9, 'Set1')
mycol_nonormal = c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666")
mycol_normal = c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666", brewer.pal(5, 'Set1'))


pdf('phylogeny.pdf')
plot.phylo(as.phylo(hc), main='Phylogenic Tree', direction = 'downwards',tip.color = c(mycol_normal[sort(hc$order)[1:length(hc$labels)-1]],'#252525'), cex=1.5, edge.width = 2)
dev.off()

pdf('phylogeny_unrooted.pdf')
plot.phylo(as.phylo(hc), main='Phylogenic Tree', direction = 'downwards', tip.color = c(mycol_normal[sort(hc$order)[1:length(hc$labels)-1]],'#252525'), cex=1.5, edge.width = 2, type = 'unrooted')
dev.off()
saveRDS(hc, file = 'hc.rds')

## Plot no normal clusters
tmp_no_normal <- tmp_matrix[1:dim(tmp_matrix)[1]-1, ]
dist <- dist(tmp_no_normal)
hc_no <- hclust(dist)
saveRDS(hc_no, file = 'hc_no_normal.rds')

