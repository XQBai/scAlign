library(dplyr)
library(Seurat)
library(viridis)
library(cowplot)
library(future)

# detach(package:plyr)

ndims=100

#epi_genes <- c("EPCAM","CLDN4","GAST","KRT7","KRT17","KRT18","KRT19","MUC1","MUC6","MUC5AC","PGC","TFF1","TFF2","TFF3",
#               "GIF", "CHGA", "LGR5")
#nonepi_genes <- c("ACTA2","CD4","CD3D","CD3E","CD3G","CD14","CD19","CD34","CD79","CD83","CD247","DCN","IGHC","IGLC","IL2RA",
#                  "MS4A1","PECAM1","PTPRC","SPARC","TRAC","TRBC","THY1","VWF",
#				  "FCGR3A","NKG7","GNLY","CD8A")

# Epithelial genes
epi_gastric <- c("MUC5AC", "PGC", "MUC6", "TFF1", "TFF2")
epi_keratin <- c("KRT7","KRT17","KRT18","KRT19")
epi_tumor <- c("EPCAM", "TFF3", "CLDN4")
epi_antrum <- c("GAST", "PDX1")
epi_intestine <- c("MUC2")
epi_stem_cells <- c("LGR5", "TROY")
epi_parietal <- c("GIF")
epi_neoendocrine <- c("CHGA","GAST","SST")

# Non-epithelial genes
fibro_myof <- c("ACTA2","DCN","SPARC","THY1","LUM","COL1A1","COL14A1","CAV1")
endothelium <- c("VWF","PECAM1","SELE","SELP","ENG")
dendritic_cells <- c("CD83","IL3RA","CLEC4C","HLA-DQA2","ID2","IRF5")
mast_cells <- c("MS4A2","TPSAB1","CPA3","TPSB2")
mono_macro <- c("CD14","FCGR3A","CD163","CD68")
pericytes <- c("PDGFRB","RGS5")

immune_marker <- c("PTPRC")
b_plasma <- c("CD19","CD79","MS4A1","IGHA","IGHC","IGLC","SDC1")
t_cells <- c("CD3D","TRAC","TRAB","CD3E","CD3G")
cd8_tcells <- c("CD8A","CD8B")
t_cytotoxic <- c("NKG7","GNLY","CTSW","GZMA","GZMB","GZMH","KLRB1","KLRD1","KLRK1","PRF1")
nk_cells <- c("NKG7","GNLY","XCL2","NCR1")
t_reg <- c("IL2RA","FOXP3")
memory_t <- c("CCR7","SELL","CD69","ITGAE","ITGA1")
t_general <- c(t_cells[1:3], cd8_tcells[1], nk_cells[1:2], t_reg[1])

# Consolidate gene markers into epithelial and non-epithelial lists for module scoring
epi_genes = c(epi_gastric, epi_keratin, epi_tumor, epi_antrum, epi_intestine, epi_stem_cells, epi_parietal, epi_neoendocrine)

# Non-epithelial list not including any mast cell or pericyte markers currently
nonepi_genes = c(fibro_myof[1:4], endothelium[1:2], dendritic_cells[1], mono_macro[1:2],
                 immune_marker, b_plasma[1:6], t_general)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=12)

seu_merged <- readRDS('seurat_aggr.rds')
      
#Subset epithelial cells based on module score for epithelial genes  
seu_merged <- AddModuleScore(seu_merged, features = list(epi_genes), name = "epithelial")
seu_merged <- AddModuleScore(seu_merged, features = list(nonepi_genes), name = "nonepithelial")
seu_merged <- AddModuleScore(seu_merged, features = list(b_plasma), name='b_plasma')
seu_merged <- AddModuleScore(seu_merged, features = list(t_general), name='t_cells')

epi_scores <- data.frame(id=seu_merged@active.ident, score=seu_merged$epithelial1, score2=seu_merged$nonepithelial1,
                         score3=seu_merged$b_plasma1, score4=seu_merged$t_cells1)

epi_scores_summary <- epi_scores %>% dplyr::group_by(id) %>% 
  dplyr::summarize(me=mean(score), mn=mean(score2), mb=mean(score3), mt=mean(score4)) %>% as.data.frame()

write.csv(epi_scores_summary, 'module_scores.csv')

# ne_cutoff is minimum positive non-epithelial module score per cluster, for clusters with epithelial module score < 0
ne_cutoff <- epi_scores_summary %>% filter(me < 0) %>% filter(mn > 0) %>% arrange(mn) %>% head(1) %>% select(mn) %>% as.numeric()
epi_clusters <- epi_scores_summary %>% filter(me > 0.01) %>% filter(mn < ne_cutoff) %>% select(id) %>% .$id %>% as.character()
  
seu_epi <- subset(seu_merged, idents=epi_clusters)
  
#seu_epi <- SCTransform(seu_epi, vars.to.regress = "percent.mt", verbose=T)
seu_epi <- SCTransform(seu_epi, vars.to.regress='nCount_RNA', return.only.var.genes=FALSE, verbose=T)
seu_epi <- RunPCA(seu_epi, verbose=T, npcs = ndims, approx=F)
seu_epi <- RunUMAP(seu_epi, dims = 1:ndims, verbose=T)
seu_epi <- FindNeighbors(seu_epi, dims=1:ndims, verbose=T)
seu_epi <- FindClusters(seu_epi, algorithm = 4)
seu_epi <- BuildClusterTree(seu_epi, reorder=T, reorder.numeric = T, dims=1:ndims)

saveRDS(seu_epi, 'seurat_epi.rds')
cellnumbers <- table(Idents(seu_epi), seu_epi@meta.data$condition)
write.csv(cellnumbers, file='cell_numbers.csv')

pdf('umap_clusters.pdf')
  DimPlot(seu_epi, label=T)
  DimPlot(seu_epi, group.by='condition', label=T)
dev.off()


