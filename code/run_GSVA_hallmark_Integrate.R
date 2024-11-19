library(Seurat)
library(dplyr) 
library(Matrix)
library(GSEABase)
library(GSVA)
library(tidyverse)
require(reshape2)
library(viridis)
library(broom)
library(future)
options(future.globals.maxSize = 6000 * 1024^2)
plan('multiprocess', workers= 18)


# Define some functions for gsva analysis
GSVAtidy <- function(GSVA_rds, seu_) {
  GSVA_tidy <- as_tibble(GSVA_rds, rownames = "pathway")
  GSVA_tidy <- melt(GSVA_tidy)
  GSVA_tidy <- setNames(GSVA_tidy, c("pathway", "cell_barcode", "value"))
  GSVA_tidy$ident <- plyr::mapvalues(GSVA_tidy$cell_barcode, from=colnames(seu_), to=as.character(Idents(seu_)))
  return(GSVA_tidy)
}


GSVAIdentCompare <- function(matGSVA) {
  GSVAtest <- matGSVA %>% 
    nest_by(pathway) %>%
    mutate(multitst = list(TukeyHSD(aov(value ~ ident, data = data)))) %>% 
    summarise(tidy(multitst))
  return(GSVAtest)
}

GSVA_scaled <- function(hallmark_test, seu){
  
  #scale, filter for significant pathways
  hallmark_test <- dplyr::filter(hallmark_test, adj.p.value < 0.05)
  
  #averaging across seurat object clusters, scaling
  cellInfo <- data.frame(seuratCluster=Idents(seu))
  hallmark_activity <-  sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                               function(cells) rowMeans(pq_hallmark_df[,cells]))
  
  #scale, filter for significant pathways
  hallmark_scaled <- as.data.frame(t(scale(t(hallmark_activity), center = T, scale=T))) %>% rownames_to_column("pathway") %>%
    dplyr::filter(pathway %in% hallmark_test$pathway) %>% column_to_rownames(var = "pathway")
  
  return(hallmark_scaled)
}

seu <- readRDS('seurat_epi.rds')
Idents(seu) <- seu@meta.data$corrlabel
# RNAmatrix <- seu@assays$RNA@data
RNAmatrix <- seu@assays$RNA@counts

hallmark <- getGmt("/mnt/ix1/Resources/scRNA_Ref/GSEA_Genesets/h.all.v6.2.symbols.gmt")
# pq_hallmark <- gsva(as.matrix(RNAmatrix), hallmark, kcdf="Gaussian",  mx.diff=T, verbose=FALSE, parallel.sz=3, min.sz=10)
pq_hallmark <- gsva(as.matrix(RNAmatrix), hallmark, kcdf="Poisson",  mx.diff=T, verbose=FALSE, parallel.sz=3, min.sz=10)
saveRDS(pq_hallmark, file="pq_hallmark.rds")

rownames(pq_hallmark) <- gsub("HALLMARK_", "", rownames(pq_hallmark)) 
pq_hallmark_df <- as.data.frame(pq_hallmark)
#tidying outputs and performing ANOVA with Tukey's HSD across clusters
hallmark_tidy <- GSVAtidy(pq_hallmark, seu)
hallmark_test <- GSVAIdentCompare(hallmark_tidy)
write.csv(hallmark_test, file="hallmark_test.csv")
hallmark_scaled <- GSVA_scaled(hallmark_test, seu)

saveRDS(hallmark_scaled, file = 'hallmark_scaled.rds')

#plotting
pdf("hallmark_gmt.pdf")
pheatmap::pheatmap(hallmark_scaled, fontsize =20, 
                   color=viridis(20),
                   fontsize_col = 15,
                   fontsize_row = 8
                   #cluster_rows = FALSE, cluster_cols = FALSE
)
dev.off()

#grabbing top_n to plot as above 
top3 <- rownames_to_column(hallmark_scaled, var = "pathway") %>% reshape2::melt()%>% group_by(variable)%>% top_n(3, value)

top3_hallmark_scaled <- as.data.frame(hallmark_scaled) %>% rownames_to_column(var="pathway")%>% dplyr::filter(pathway %in% top3$pathway)%>% column_to_rownames(var = "pathway")

pdf("hallmark_gmt_top3.pdf")
pheatmap::pheatmap(top3_hallmark_scaled, fontsize =15, 
                   color=viridis(20)
                   #cluster_rows = FALSE, cluster_cols = FALSE
)
dev.off()





# #############################################################################
# # Use sepcific csv file
# # Use reads counts matrix and DE genes as GSVA input
# RNAmatrix <- GetAssayData(object = seu[["RNA"]], slot = "counts")
# seu <- FindVariableFeatures(seu, assay = 'RNA', selection.method = 'vst', nfeatures = 2000)
# var_genes <- seu[['RNA']]@var.features
# RNAmatrix <- RNAmatrix[var_genes, ]
# 
# metabolic_genes <- read.table('/mnt/ix1/Resources/scRNA_Ref/GSEA_Genesets/cd8_exhaustion_cellcycle_cytotoxic.csv', header=T, sep=",", stringsAsFactors=F)
# metabolic_pways <- unique(metabolic_genes$Pathway)
# metabolic_pways_dtls <- lapply(metabolic_pways, function(x) metabolic_genes[metabolic_genes$Pathway==x,1])
# names(metabolic_pways_dtls) <- metabolic_pways
# pq_metabolic <- gsva(as.matrix(RNAmatrix), metabolic_pways_dtls, kcdf="Gaussian",  mx.diff=T, verbose=FALSE, parallel.sz=3, min.sz=10)
# saveRDS(pq_metabolic, file="pq_metabolic.rds")
# 
# pq_metabolic_df <- as.data.frame(pq_metabolic)
# #tidying outputs and performing ANOVA with Tukey's HSD across clusters
# goldrath_tidy <- GSVAtidy(pq_metabolic, seu)
# saveRDS(goldrath_tidy, file = 'goldrath_tidy.rds')
# goldrath_test <- GSVAIdentCompare(goldrath_tidy)
# write.csv(goldrath_test, file="goldrath_test.csv")
# 
# #scale, filter for significant pathways
# goldrath_test <- dplyr::filter(goldrath_test, adj.p.value < 0.05)
# 
# #averaging across seurat object clusters, scaling
# #Idents(seu) <- seu@meta.data$new.cluster
# cellInfo <- data.frame(seuratCluster=Idents(seu))
# goldrath_activity <-  sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
#                              function(cells) rowMeans(pq_metabolic_df[,cells]))
# 
# #scale, filter for significant pathways
# goldrath_scaled <- as.data.frame(t(scale(t(goldrath_activity), center = T, scale=T))) %>% rownames_to_column("pathway") %>%
#   dplyr::filter(pathway %in% goldrath_test$pathway) %>% column_to_rownames(var = "pathway")
# 
# #plotting
# pdf("goldrath_cd8.pdf")
# pheatmap::pheatmap(goldrath_scaled, fontsize =15, 
#                    color=viridis(20),
#                    cluster_rows = FALSE, cluster_cols = FALSE)
# dev.off()
# 
# #grabbing top_n to plot as above 
# top3 <- rownames_to_column(goldrath_scaled, var = "pathway") %>% reshape2::melt()%>% group_by(variable)%>% top_n(3, value)
# top3_goldrath_scaled <- as.data.frame(goldrath_scaled) %>% rownames_to_column(var="pathway")%>% dplyr::filter(pathway %in% top3$pathway)%>% column_to_rownames(var = "pathway")
# 
# pdf("goldrath_cd8_top3.pdf")
# pheatmap::pheatmap(top3_goldrath_scaled, fontsize =15, 
#                    color=viridis(20),
#                    cluster_rows = FALSE, cluster_cols = FALSE)
# dev.off()


