# run infercnv when use one of clutser in seurat object as normal cells (reference)
library(Seurat)
library(optparse)
library(dplyr)
library(viridis)
library(cowplot)
library(future)
library(infercnv)

seurat_obj <- readRDS('seurat_epi.rds')
proj_alignment <- seurat_obj$pro_alignment

normal_cells <- colnames(seurat_obj)[which(seurat_obj$pro_alignment == 'normal')]
normal_select <- sample(normal_cells, 500)

cell_select <- c(normal_select, colnames(seurat_obj)[which(seurat_obj$pro_alignment != 'normal')])

seurat_obj <- subset(seurat_obj@assays$RNA, cells = cell_select)

counts_matrix = GetAssayData(seurat_obj, assay="RNA", slot='counts')

# # a raw counts matrix of single-cell RNA-Seq expression
# DefaultAssay(seurat_obj) <- 'RNA'
# counts_matrix = GetAssayData(seurat_obj, slot='counts')

label <- as.character(proj_alignment[colnames(seurat_obj)])

# label <- as.character(seurat_obj$pro_alignment)
# # Cell annotation (use cluster# from original sample)
cell_annot <- as.data.frame(label)
cell_annot$barcode <- colnames(seurat_obj)
write.table(cell_annot[,c(2,1)], 'idents.txt', sep="\t", row.names=F, col.names=F)

if ('normal' %in% unique(label)){
  # create the infercnv object 
  infercnv_obj = CreateInfercnvObject(counts_matrix, 
                                      annotations_file= "idents.txt", 
                                      delim="\t", gene_order_file="../gene_locs.sorted.bed", 
                                      ref_group_names='normal')
}else{
  # create the infercnv object 
  infercnv_obj = CreateInfercnvObject(counts_matrix, 
                                      annotations_file= "idents.txt", 
                                      delim="\t", gene_order_file="../gene_locs.sorted.bed", 
                                      ref_group_names=NULL)
}

# perform infercnv operations to reveal cnv signal
infercnv_obj_final = infercnv::run(infercnv_obj, 
                                   cutoff=.1, 
                                   out_dir="test", 
                                   cluster_by_groups=T, 
                                   num_threads=40, 
                                   denoise=T,
                                   HMM=T)
