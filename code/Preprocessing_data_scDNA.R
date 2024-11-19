library(Seurat)
library(data.table)
library(dplyr)
library(Matrix)
library(mixtools)
library(optparse)
library(future)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=24)

## The folder which scDNA-seq matrix located 
project_path <- '../example/P5931/scDNA'

# input the bins genome reference
genome_reference <- fread(file.path(project_path, 'GRCh38_cellranger_20k.canonical.rownumbers.bed'), header=F, data.table=F)

# input the cells by bins matrix 
scdna_matrix <- fread(file.path(project_path, 'P5931_801_tumor.tsv'))

# input the cell barcode list 
barcodes <- read.csv(file.path(project_path, 'P5931_801_tumor.barcodes.txt'), sep='\n', header = FALSE)
colnames(scdna_matrix) <- as.character(barcodes$V1)

# number of 20kb bins, the convert the neighborhood 
# 50 bins as 1Mb segments  
bin_size <- 50 

# Merge 50 bins 
genome_reference <- genome_reference %>% group_by(V1) %>% mutate(bin=floor(row_number(V1)/bin_size))
genome_reference$bin_corrected <- cumsum(c(0,as.numeric(diff(genome_reference$bin))!=0))

scdna_matrix$chr <- genome_reference$V1
scdna_matrix$bin <- genome_reference$bin_corrected

# Apply the median of CNVs as the CNV value on 1Mb segment
scdna_matrix_merge <- scdna_matrix %>% group_by(chr,bin) %>% summarise_all(list(median))
scdna_matrix_locs <- scdna_matrix_merge[,c("chr","bin")]
scdna_matrix_locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY')

scdna_matrix_merge <- scdna_matrix_merge %>% filter(chr != 'chrX' & chr != 'chrY')

scdna_matrix_merge$chr <- NULL
scdna_matrix_merge$bin <- NULL

row.names(scdna_matrix_merge) <- paste0('segement', 1:dim(scdna_matrix_merge)[1])
# Run Seurat to initial merged data 
dims.reduce <- 50
seurat_scDNA <- CreateSeuratObject(counts = as.matrix(scdna_matrix_merge), project='seurat-v3', min.cells = -Inf, min.features = -Inf)

seurat_scDNA <- seurat_scDNA %>% ScaleData(do.center=F, do.scale = F) %>% RunPCA(verbose=T, features=rownames(seurat_scDNA)) %>% RunUMAP(umap.method = 'uwot',n.neighbors = min(20, ncol(seurat_scDNA)), dims=1:dims.reduce) %>% 
  FindNeighbors(reduction='umap', dims=1:2)

scdna_matrix_merge_allbarcodes <- scdna_matrix_merge

# Save the objects for all downstream analysis
# The CNV matrix for merged 1Mb segments
saveRDS(scdna_matrix_merge_allbarcodes, file='scdna_matrix_all_barcodes.rds')
# The genome annotation/reference for 1Mb segments on Whole genome sequencing 
saveRDS(scdna_matrix_locs, file='scdna_matrix_locs.rds')
# The seurat object
saveRDS(seurat_scDNA, file='seurat_scDNA.rds')






