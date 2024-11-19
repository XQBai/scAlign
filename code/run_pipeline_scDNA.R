library(Seurat)

## read the cell metrics file conducted by cellranger-dna
per_cell_metrics <- read.csv(file.path(project_path, 'per_cell_summary_metrics.csv'), sep=",", header=T)
project_path <- '../example/P5931/scDNA'
seurat_scDNA <- readRDS('seurat_scDNA.rds')
scdna_matrix_locs <- readRDS('scdna_matrix_locs.rds')

seurat_scDNA <- Identify_cellranger_noise(per_cell_metrics, seurat_scDNA)
seurat_scDNA <- Identify_technoise(seurat_scDNA)
seurat_scDNA <- ReIdentify_replicates(seurat_scDNA)
seurat_scDNA <- Identify_normal(seurat_scDNA, scdna_matrix_locs)
seurat_scDNA <- Identify_replicates(seurat_scDNA, scdna_matrix_locs)

seurat_scDNA <- Construct_subclones(seurat_scDNA)

saveRDS(seurat_scDNA, file = 'seurat_scDNA.rds')

