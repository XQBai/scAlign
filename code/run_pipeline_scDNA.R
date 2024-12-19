library(Seurat)
library(optparse)
source('/mnt/ix1/Projects/M070_200622_GI_multiomics/GithubCode_202411/code/Filter_noise_scDNA.R')
source('/mnt/ix1/Projects/M070_200622_GI_multiomics/GithubCode_202411/code/Identify_cell_components_scDNA.R')
source('/mnt/ix1/Projects/M070_200622_GI_multiomics/GithubCode_202411/code/Construct_subclone_scDNA.R')

option_list = list(
  make_option(c("-p", "--project_path"), action="store", default=NA, type='character',
              help="Project name, default= sample name")
)
opt = parse_args(OptionParser(option_list=option_list))

## read the cell metrics file conducted by cellranger-dna
per_cell_metrics <- read.csv(file.path(opt$project_path, 'per_cell_summary_metrics.csv'), sep=",", header=T)
seurat_scDNA <- readRDS('seurat_scDNA.rds')
scdna_matrix_locs <- readRDS('scdna_matrix_locs.rds')

seurat_scDNA <- Identify_cellranger_noise(per_cell_metrics, seurat_scDNA)
seurat_scDNA <- Identify_technoise(seurat_scDNA)
seurat_scDNA <- ReIdentify_replicates(seurat_scDNA)
seurat_scDNA <- Identify_normal(seurat_scDNA, scdna_matrix_locs)
seurat_scDNA <- Identify_replicates(seurat_scDNA, scdna_matrix_locs)

seurat_scDNA <- Construct_subclones(seurat_scDNA)

saveRDS(seurat_scDNA, file = 'seurat_scDNA.rds')

