# scAlign

scAlign: Single-cell clonal alignment for integrating scRNA-seq and scDNA-seq from the same specimen. scAlign is a single-cell multi-omics integration framework designed to resolve clonal architecture and phenotypic diversity. First, scAlign constructs subclone architecture using scDNA-seq data and characterizes large-scale copy number variations across subclones. Second, it assigns scRNA-seq cells to matched subclones based on gene dosage effects. Overall, it enables integrated analysis of gene expression and genomic alterations to define clonal phenotypes and identify differential pathways among subclones. We applied scAlign to thousands of single cells from primary gastric cancer and metastatic colon cancer specimens, demonstrating its ability to resolve subclonal genomic structure, assign transcriptomic profiles to subclones, and characterize clonal phenotypic heterogeneity.

## Installation

```
install.packages('devtools')
devtools::install_github("XQBai/scAlign")
```

## Download the example data 
Download gastric cancer sample P5931 scDNA-seq and scRNA-seq data from [P5931](https://github.com/XQBai/Single-cell-multi-omic-integration/releases/tag/P5931). 
ScDNA-seq and scRNA-seq of all samples can be downloaded from [Ji Research Group](https://dna-discovery.stanford.edu/research/datasets/).

## Example Workflow for sample P5931 
```
library(scAlign)
library(Seurat)
library(tidyverse)

# create output directories
dir.create('./output/scDNA', showWarnings = FALSE, recursive = TRUE)
dir.create('./output/scRNA', showWarnings = FALSE, recursive = TRUE)

# define input paths
data_dir <- './data/P5931'
scdna_dir <- file.path(data_dir, 'scDNA')
scrna_dir <- file.path(data_dir, 'scRNA')

genome_reference_file <- file.path(scdna_dir, 'GRCh38_cellranger_20k.canonical.rownumbers.bed')
scdna_tsv_file <- file.path(scdna_dir, 'P5931_801_scdna_matrix.tsv')
barcode_txt_file <- file.path(scdna_dir, 'P5931_801_tumor.barcodes.txt')
per_cell_metrics_file <- file.path(scdna_dir, 'per_cell_summary_metrics.csv')
genes_bin_file <- './data/genes_bin.RData'
gene_locs_path <- './data/gene_locs.sorted.bed'

```

### scAlign construct subclones of scDNA-seq data:
```
# merge genome bins as 1MB segments
scdna_res <- scdna_preprocessing(genome_reference_file, scdna_tsv_file, barcode_txt_file)
scdna_matrix_merge <- scdna_res[[1]]
scdna_matrix_locs <- scdna_res[[2]]

saveRDS(scdna_matrix_merge, './output/scDNA/scdna_matrix_all_barcodes.rds')
saveRDS(scdna_matrix_locs, './output/scDNA/scdna_matrix_locs.rds')

# Identify cellular components
seurat_scDNA <- create_seurat_scdna(scdna_matrix_merge, scdna_matrix_locs)
seurat_scDNA <- identify_cellranger_noise(per_cell_metrics_file, seurat_obj = seurat_scDNA)
seurat_scDNA <- identify_technical_noise(seurat_obj = seurat_scDNA)
seurat_scDNA <- identify_normal_cells(seurat_scDNA, scdna_matrix_locs)
seurat_scDNA <- identify_replication_cells(seurat_scDNA, scdna_matrix_locs)

# Construct subclones and visualize subclones 
seurat_scDNA <- construct_subclones(seurat_scDNA)
plot_subclonal_heatmap(seurat_scDNA, scdna_matrix_locs, celltype = 'G0G1',
                       output_path = './output/scDNA/G0G1_subclone_heatmap.pdf')
# Generate gene-level subclone matrix for preparation of scAlign integration
scdna_gene_subclones <- generate_subclone_cnv_gene_matrix(
  seurat_scDNA,
  genes_bin_file,
  scdna_tsv_file,
  barcode_txt_file
)
saveRDS(scdna_gene_subclones, './output/scDNA/scdna_gene_subclones.rds')
```
## scRNA-seq analysis: 

 ```
# Seurat preprocessing each sample for quality control, doublets filtering, and cell cycle assignment
Rscript ./script/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_normal_1/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_1' -r './data/10x_multiplet_rate.csv'
Rscript ./script/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_normal_2/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_2' -r './data/10x_multiplet_rate.csv'
Rscript ./script/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_tumor_1/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_1' -r './data/10x_multiplet_rate.csv'
Rscript ./script/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_tumor_2/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_2' -r './data/10x_multiplet_rate.csv'

# Merge scRNA-seq data and subset epithelial G0/G1 cells
rds_files <- list.files('./output/scRNA', pattern = '*seurat_dfx.rds', full.names = TRUE)
sample_list <- list('P5931_normal_1', 'P5931_normal_2', 'P5931_tumor_1', 'P5931_tumor_2')
condition_list <- list('normal', 'normal', 'tumor', 'tumor')

seu_epi <- run_scrna_pipeline(rds_files, sample_list, condition_list)
saveRDS(seu_epi, './output/scRNA/seurat_epi.rds')

# Normalize scRNA-seq matrix for preparing input for scAlign integration
scrna_normalized_matrix <- run_scrna_normalization(seu_epi, gene_locs_path)
saveRDS(scrna_normalized_matrix, './output/scRNA/RNA_normalized_matrix.rds')
 ```
 
 2. Merged normal and tumor objects
 ```
 Rscript merge_seurat_scRNA.R
 Rscript merge_sctransform_scRNA.R
 Rscript find_markers_scRNA.R -r 'seurat_aggr.rds'
 ```
 
 3. Extracted the epithelial cell and selected the G0/G1 phase cells
 ```
 Rscript epi_subset_scRNA.R
 Rscript split_Gphase_scRNA.R
 ```
 
## Integration analysis 
 1. Normalize the gene expression of scRNA-seq in G0/G1 Phase
 2. Assigned epithelial single cells (G0G1 phase) from scRNA-seq into subclones of scDNA-seq
 ```
 Rscript scRNA_pre_processing_Integrate.R
 Rscript Avg_CNV_subclones.R 
 Rscript Projection_alignment_Integrate.R
 ```
 3. Evaluated the subclone assignment by inferCNV package 
 ```
 Rscript run_infercnv_Integrate.R
 ```
 4. Discovered phenotype biology of subclones  
 ```
 Rscript run_pathway_Integrate.R
 ```
## Data
The scDNA-seq and scRNA-seq datasets generated for this study are available in NCBI's dbGAP repositories, accession numbers phs001711 and phs001818. 
<!--
## Reference
[Single cell multi-omic mapping of subclonal architecture and pathway phenotype in primary gastric and metastatic colon cancers, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.03.498616v1)
-->

## Contact 
### xiangqi@stanford.edu 
