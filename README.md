# scAlign

scAlign: Single-cell clonal alignment for integrating scRNA-seq and scDNA-seq from the same specimen. scAlign is a single-cell multi-omics integration framework designed to resolve clonal architecture and phenotypic diversity. First, scAlign constructs subclone architecture using scDNA-seq data and characterizes large-scale copy number variations across subclones. Second, it assigns scRNA-seq cells to matched subclones based on gene dosage effects. Overall, it enables integrated analysis of gene expression and genomic alterations to define clonal phenotypes and identify differential pathways among subclones. We applied scAlign to thousands of single cells from primary gastric cancer and metastatic colon cancer specimens, demonstrating its ability to resolve subclonal genomic structure, assign transcriptomic profiles to subclones, and characterize clonal phenotypic heterogeneity. The workflow in scAlign is illustrated in the figure below. ![](https://github.com/XQBai/scAlign/blob/main/inst/extdata/scAlign_Workflow.png)

## Installation

```
install.packages('devtools')
devtools::install_github("XQBai/scAlign")
```

## Download the example data 
Download gastric cancer sample P5931 scDNA-seq and scRNA-seq data from [P5931](https://github.com/XQBai/Single-cell-multi-omic-integration/releases/tag/P5931). 
ScDNA-seq and scRNA-seq data for all samples can be downloaded from the [Ji Research Group](https://dna-discovery.stanford.edu/research/datasets/).

## Example Workflow for sample P5931 
```
library(scAlign)
library(Seurat)
library(tidyverse)

# create output directories
dir.create('./output/scDNA', showWarnings = FALSE, recursive = TRUE)
dir.create('./output/scRNA', showWarnings = FALSE, recursive = TRUE)

# download P5931's data and define input paths
data_dir <- './data/P5931'
genome_ref_dir <- './inst/genome_reference'
scdna_dir <- file.path(data_dir, 'scDNA')
scrna_dir <- file.path(data_dir, 'scRNA')

genome_reference_file <- file.path(genome_ref_dir, 'GRCh38_cellranger_20k.canonical.rownumbers.bed')
scdna_tsv_file <- file.path(scdna_dir, 'P5931_801_scdna_matrix.tsv')
barcode_txt_file <- file.path(scdna_dir, 'P5931_801_tumor.barcodes.txt')
per_cell_metrics_file <- file.path(scdna_dir, 'per_cell_summary_metrics.csv')
genes_bin_file <- file.path(genome_ref_dir, 'genes_bin.RData')
gene_locs_path <- file.path(genome_ref_dir, 'gene_locs.sorted.bed')

# define output folders for scDNA and scRNA
output_dir_scdna <- './output/scDNA'
output_dir_scrna <- './output/scRNA'
```

### scAlign construct subclones of scDNA-seq data:
```
# merge genome bins as 1MB segments
scdna_res <- scdna_preprocessing(genome_reference_file, scdna_tsv_file, barcode_txt_file, output_dir = output_dir_scdna)
scdna_matrix_merge <- scdna_res[[1]]
scdna_matrix_locs <- scdna_res[[2]]

# Identify cellular components
seurat_scDNA <- create_seurat_scdna(scdna_matrix_merge, scdna_matrix_locs, output_dir = output_dir_scdna)
seurat_scDNA <- identify_cellranger_noise(per_cell_metrics_file, seurat_obj = seurat_scDNA)
seurat_scDNA <- identify_technical_noise(seurat_obj = seurat_scDNA, output_dir = output_dir_scdna)
seurat_scDNA <- identify_normal_cells(seurat_scDNA, scdna_matrix_locs)
seurat_scDNA <- identify_replication_cells(seurat_scDNA, scdna_matrix_locs, output_dir = output_dir_scdna)

# Construct subclones and visualize subclones 
seurat_scDNA <- construct_subclones(seurat_scDNA, output_dir = output_dir_scdna)
plot_subclonal_heatmap(seurat_scDNA, scdna_matrix_locs, celltype = 'G0G1',
                       output_pdf = 'subclone_heatmap.pdf', output_dir = output_dir_scdna)
# Generate gene-level subclone matrix for preparation of scAlign integration
scdna_gene_subclones <- generate_subclone_cnv_gene_matrix(
  seurat_scDNA,
  genes_bin_file,
  scdna_tsv_file,
  barcode_txt_file,
  output_dir = output_dir_scdna
)
```
### scAlign analyze scRNA-seq data
 ```
# Seurat preprocessing each sample for quality control, doublets filtering, and cell cycle assignment
Rscript ./inst/scripts/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_normal_1/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_1' -r './data/10x_multiplet_rate.csv'
Rscript ./inst/scripts/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_normal_2/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_2' -r './data/10x_multiplet_rate.csv'
Rscript ./inst/scripts/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_tumor_1/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_1' -r './data/10x_multiplet_rate.csv'
Rscript ./inst/scripts/Seurat_scRNA.R -m './data/P5931/scRNA/P5931_tumor_2/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_2' -r './data/10x_multiplet_rate.csv'

# Merge scRNA-seq data and subset epithelial G0/G1 cells
rds_files <- list.files('./output/scRNA', pattern = '*seurat_dfx.rds', full.names = TRUE)
sample_list <- list('P5931_normal_1', 'P5931_normal_2', 'P5931_tumor_1', 'P5931_tumor_2')
condition_list <- list('normal', 'normal', 'tumor', 'tumor')

seu_epi <- run_scrna_pipeline(rds_files, sample_list, condition_list, output_dir = output_dir_scrna)

# Normalize scRNA-seq matrix for preparing input for scAlign integration
scrna_normalized_matrix <- run_scrna_normalization(seu_epi, gene_locs_path, output_dir = output_dir_scrna)

 ```
### scAlign assigns scRNA-seq cells into subclones detected by scDNA-seq 
Assigned epithelial single cells (G0G1 phase) from scRNA-seq into subclones of scDNA-seq
 ```
# load inputs 
scdna_gene_matrix <- readRDS('./output/scDNA/scdna_gene_matrix.rds')
scrna_normalized_matrix <- readRDS('./output/scRNA/RNA_normalized_matrix.rds')
scdna_subclones <- readRDS('./output/scDNA/scdna_gene_subclones.rds')
seu_epi <- readRDS('./output/scRNA/seurat_epi.rds')
seurat_obj <- readRDS('./output/scDNA/seurat_scDNA.rds')

# Run scAlign integration
# define signal_chr parameter according to chromosomes with significant CNV changes of scDNA-seq subclones
seu_epi <- run_scalign_pipeline(
  scrna_normalized_matrix,
  scdna_gene_matrix,
  scdna_subclones,
  seurat_obj,
  seu_epi,
  gene_locs_path,
  signal_chr = c(3, 7, 8, 21), 
  output_seurat_rds = 'seurat_epi.rds',
  output_alignment_csv = 'project_alignment.csv',
  method = 'euclidean',
  output_dir = output_dir_scrna
)
print(table(seu_epi$pro_alignment))

# Evaluated the subclone assignment by the inferCNV package
Rscript ./inst/scripts/run_infercnv_Integrate.R
 ```
## Data
The scDNA-seq and scRNA-seq datasets generated for this study are available in NCBI's dbGAP repositories, accession numbers phs001711 and phs001818. 

## Reference
Single-cell aneuploidy and chromosomal arm imbalances define subclones with divergent transcriptomic phenotypes, under revision

<!--
## Reference
[Single cell multi-omic mapping of subclonal architecture and pathway phenotype in primary gastric and metastatic colon cancers, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.03.498616v1)
-->

## Contact 
xiangqi@stanford.edu 
