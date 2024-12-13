# Single-cell-multi-omics-integration

Single-cell multi-omics data analysis pipeline for defining the subclone architecture and clonal phenotype from the primary gastric and metastatic colon cancers. We conducted scDNA-seq and scRNA-seq assays on thousands of single cells derived from the same specimen. We constructed the subclone architecture by scDNA-seq data and characterized the distinct large CNVs among subclones. Then we analyzed the paired scRNA-seq gene expression profile to determine cell types and transcriptome features. Utilizing both results, we assigned the scRNA-seq gene expression to subclones derived in scDNA-seq based on gene dosage effect. Finally, we identified the differential phenotype pathways among subclones by GSVA. There are three major parts included in the whole pipeline:
## scDNA-seq analysis:
 1. Cell quality control (QC)
 2. Identified cellular components 
 3. Constructed subclones 
 
 ```
 Rscript Preprocessing_data_scDNA.R
 Rscript run_pipeline_scDNA.R
 Rscript Plot_umap_scDNA.R
 Rscript Plot_evolution_scDNA.R
 Rscript Plot_subclones_heatmap_scDNA.R
 ```
## scRNA-seq analysis: 
 1. Cell and gene quality control for each sample by Seurat
 ```
 Rscript Seurat_scRNA.R -m './example/P5931/scRNA/P5931_normal_1/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_1' -r './example/10x_multiplet_rate.csv'
 Rscript Seurat_scRNA.R -m './example/P5931/scRNA/P5931_normal_2/outs/filtered_feature_bc_matrix/' -p 'P5931_normal_2' -r './example/10x_multiplet_rate.csv'
 Rscript Seurat_scRNA.R -m './example/P5931/scRNA/P5931_tumor_1/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_1' -r './example/10x_multiplet_rate.csv'
 Rscript Seurat_scRNA.R -m './example/P5931/scRNA/P5931_tumor_2/outs/filtered_feature_bc_matrix/' -p 'P5931_tumor_2' -r './example/10x_multiplet_rate.csv'
 ```
 
 2. Merged normal and tumor objects
 ```
 Rscript merge_seurat_scRNA.R
 Rscript merge_sctransform_scRNA.R
 Rscript find_markers_scRNA.R
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
 Rscript Projection_alignment_Integrate.R
 ```
 2. Evaluated the subclone assignment by inferCNV package 
 ```
 Rscript run_infercnv_Integrate.R
 ```
 3. Discovered phenotype biology of subclones  
 ```
 Rscript run_GSVA_hallmark_Integrate.R
 ```
## Data
The scDNA-seq and scRNA-seq datasets generated for this study are available in NCBI's dbGAP repositories, accession numbers phs001711 and phs001818. 
<!--
## Reference
[Single cell multi-omic mapping of subclonal architecture and pathway phenotype in primary gastric and metastatic colon cancers, bioRxiv](https://www.biorxiv.org/content/10.1101/2022.07.03.498616v1)
-->

## Contact 
### xiangqi@stanford.edu 
