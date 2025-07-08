# Desc: Script reads saved Seurat R object and performs PCA dimension reduction #
# Load packages and custom helper functions; parse config.ini parameters if provided      #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
library(Matrix)
library(optparse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(pracma)
library(stringr)
library(ini) 
library(DoubletFinder)
library(viridis)

option_list = list(
  make_option(c("-m", "--fn_prefix"), action="store", default=NA, type='character',
              help="Prefix for cellranger matrix file (<mprefix>.matrix.mtx)"), 
  make_option(c("-p", "--sample"), action="store", default=NA, type='character',
              help="Project name, default= sample name"),
  make_option(c('-r', "--multiplet_rate", action='store', default = NA, type='character',
                help='10X multiplet rate table, default = 10x_multiplet_rate.csv'))
)
opt = parse_args(OptionParser(option_list=option_list))


# Read and write 10X sparse matrix trios
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
read10Xtrio <- function(fn_prefix) {
  
  barcode.path <- paste0(fn_prefix, "barcodes.tsv")
  matrix.path <- paste0(fn_prefix, "matrix.mtx")
  expr_matrix <- readMM(matrix.path)
  
  features_tsv <- paste0(fn_prefix, "features.tsv")
  genes_tsv    <- paste0(fn_prefix, "genes.tsv")
  features.path <- ifelse(file.exists(features_tsv), features_tsv, genes_tsv)
  feature.names = read.delim(features.path, sep="\t", header=F, stringsAsFactors=F)
  rownames(expr_matrix)= feature.names$V2
  
  barcode.names = read.delim(barcode.path, sep="\t", header=F, stringsAsFactors=F) 
  colnames(expr_matrix) = barcode.names$V1
  
  return(expr_matrix)
}

##Annotate with percent mitochondrial gene/UMI counts
annotateMito <- function(seu_, subsetMito=F, max.mito=25, geneID='hugo', organism='Homo sapiens') {
  if (organism == 'Mus musculus') {
    mito_prefix = 'mt-'
  } else {
    mito_prefix = 'MT-'
  }	
  if (geneID == 'both') {
    mito_regex = paste0("^ENS.{12}::", mito_prefix)
  } else {
    mito_regex = paste0("^", mito_prefix)
  }
  
  seu_[["percent.mito"]] <- PercentageFeatureSet(seu_, pattern=mito_regex)
  if (subsetMito) {
    seu_ <- subset(seu_, subset = percent.mito < max.mito)
  }
  return(seu_)
}  

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Set standard default parameters                                             #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
params = list(min_cells=3, min_genes=200, max_genes=5000, max_mito=30, 
              nr_pcs = 20, do_approx = TRUE, cell_cycle='Y', 
              sct_features=2000, sct_regress=NULL, only_var_genes=FALSE,
              dims_reduce=20, cluster_res=0.8, do_reorder='Y', algorithm=1, vst_features=2000,
              remove_doublets ='Y')
metadata = FALSE

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Initialize Seurat object, and write as R object for subsequent analysis     #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
expr_mtx <- read10Xtrio(opt$fn_prefix)
seurat_obj <- CreateSeuratObject(counts=expr_mtx, min.cells=params$min_cells, 
                                 min.features=params$min_genes, project=opt$sample)
length(Idents(seurat_obj))
seurat_obj <- annotateMito(seurat_obj, geneID='hugo')

#Initial QC plots
pdf(paste0(seurat_obj@project.name, '_qc_plots.pdf'))
VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)
plot1 <- FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="percent.mito")
plot2 <- FeatureScatter(seurat_obj, feature1="nCount_RNA", feature2="nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# Filter out low quality cells
seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > params$min_genes & nFeature_RNA < params$max_genes & percent.mito < params$max_mito)
length(Idents(seurat_obj))

### Perform the SCTransform 
seurat_obj <- SCTransform(seurat_obj, vars.to.regress=params$sct_regress, variable.features.n=params$sct_features, 
                          return.only.var.genes=params$only_var_genes)

# Annotate cell cycle
if ( params$cell_cycle == 'Y' ) {
  seurat_obj <- CellCycleScoring(seurat_obj, s.features=cc.genes$s.genes, g2m.features=cc.genes$g2m.genes, set.ident=F)
  seurat_obj$CC.Difference <- seurat_obj$S.Score - seurat_obj$G2M.Score
  cell_cycle_cts <- as.data.frame(table(seurat_obj@meta.data$Phase))
  colnames(cell_cycle_cts) <- c('Phase', 'Cells')
  write.csv(cell_cycle_cts, file=paste0(seurat_obj@project.name,'.cell_cycle.csv'), row.names=F, quote=F)
}

# Run principal component analysis using variable genes
# defaults: npcs: 20; do_approx=T
seurat_pca <- RunPCA(seurat_obj, npcs=params$nr_pcs, approx=params$do_approx)
# saveRDS(seurat_pca, file=paste0(seurat_pca@project.name, '.seurat_pca.rds'))

pdf(paste0(seurat_pca@project.name, '.PCAplot.pdf'))
PCAPlot(seurat_pca)
cells_for_pca = min(dim(seurat_pca[["RNA"]])[2], 500)
DimHeatmap(seurat_pca, dims=1:10, cells=cells_for_pca, balanced=TRUE)
dev.off()


#UMAP dimension reduction and associated plots
seurat_umap <- FindNeighbors(seurat_pca, dims=1:params$dims_reduce)
seurat_umap <- FindClusters(seurat_umap, resolution=params$cluster_res, algorithm=params$algorithm)
if ( params$do_reorder == 'Y' ) {
  seurat_umap <- BuildClusterTree(seurat_umap, reorder=T, reorder.numeric=T, dims=1:params$dims_reduce)
}	
seurat_umap <- RunUMAP(seurat_umap, dims=1:params$dims_reduce, do.label=T)

pdf(paste0(seurat_umap@project.name, '_umap_clusters.pdf'))
DimPlot(seurat_umap, reduction='umap', label=T)
dev.off()

# cellnumbers <- table(Idents(seurat_umap), seurat_umap@meta.data$orig.ident)
# write.csv(cellnumbers, file=paste0(seurat_umap@project.name,'.cellnumbers.csv'))

# Normalize, scale RNA slot in preparation for find markers
seurat_umap <- NormalizeData(seurat_umap, assay="RNA")
seurat_umap <- FindVariableFeatures(seurat_umap, assay="RNA", selection.method="vst", nfeatures=params$vst_features)
seurat_umap <- ScaleData(seurat_umap, assay="RNA", features=rownames(seurat_umap), vars.to.regress="nCount_RNA")

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
# Read Seurat object, annotate doublets, save for further analysis            #
#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
#get pK value
sweep.res.list <- paramSweep_v3(seurat_umap, PCs = 1:20, sct = TRUE)
sweep.stats <- summarizeSweep(sweep.res.list, GT=F)
bcmvn <- find.pK(sweep.stats)
tmp <- bcmvn[with(bcmvn,order(-BCmetric)),]
tmp$pK <- as.character(tmp$pK)
pk_val= as.numeric(tmp[1,2])

# pN (proportion of doublets vs merged real-artificial data) defaults to 0.25; results are largely pN invariant
# pK (PC neighborhood size as a proportion of merged real-artificial data); from find.pK above
#two thresholds for doublets based on cell loading and homotypic doublets. First with cell loading
#only cell loading used
cell_loading <- read.csv(opt$multiplet_rate)
seurat_nr_cells = length(Idents(seurat_umap))
doublet_rate_for_seurat_cells <- cell_loading[which(abs(cell_loading$cells_recovered-seurat_nr_cells)==min(abs(cell_loading$cells_recovered-seurat_nr_cells))),]
droplet_doublet_rate <- doublet_rate_for_seurat_cells[1,1]
nExp_poi <- round(droplet_doublet_rate*seurat_nr_cells)
seurat_umap <- doubletFinder_v3(seurat_umap, PCs = 1:20, pN=0.25, pK=pk_val, nExp=nExp_poi, reuse.pANN=F, sct = TRUE)
names(seurat_umap@meta.data) <- gsub("^.*DF.*$","DF", names(seurat_umap@meta.data))
names(seurat_umap@meta.data) <- gsub("^.*pANN.*$","pANN", names(seurat_umap@meta.data))


if (params$remove_doublets == 'Y' ) {
  seurat_umap <- subset(seurat_umap, subset= DF!='Doublet')
  length(Idents(seurat_umap))
  saveRDS(seurat_umap, file=paste0(seurat_umap@project.name,'.seurat_dfx.rds'))
} else {
  saveRDS(seurat_umap, file=paste0(seurat_umap@project.name,'.seurat_df.rds'))
}


#Find cluster markers from RNA slot
DefaultAssay(seurat_umap) <- "RNA"
all.markers <- FindAllMarkers(seurat_umap, assay="RNA", min.pct=0.25, logfc.threshold=0.25)
top20.markers <- all.markers %>% group_by(cluster) %>% top_n(n=20, wt=avg_log2FC)

write.csv(all.markers, file=paste0(seurat_umap@project.name, '.all_markers.csv'))
write.csv(top20.markers, file=paste0(seurat_umap@project.name, '.top20_markers.csv'))

top5.markers <- top20.markers %>% group_by(cluster) %>% top_n(n=5, wt=avg_log2FC)
pdf(paste0(seurat_umap@project.name, '.top5_markers.pdf'))
DoHeatmap(seurat_umap, features=as.vector(top5.markers$gene), label=F) + theme(axis.text.y = element_text(size=5))
dev.off()

cell_markers<- c("EPCAM", "TFF3", "ALB", "DCN", "ACTA2", "VWF", "CD14", "CD68", "CD3D", "CD8A","NKG7","CD79A","CD19")
immune_markers <- c("CD3D", "CD4","IL7R", "FOXP3", "IL2RA","CD8A", "CD8B","GNLY", "NKG7","TRDC", "SELL", "CCR7", "CCR6", "KLRB1", "CD19", "MS4A1", "SDC1")

pdf("cell_markers.pdf")
DoHeatmap(seurat_umap, features=cell_markers, label=F) + scale_fill_viridis(option="D")
dev.off()

pdf("immune_markers.pdf")
DoHeatmap(seurat_umap, features=immune_markers, label=F) + scale_fill_viridis(option="D")
dev.off()



