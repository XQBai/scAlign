utils::globalVariables(c("id", "score", "score2", "score3", "score4", "me", "mn", ".", "cluster", "avg_log2FC"))
#' Run full scRNA-seq integration and analysis pipeline with step control
#' @param rds_files Vector of Seurat RDS file paths
#' @param sample_list List of sample names for each Seurat object (optional)
#' @param condition_list List of condition names for each Seurat object (optional)
#' @param do_merge Whether to merge Seurat objects (default TRUE)
#' @param do_process Whether to run SCTransform/PCA/UMAP/clustering (default TRUE)
#' @param do_marker Whether to find and plot markers (default TRUE)
#' @param do_epithelial Whether to extract epithelial cells and rerun pipeline (default TRUE)
#' @param output_dir Directory to save output files, default is current working directory
#' @param dims.reduce Number of PCA，default 50
#' @param cluster.res Clustering resoultion, default 0.8
#' @param ... Other parameters for downstream functions
#' @return List of all intermediate and final objects
#' @export
run_scrna_pipeline <- function(
    rds_files,
    sample_list = NULL,
    condition_list = NULL,
    do_merge = TRUE,
    do_process = TRUE,
    do_marker = TRUE,
    do_epithelial = TRUE,
    output_dir = ".",
    dims.reduce = 50,
    cluster.res = 0.8,
    ...
) {
  seu_merged <- NULL
  aggr <- NULL
  seu_epi <- NULL

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (length(rds_files) > 1 && do_merge) {
    seu_merged <- merge_seurat_objects(rds_files, sample_list, condition_list)
    saveRDS(seu_merged, file = file.path(output_dir, 'seurat_merged.rds'))
  } else if (length(rds_files) == 1) {
    seu_merged <- readRDS(rds_files)
  } else {
    stop("No input RDS files provided!")
  }

  if (do_process) {
    aggr <- process_merged_seurat(seu_merged, dims.reduce = dims.reduce, cluster.res = cluster.res)
  } else {
    aggr <- seu_merged
  }

  if (do_marker) {
    find_and_plot_markers(aggr, output_dir = output_dir)
  }

  if (do_epithelial) {
    seu_epi <- extract_epithelial_cells(aggr)
    if (!is.null(seu_epi) && "Phase" %in% names(seu_epi@meta.data)) {
      seu_epi <- subset(seu_epi, cells = which(seu_epi$Phase != 'S'))
    }
    saveRDS(seu_epi, file = file.path(output_dir, 'seurat_epi.rds'))
  }

  return(seu_epi=seu_epi)
}

#' Merge multiple Seurat objects and harmonize barcodes
#' @param rds_files Character vector of Seurat RDS file paths
#' @param sample_list List of sample names for each Seurat object (optional)
#' @param condition_list List of condition names for each Seurat object (optional)
#' @return Merged Seurat object
#' @export
merge_seurat_objects <- function(rds_files, sample_list=NULL, condition_list = NULL) {

  modifyBarcodeSuffix <- function(barcodes, cell.suffix) {
    bc_has_suffix = (length(grep("-", barcodes)) > 0)
    if (cell.suffix == 'N' && bc_has_suffix) {
      new_barcodes = gsub("-\\d", '', barcodes)
    } else if (cell.suffix %in% c('K','N')) {
      new_barcodes = barcodes
    } else if (bc_has_suffix) {
      new_bc = paste0('-',cell.suffix)
      new_barcodes = gsub("-\\d", new_bc, barcodes)
    } else {
      new_barcodes = paste(barcodes, paste0('-',cell.suffix), sep="")
    }
    return(new_barcodes)
  }
  seu <- vector("list", length(rds_files))

  for (i in seq_along(rds_files)) {
    seu[[i]] <- readRDS(rds_files[i])
    barcodes <- colnames(seu[[i]])
    seu[[i]] <- RenameCells(seu[[i]], new.names=modifyBarcodeSuffix(barcodes, i))

    # Add sample and condition metadata if provided
    if (!is.null(sample_list) && length(sample_list) >= i) {
      seu[[i]]$sample <- sample_list[[i]]
    }
    if (!is.null(condition_list) && length(condition_list) >= i) {
      seu[[i]]$condition <- condition_list[[i]]
    }
  }

  seu_merged <- merge(x=seu[[1]], y=seu[-1])
  return(seu_merged)
}

#' SCTransform, PCA, clustering, UMAP for merged Seurat object
#' @param seu_merged Merged Seurat object
#' @param dims.reduce Number of dimensions for reduction (default 20)
#' @param cluster.res Clustering resolution (default 1)
#' @return Processed Seurat object
#' @export
process_merged_seurat <- function(seu_merged, dims.reduce=20, cluster.res=1) {

  aggr <- seu_merged %>%
    SCTransform(return.only.var.genes=FALSE) %>%
    RunPCA(algorithm=4, approx=F) %>%
    FindNeighbors(dims=1:dims.reduce) %>%
    FindClusters(resolution=cluster.res) %>%
    BuildClusterTree(reorder=T, reorder.numeric=T, dims=1:dims.reduce) %>%
    RunUMAP(dims=1:dims.reduce)
  return(aggr)
}

#' Find markers and plot heatmaps for merged Seurat object
#' @importFrom Seurat FindAllMarkers
#' @importFrom utils write.csv
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by top_n
#' @importFrom grDevices pdf dev.off
#' @importFrom Seurat DoHeatmap
#' @importFrom viridis scale_fill_viridis
#' @importFrom ggplot2 theme element_text
#' @param aggr Processed Seurat object
#' @param output_dir Directory to save output files, default is current working directory
#' @return NULL
#' @export
find_and_plot_markers <- function(aggr, output_dir = ".") {

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  markers_RNA <- FindAllMarkers(aggr, assay = "RNA", min.pct = 0.25, logfc.threshold = 0.25)
  write.csv(markers_RNA, file = file.path(output_dir, "markers_RNA.csv"))
  top20 <- markers_RNA %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
  write.csv(top20, file = file.path(output_dir, "top20_markers_RNA.csv"))
  top5 <- markers_RNA %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
  pdf(file.path(output_dir, "top5_cluster.pdf"))
  DoHeatmap(aggr, features = top5$gene) + scale_fill_viridis(option = "D") + theme(axis.text.y = element_text(size = 6))
  dev.off()
}

#' Extract epithelial cells and rerun Seurat pipeline
#' @importFrom Seurat AddModuleScore
#' @importFrom dplyr filter arrange select %>%
#' @importFrom Seurat SCTransform RunPCA RunUMAP FindNeighbors FindClusters BuildClusterTree
#' @param seu_merged Processed Seurat object
#' @param ndims Number of dimensions for PCA/UMAP (default 100)
#' @return Epithelial Seurat object
#' @export
extract_epithelial_cells <- function(seu_merged, ndims=100) {
  # ..epi_genes, nonepi_genes, b_plasma, t_general
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

  seu_merged <- AddModuleScore(seu_merged, features = list(epi_genes), name = "epithelial")
  seu_merged <- AddModuleScore(seu_merged, features = list(nonepi_genes), name = "nonepithelial")
  seu_merged <- AddModuleScore(seu_merged, features = list(epi_genes), name = "epithelial")
  seu_merged <- AddModuleScore(seu_merged, features = list(nonepi_genes), name = "nonepithelial")
  seu_merged <- AddModuleScore(seu_merged, features = list(b_plasma), name='b_plasma')
  seu_merged <- AddModuleScore(seu_merged, features = list(t_general), name='t_cells')

  epi_scores <- data.frame(id=seu_merged@active.ident, score=seu_merged$epithelial1, score2=seu_merged$nonepithelial1,
                           score3=seu_merged$b_plasma1, score4=seu_merged$t_cells1)

  epi_scores_summary <- epi_scores %>% dplyr::group_by(id) %>%
    dplyr::summarize(me=mean(score), mn=mean(score2), mb=mean(score3), mt=mean(score4)) %>% as.data.frame()

  # write.csv(epi_scores_summary, 'module_scores.csv')

  # ne_cutoff is minimum positive non-epithelial module score per cluster, for clusters with epithelial module score < 0
  ne_cutoff <- epi_scores_summary %>% filter(me < 0) %>% filter(mn > 0) %>% arrange(mn) %>% head(1) %>% select(mn) %>% as.numeric()
  epi_clusters <- epi_scores_summary %>% filter(me > 0.01) %>% filter(mn < ne_cutoff) %>% select(id) %>% .$id %>% as.character()


  # ...epi_scores_summary, epi_clusters...
  seu_epi <- subset(seu_merged, idents=epi_clusters)
  seu_epi <- SCTransform(seu_epi, vars.to.regress='nCount_RNA', return.only.var.genes=FALSE, verbose=T)
  seu_epi <- RunPCA(seu_epi, verbose=T, npcs = ndims, approx=F)
  seu_epi <- RunUMAP(seu_epi, dims = 1:ndims, verbose=T)
  seu_epi <- FindNeighbors(seu_epi, dims=1:ndims, verbose=T)
  seu_epi <- FindClusters(seu_epi, algorithm = 1)
  seu_epi <- BuildClusterTree(seu_epi, reorder=T, reorder.numeric = T, dims=1:ndims)
  return(seu_epi)
}




