#' Construct subclones in G0G1 tumor cells using Seurat clustering
#'
#' This function identifies subclones in G0G1 tumor cells by running Seurat clustering and assigning subclone labels.
#'
#' @param seurat_obj A Seurat object with celltype annotation
#' @param dims.reduce Number of dimensions for reduction (default 50)
#' @param resolution Clustering resolution parameter (default: 0.1 if <=600 cells, else 0.5)
#' @param merge_clusters Whether to merge overclustered subclones using Jaccard similarity (default TRUE)
#' @param jaccard_threshold Jaccard similarity threshold for merging clusters (default 0.95)
#' @param binarize_threshold Numeric threshold for binarizing CNV values (default 3)
#' @param output_dir Directory to save output files, default is current working directory
#' @return A Seurat object with subclone labels assigned
#' @export
construct_subclones <- function(seurat_obj, dims.reduce = 50, resolution = 0.5, merge_clusters = TRUE, jaccard_threshold = 0.99, binarize_threshold = 2, output_dir = ".") {

  celltype <- seurat_obj$celltype
  if (any(celltype == 'G0G1')){
    seurat_obj_subset <- subset(seurat_obj, cells = which(celltype == 'G0G1'))
  }else{
    seurat_obj_subset <- seurat_obj
  }
  seurat_obj_subset <- run_seurat(seurat_obj_subset, dims.reduce = dims.reduce, resolution = resolution)
    seurat_obj_subset$subclones <- paste0('C', as.numeric(seurat_obj_subset$seurat_clusters))
  # merge the overclustering to the same subclone
  if (merge_clusters) {
    message('Merging overclustered subclones using Jaccard similarity...')
    labels <- paste0('C', as.numeric(seurat_obj_subset$seurat_clusters))
    mat <- seurat_obj_subset@assays$RNA@layers$counts
    mat <- as.matrix(mat)
    res <- merge_clusters_by_jaccard(t(mat), labels, threshold = jaccard_threshold, binarize_threshold = binarize_threshold)
    labels_merge <- res$new_labels
    seurat_obj_subset$subclones <- labels_merge
  }

  # update the subclone in Seurat meta data (only if 'G0G1' exists)
  celltype <- seurat_obj$celltype
  if (any(celltype == 'G0G1')) {
    celltype[which(celltype == 'G0G1')] <- seurat_obj_subset$subclones
    seurat_obj$subclones <- celltype
  }
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(seurat_obj_subset, file = file.path(output_dir, 'seurat_scDNA_tumor.rds'))
  saveRDS(seurat_obj, file = file.path(output_dir, 'seurat_scDNA.rds'))
  return(seurat_obj)
}
#' Run Seurat clustering on scDNA-seq matrix
#'
#' This function performs feature selection, dimensionality reduction, and clustering on a scDNA-seq matrix.
#' @importFrom Seurat Idents Idents<-
#' @param seurat_obj A Seurat object
#' @param dims.reduce Number of dimensions for reduction (default 50)
#' @param resolution Clustering resolution parameter (default: 0.1 if <=600 cells, else 0.5)
#' @return A Seurat object with clustering results
#' @export
run_seurat <- function(seurat_obj, dims.reduce = 50, resolution = 0.5) {


  scdna_matrix = seurat_obj@assays$RNA@layers$counts
  # row.names(scdna_matrix) <- paste0('segment', 1:dim(scdna_matrix)[1])

  # Feature selection
  y <- tryCatch({
    Iden_signal_segments(scdna_matrix)
  }, warning = function(w) {
    return(0)
  }, error = function(e) {
    return(0)
  })

  if(length(y) == 1 && y == 0) {
    index_segment <- seq(1, dim(scdna_matrix)[1])
  } else {
    index_segment <- Iden_signal_segments(scdna_matrix)
  }

  signal_segment <- rownames(seurat_obj)[index_segment]
  seurat_obj <- subset(seurat_obj, features = signal_segment)

  seurat_obj <- seurat_obj %>%
    RunPCA(verbose = TRUE, features = rownames(seurat_obj), npcs = min(dims.reduce, ncol(seurat_obj) - 1)) %>%
    RunUMAP(umap.method = 'uwot', n.neighbors = min(20, ncol(seurat_obj)), dims = 1:min(dims.reduce, ncol(seurat_obj) - 1), verbose = FALSE) %>%
    FindNeighbors(reduction = 'umap', dims = 1:2) %>%
    FindClusters(reduction = 'umap', resolution = resolution)

  Idents(seurat_obj) <- paste0('C', as.numeric(Idents(seurat_obj)))
  return(seurat_obj)
}


#' Identify informative segments in scDNA-seq matrix using mixture model
#'
#' This function selects informative segments based on the standard deviation of CNVs across segments,
#' using a mixture model to distinguish signal from background.
#' @importFrom stats sd
#' @importFrom grDevices pdf dev.off
#' @importFrom mixtools normalmixEM
#' @param scdna_matrix A matrix of scDNA-seq CNV values (segments x cells)
#' @param output_dir Directory to save output files, default is current working directory
#' @return A vector of indices for informative segments
#' @export
Iden_signal_segments <- function(scdna_matrix, output_dir = ".") {
 
  arm_sd <- apply(scdna_matrix, 1, sd)
  arm_mean <- apply(scdna_matrix-2, 1, mean)
  arm_sd <- arm_sd[union(which(arm_sd != 0), which(arm_mean != 0))]

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  pdf(file.path(output_dir, 'Segments_density_curve.pdf'))
  mixtools::plot.mixEM(normalmixEM(arm_sd), whichplots = 2, main2 = 'Selection of informative segments', xlab2 = 'Standard deviation of CNVs on segments',
                      lwd2 = 3, marginal = TRUE)
  dev.off()

  estimate_mix_model <- normalmixEM(arm_sd)
  index_seg_cluster <- which(estimate_mix_model$lambda == max(estimate_mix_model$lambda))
  index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
  index_segment <- which(index_label != index_seg_cluster)

  index_segment_tmp <- index_segment
  index_seg_cluster_tmp <- index_seg_cluster

  # Stop select when length of selected segments is less than 20% of total
  while(length(index_segment) <= round(0.2 * dim(scdna_matrix)[1])){
    if (length(index_segment) >= round(0.21 * dim(scdna_matrix)[1])){
      break
    }
    segment_index_remain <- which(index_label == index_seg_cluster_tmp)
    arm_sd_tmp <- arm_sd[segment_index_remain]
    estimate_mix_model_tmp <- normalmixEM(arm_sd_tmp)
    index_seg_cluster_tmp <- which(estimate_mix_model_tmp$lambda == max(estimate_mix_model_tmp$lambda))
    index_label_tmp <- apply(estimate_mix_model_tmp$posterior, 1, function(x){which(x == max(x))})
    index_segment_tmp <- segment_index_remain[which(index_label_tmp != index_seg_cluster_tmp)]
    index_segment <- c(index_segment, index_segment_tmp)
  }
  return(index_segment)
}

#' Merge clusters by Jaccard similarity of CNV profiles
#'
#' This function merges clusters based on the Jaccard similarity of their binarized CNV profiles.
#' Clusters with Jaccard similarity above the threshold are merged into the same group.
#' @importFrom Seurat RenameCells
#' @importFrom stats setNames
#' @importFrom igraph graph_from_adjacency_matrix components
#' @param data A matrix or data.frame of CNV values (cells x features)
#' @param labels A vector or factor of cluster labels for each cell
#' @param threshold Jaccard similarity threshold for merging clusters (default 0.95)
#' @param binarize_threshold Numeric threshold for binarizing CNV values (default 3)
#' @return A list containing the Jaccard matrix, cluster graph, new labels, and cluster map
#' @export
merge_clusters_by_jaccard <- function(data, labels, threshold = 0.95, binarize_threshold = 3) {

  labels <- as.factor(labels)
  clusters <- levels(labels)

  # average CNV of clusters
  avg_mat <- sapply(sort(unique(labels)), function(rl) {
    rows <- which(labels == rl)
    round(colMeans(data[rows, , drop = FALSE]))
  })

  colnames(avg_mat) <- paste0("C", clusters)

  if(ncol(avg_mat) > 1){
    # Binarize CNV matrix for Jaccard calculation
    binary_mat <- apply(avg_mat, MARGIN = 1:2, function(x) {
      if (x > binarize_threshold) {
        return("amp")
      } else if (x == binarize_threshold) {
        return("natural")
      } else {
        return("del")
      }
    })

    n <- ncol(binary_mat)
    jaccard_mat <- matrix(0, n, n)
    rownames(jaccard_mat) <- colnames(jaccard_mat) <- colnames(avg_mat)

    # pairwise Jaccard index
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        inter <- sum(binary_mat[, i] == binary_mat[, j])
        union <- length(binary_mat[, i])
        jaccard_mat[i, j] <- inter / union
        jaccard_mat[j, i] <- jaccard_mat[i, j]
      }
    }
  }else{
    jaccard_mat =  matrix(1, nrow = 1, ncol = 1)
  }

  # Build a graph where edge = high similarity, and find connected components
  adj <- jaccard_mat > threshold

  if(ncol(adj) > 1){
    diag(adj) <- FALSE
  }else{
    diag(adj) <- TRUE
  }

  g <- graph_from_adjacency_matrix(adj, mode = "undirected")
  comps <- components(g)$membership

  # Map old cluster labels (C1, C2, ...) to new cluster groupings (G1, G2, ...)
  cluster_map <- paste0("C", comps)
  original_labels <- colnames(avg_mat)
  new_cluster_labels <- setNames(cluster_map, original_labels)

  # Assign new cluster labels for each original observation
  labels_new <- new_cluster_labels[paste0("C", labels)]
  names(labels_new) <- names(labels)

  # Return results
  return(list(
    jaccard_matrix = jaccard_mat,
    cluster_graph = g,
    new_labels = labels_new,
    cluster_map = new_cluster_labels
  ))
}


#' Generate subclone-level and tumor-level gene CNV matrices from scDNA-seq data
#'
#' This function generates a gene-level CNV matrix for all tumor cells and a subclone-averaged CNV matrix
#' based on a Seurat object annotated with subclones and cell types. It supports flexible input of gene-bin mapping,
#' raw scDNA-seq matrix, and barcode files.
#' @importFrom dplyr across where
#' @importFrom magrittr %>%
#' @param seurat_obj Seurat object (or path to RDS file) containing scDNA-seq data, annotated with subclones and cell types
#' @param genes_bin  gene-bin mapping RData file (should contain a variable genes_bin with $bins, $entrezgene, $symbolgene)
#' @param scdna_tsv  raw scDNA-seq matrix file (e.g., .tsv or .mtx)
#' @param barcode_txt barcode file (e.g., .txt)
#' @param subclone_matrix_out save subclone-averaged gene CNV matrix (default 'scdna_gene_subclones.rds')
#' @param output_dir Directory to save output files, default is current working directory
#' @return A list with:
#'   \item{scdna_gene_matrix_tumor}{Gene-level CNV matrix for all tumor cells (genes x cells)}
#'   \item{scdna_gene_matrix_subclone}{Subclone-averaged gene CNV matrix (subclones x genes)}
#' @export
generate_subclone_cnv_gene_matrix <- function(
  seurat_obj,
  genes_bin,
  scdna_tsv,
  barcode_txt,
  # tumor_matrix_out = 'scdna_tumor_matrix_gene.rds',
  subclone_matrix_out = 'scdna_gene_subclones.rds',
  output_dir = "."
) {

  # Generate gene-level CNV matrix
  scdna_gene_matrix <- Convert_scDNA_to_gene_matrix(
    genes_bin = genes_bin,
    scdna_tsv = scdna_tsv,
    barcode_txt = barcode_txt, 
    output_dir = output_dir
  )

  # Check column name consistency between gene matrix and Seurat object
  gene_list <- rownames(scdna_gene_matrix)
  common_cells <- intersect(colnames(scdna_gene_matrix), colnames(seurat_obj))
  if (length(common_cells) == 0) {
    stop("No common cells between scDNA gene matrix and Seurat object. Please check barcodes and Seurat object!")
  }
  scdna_gene_matrix <- scdna_gene_matrix[, common_cells]
  seurat_obj <- subset(seurat_obj, cells = common_cells)

  # Extract G0G1 (non-noise) cells
  cell_nonoise <- which(seurat_obj$celltype %in% c('G0G1'))
  seurat_obj_subset <- subset(seurat_obj, cells = cell_nonoise)
  scdna_matrix_tumor_gene <- scdna_gene_matrix[, colnames(seurat_obj_subset)]
  rownames(scdna_matrix_tumor_gene) <- rownames(scdna_gene_matrix)
  # saveRDS(scdna_matrix_tumor_gene, file = tumor_matrix_out)

  # Calculate subclone-averaged CNV matrix
  scdna_matrix_subset <- as.data.frame(t(scdna_gene_matrix[, cell_nonoise]))
  scdna_matrix_subset$subclone <- seurat_obj_subset$subclones
  CNV_subset_mean <- scdna_matrix_subset %>%
    dplyr::group_by(subclone) %>%
    dplyr::summarise(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
  subclone <- CNV_subset_mean$subclone
  CNV_subset_mean$subclone <- NULL
  CNV_subset_mean <- round(CNV_subset_mean)
  colnames(CNV_subset_mean) <- gene_list
  rownames(CNV_subset_mean) <- subclone

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(CNV_subset_mean, file = file.path(output_dir, subclone_matrix_out))

  return(CNV_subset_mean)
}

#' Convert scDNA-seq folder data to gene-level CNV matrix
#'
#' This function loads genome bin-gene mapping, DNA matrix, and barcode files from folders,
#' merges bins to gene-level, and outputs a gene-level CNV matrix suitable for downstream analysis.
#' @importFrom utils read.csv read.table write.csv head
#' @importFrom data.table fread
#' @importFrom dplyr filter %>%
#' @importFrom stats median
#' @param genes_bin Path to gene_bin.RData (should contain a variable genes_bin with $bins)
#' @param scdna_tsv Folder path containing DNA matrix (e.g., .tsv or .mtx file)
#' @param barcode_txt Folder path containing barcode file (e.g., .txt file)
#' @param output_dir Directory to save output files, default is current working directory
#' @return A list with gene-level CNV matrix and used gene-bin mapping
#' @export 
Convert_scDNA_to_gene_matrix <- function(
  genes_bin,
  scdna_tsv,
  barcode_txt,
  output_dir = "."
) {

  # Read genome reference
  load(genes_bin)
  # Read scDNA-seq matrix
  scdna_matrix <- fread(scdna_tsv)
  # Read barcode
  barcodes <- read.csv(barcode_txt, sep='\n', header = FALSE)
  colnames(scdna_matrix) <- as.character(barcodes$V1)

  # Merge bins to gene-level matrix
  # (Each gene may correspond to multiple bins)
  tmp_sizes <- lapply(genes_bin$bins, length)
  bin_index <- rep(1:length(tmp_sizes), unlist(tmp_sizes))
  scdna_matrix_bins <- as.data.frame(scdna_matrix[unlist(genes_bin$bins), ])
  scdna_matrix_bins$index <- bin_index

  scdna_gene_matrix <- scdna_matrix_bins %>%
    dplyr::group_by(scdna_matrix_bins$index) %>%
    dplyr::summarise_all(list(mean), na.rm = TRUE)

  scdna_gene_matrix$index <- NULL

  ### remove the duplicated genes
  scdna_gene_matrix$entrezgene <- genes_bin$entrezgene
  scdna_gene_matrix$symbolgene <- genes_bin$symbolgene

  #remove duplicated symbolgenes
  scdna_gene_matrix <-scdna_gene_matrix %>% filter(!duplicated(symbolgene))

  symbolgene <-scdna_gene_matrix$symbolgene

  scdna_gene_matrix$entrezgene <- NULL
  scdna_gene_matrix$symbolgene <- NULL

  rownames(scdna_gene_matrix) <- symbolgene
  colnames(scdna_gene_matrix) <- barcodes$V1

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(scdna_gene_matrix, file = file.path(output_dir, 'scdna_gene_matrix.rds'))
  return(scdna_gene_matrix)
}

