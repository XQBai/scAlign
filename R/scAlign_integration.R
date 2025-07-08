
utils::globalVariables("V2")
#' Main scAlign pipeline for scDNA and scRNA integration
#' @param scrna_normalized_matrix Normalized scRNA matrix
#' @param scdna_gene_matrix scDNA gene matrix
#' @param scdna_subclones scDNA subclones matrix
#' @param seurat_scdna scDNA Seurat object
#' @param seurat_epi epithelial Seurat object
#' @param gene_locs_path Path to gene locations bed file
#' @param signal_chr Numeric vector of signal chromosomes (e.g., c(1, 2, 4))
#' @param output_seurat_rds Path to save updated Seurat object (default 'seurat_epi.rds')
#' @param output_alignment_csv Path to save alignment results (default 'project_alignment.csv')
#' @param method Distance method: "euclidean" or "correlation" (default "euclidean")
#' @param output_dir Directory to save output files, default is current working directory
#' @return Updated Seurat object with alignment labels
#' @export
run_scalign_pipeline <- function(
    scrna_normalized_matrix,
    scdna_gene_matrix,
    scdna_subclones,
    seurat_scdna,
    seurat_epi,
    gene_locs_path,
    signal_chr = c(1, 2, 3, 4, 5),
    output_seurat_rds = 'seurat_epi.rds',
    output_alignment_csv = 'project_alignment.csv',
    method = "euclidean",
    output_dir = "."
) {

  # Get signal chromosomes
  gene_locs <- read.table(gene_locs_path)

  # Process signal chromosomes - check if they need 'chr' prefix
  if (all(signal_chr %in% gene_locs$V2)) {
    signal_chr_formatted <- signal_chr
  } else {
    signal_chr_formatted <- paste0('chr', signal_chr)
  }

  gene_locs_sub <- gene_locs %>% dplyr::filter(V2 %in% signal_chr_formatted)

  # Find common genes between scRNA and signal chromosomes
  common_genes <- intersect(rownames(scrna_normalized_matrix), gene_locs_sub$V1)
  common_genes <- intersect(common_genes, colnames(scdna_subclones))

  # Check if G0G1 cells exist in scDNA data
  if ("celltype" %in% names(seurat_scdna@meta.data) && any(seurat_scdna$celltype == 'G0G1')) {
    subclone_label <- colnames(seurat_scdna)[which(seurat_scdna$celltype == 'G0G1')]
    scdna_gene_matrix_tumor <- scdna_gene_matrix[, subclone_label]
    rownames(scdna_gene_matrix_tumor) <- rownames(scdna_gene_matrix)
  } else {
    # If no G0G1 cells or celltype column doesn't exist, use all cells
    message("No G0G1 cells found or celltype column missing, using all scDNA cells")
    scdna_gene_matrix_tumor <- scdna_gene_matrix
  }

  # Align matrices to common gene set
  aligned_data <- align_scdna_scrna_genes(
    scdna_gene_matrix_tumor,
    scrna_normalized_matrix[common_genes, ],
    scdna_subclones
  )
  scrna_common <- aligned_data[[1]]
  scdna_common <- aligned_data[[2]]
  scdna_subclones_common <- aligned_data[[3]]

  # Project scRNA to scDNA space
  scrna_proj_matrix <- project_scrna_to_scdna(scrna_common, scdna_common)

  # Find common genes for subclone assignment
  common_subclone_genes <- intersect(colnames(scdna_subclones_common), gene_locs_sub$V1)
  colnames(scrna_proj_matrix) <- colnames(scdna_subclones_common)
  rownames(scdna_subclones_common) <- rownames(scdna_subclones)

  # Assign scRNA cells to scDNA subclones
  assign_clone <- assign_scrna_subclones(
    scdna_subclones_common[, common_subclone_genes],
    scrna_proj_matrix[, common_subclone_genes],
    method = method
  )

  # Update Seurat object with alignment labels
  if ("condition" %in% names(seu_epi@meta.data)) {
    label_new <- seu_epi$condition
  } else {
    # If neither condition nor sample exists, create a default condition
    seu_epi$condition <- "tumor"
    label_new <- seu_epi$condition
  }

  label_new[which(label_new == 'tumor')] <- assign_clone

  seu_epi$pro_alignment <- label_new

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save results
  saveRDS(seu_epi, file = file.path(output_dir, output_seurat_rds))

  df <- as.data.frame(table(label_new))
  write.csv(df, file = file.path(output_dir, output_alignment_csv))

  return(seu_epi)
}

#' Align scDNA and scRNA matrices to common gene set
#' @param scdna_gene_matrix scDNA gene matrix
#' @param scrna_normalized_matrix Normalized scRNA matrix
#' @param scdna_gene_subclone scDNA clone matrix
#' @return List of aligned matrices (scrna_normalized_matrix, scdna_gene_matrix, scdna_gene_subclone)
#' @export
align_scdna_scrna_genes <- function(scdna_gene_matrix, scrna_normalized_matrix, scdna_gene_subclone){

  if (identical(rownames(scdna_gene_matrix), rownames(scrna_normalized_matrix)) && identical(rownames(scdna_gene_matrix), colnames(scdna_gene_subclone))){
    message("The genes are consistent between scRNA and scDNA")
  } else{
    mean_gene <- apply(scrna_normalized_matrix, 1, sum)
    scrna_normalized_matrix <- scrna_normalized_matrix[which(mean_gene != 0), ]
    common <- intersect(rownames(scdna_gene_matrix), rownames(scrna_normalized_matrix))
    scrna_normalized_matrix <- scrna_normalized_matrix[common, ]
    scdna_gene_matrix <- scdna_gene_matrix[common, ]
    scdna_gene_subclone <- as.matrix(scdna_gene_subclone[, common])
  }
  return(list(scrna_normalized_matrix, scdna_gene_matrix, scdna_gene_subclone))
}


#' Project scRNA expression into scDNA space using linear mapping
#' @param scrna_normalized_matrix scRNA expression matrix (genes x cells)
#' @param scdna_gene_matrix scDNA CNV matrix (genes x cells)
#' @return Projected scRNA matrix in scDNA space
#' @export
project_scrna_to_scdna <- function(scrna_normalized_matrix, scdna_gene_matrix){

  beta <- matrix(0, ncol = 1, nrow = dim(scrna_normalized_matrix)[1])
  X <- t(scrna_normalized_matrix)
  Y <- t(scdna_gene_matrix)

  for (i in 1:dim(beta)[1]){
    beta[i] <- (sum(X[, i])*sum(Y[, i]))/dim(Y)[1]/(norm(X[, i], type = '2')^2)
  }

  ### Projection space
  scrna_proj_matrix <- matrix(NA, nrow = dim(X)[1], ncol = dim(X)[2])

  for (i in 1:dim(scrna_proj_matrix)[1]){
    scrna_proj_matrix[i, ] <- X[i, ] * beta
  }
  return(scrna_proj_matrix)
}


#' Assign scRNA cells to scDNA subclones using distance or correlation
#' @importFrom stats dist cor
#' @param scdna_gene_subclone scDNA subclone matrix (subclones x genes)
#' @param scrna_proj_matrix Projected scRNA matrix in scDNA space (cells x genes)
#' @param method Distance method: "euclidean" or "correlation" (default "euclidean")
#' @return Named vector of subclone assignments for each scRNA cell
#' @export
assign_scrna_subclones <- function(scdna_gene_subclone, scrna_proj_matrix, method = "euclidean"){

  if (method == "euclidean") {
    dis.mat <- apply(scdna_gene_subclone, 1, function(col1) {
      apply(scrna_proj_matrix, 1, function(col2) {
        dist(rbind(col1,col2), method = "Euclidean")
      })
    })
    assign_clone <- unlist(apply(dis.mat, 1, function(x){which(x == min(x))}))
  } else if (method == "correlation") {
    dis.mat <- apply(scdna_gene_subclone, 1, function(col1) {
      apply(scrna_proj_matrix, 1, function(col2) {
        cor(col1, col2, method = "pearson")
      })
    })
    assign_clone <- unlist(apply(dis.mat, 1, function(x){which(x == max(x))}))
  } else {
    stop("Method must be either 'euclidean' or 'correlation'")
  }

  lookup_table <- sort(rownames(scdna_gene_subclone))
  names(lookup_table) <- as.character(seq(1, length(lookup_table), 1))

  # Replace numeric values with characters using the lookup table
  assign_clone <- lookup_table[as.character(assign_clone)]
  names(assign_clone) <- colnames(scrna_proj_matrix)

  return(assign_clone)
}

