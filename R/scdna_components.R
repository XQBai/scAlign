utils::globalVariables(c(
  "chr",
  "variable",
  "cell_index",
  "dropout_proportion", 
  "value"
))
#' Create Seurat object for scDNA-seq data and perform dimensionality reduction
#' @importFrom Seurat CreateSeuratObject SetAssayData GetAssayData
#' @importFrom Seurat ScaleData RunPCA RunUMAP FindNeighbors
#' @importFrom magrittr %>%
#' @param scdna_matrix_merge Merged scDNA-seq matrix (data.frame or matrix)
#' @param scdna_matrix_locs Data frame of segment locations (not used in this function, but kept for interface consistency)
#' @param dims_reduce Number of dimensions for reduction, default 50
#' @param output_file Output file name for saving Seurat object, default 'seurat_scDNA.rds'
#' @param output_dir Directory to save output files, default is current working directory
#' @return Seurat object after PCA, UMAP, and neighbor finding
#' @export
create_seurat_scdna <- function(scdna_matrix_merge,
                               scdna_matrix_locs,
                               dims_reduce = 50,
                               output_file = 'seurat_scDNA.rds',
                               output_dir = ".") {

  # Create Seurat object from input matrix
  seurat_obj <- CreateSeuratObject(
    counts = as.matrix(scdna_matrix_merge),
    project = 'scDNA',
    min.cells = -Inf,
    min.features = -Inf
  )

  seurat_obj <- SetAssayData(
    object = seurat_obj,
    layer= "data",
    new.data = GetAssayData(seurat_obj, layer = "counts")
  )

  # Perform dimensionality reduction pipeline
  seurat_obj <- seurat_obj %>%
    ScaleData(do.center = FALSE, do.scale = FALSE) %>%
    RunPCA(verbose = TRUE, features = rownames(seurat_obj)) %>%
    RunUMAP(
      umap.method = 'uwot',
      n.neighbors = min(20, ncol(seurat_obj)),
      dims = 1:dims_reduce
    ) %>%
    FindNeighbors(reduction = 'umap', dims = 1:2)

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save Seurat object
  saveRDS(seurat_obj, file = file.path(output_dir, output_file))
  return(seurat_obj)
}

#' Identify cellranger noise cells based on cellranger-dna metrics
#' @importFrom utils read.csv
#' @param per_cell_metrics_file Path to per_cell_summary_metrics.csv file
#' @param seurat_obj Seurat object or path to seurat_scDNA.rds file
#' @return Seurat object with celltype annotation
#' @export
identify_cellranger_noise <- function(per_cell_metrics_file = 'per_cell_summary_metrics.csv',
                                     seurat_obj = 'seurat_scDNA.rds') {

  # Load data
  per_cell_metrics <- read.csv(per_cell_metrics_file, sep = ",", header = TRUE)

  # Handle seurat_obj input (can be file path or Seurat object)
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }

  # Initialize celltype vector
  celltype <- rep(0, ncol(seurat_obj))

  # Mark noisy cells based on cellranger metrics
  noisy_indices <- which(per_cell_metrics$is_noisy == 1)
  celltype[noisy_indices] <- 'cellranger noise'

  # Add celltype annotation to Seurat object
  seurat_obj$celltype <- celltype

  return(seurat_obj)
}

#' Identify technical noise cells based on CNV dropout proportion
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab scale_y_continuous geom_hline theme element_text ggsave
#' @param seurat_obj Seurat object with scDNA-seq data
#' @param dropout_threshold Threshold for CNV dropout proportion, default 0.1
#' @param output_plot Whether to save dropout plot, default TRUE
#' @param plot_filename Output plot filename, default 'Miss_val_proportion.pdf'
#' @param output_dir Directory to save output files, default is current working directory
#' @return Seurat object with updated celltype annotation
#' @export
identify_technical_noise <- function(seurat_obj = 'seurat_scDNA.rds',
                                   dropout_threshold = 0.1,
                                   output_plot = TRUE,
                                   plot_filename = 'Miss_val_proportion.pdf',
                                   output_dir = ".") {

  # Handle seurat_obj input (can be file path or Seurat object)
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }

  # Extract CNV matrix
  cnv_matrix <- seurat_obj@assays$RNA@layers$counts

  # Calculate dropout proportion for each cell
  num_zeros <- apply(cnv_matrix, 2, function(x) sum(x == 0))
  dropout_prop <- num_zeros / nrow(cnv_matrix)

  # Create dropout proportion plot
  if (output_plot) {
    dropout_df <- data.frame(
      cell_index = seq_along(dropout_prop),
      dropout_proportion = dropout_prop
    )

    p <- ggplot(data = dropout_df, aes(x = cell_index, y = dropout_proportion)) +
      geom_point(size = 1, alpha = 0.5) +
      xlab("Cell index") +
      ylab("CNV dropout proportion") +
      scale_y_continuous(breaks = seq(0, max(dropout_df$dropout_proportion), 0.1)) +
      geom_hline(
        aes(yintercept = dropout_threshold, colour = '#990000', linetype = 'dashed'),
        show.legend = FALSE
      ) +
      theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
      theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
      theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
      theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) +
      theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold'))

    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    ggsave(filename = file.path(output_dir, plot_filename), plot = p, device = "pdf", width = 8, height = 6)
  }

  # Check if celltype column exists, if not, initialize as 0
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- rep(0, ncol(seurat_obj))
  }
  seurat_obj$celltype[which(dropout_prop > dropout_threshold)] <- 'noise'

  return(seurat_obj)
}


#' Identify normal cells based on chromosome-wise CNV analysis
#' @importFrom dplyr group_by summarise_all
#' @importFrom magrittr %>%
#' @importFrom reshape2 melt
#' @param seurat_obj Seurat object with celltype annotation
#' @param scdna_matrix_locs Data frame with chromosome information for segments
#' @param cnv_threshold Threshold for CNV value to identify tumor cells, default 2.5
#' @return Seurat object with updated celltype annotation
#' @export
identify_normal_cells <- function(seurat_obj='seurat_scDNA.rds',
                                 scdna_matrix_locs,
                                 cnv_threshold = 2.5) {

  # Handle seurat_obj input (can be file path or Seurat object)
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }

  # Check if celltype column exists, if not, initialize as 0
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- rep(0, ncol(seurat_obj))
  }

  # Filter out noise cells
  non_noise_indices <- which(seurat_obj$celltype == 0)

  if (length(non_noise_indices) == 0) {
    warning("No non-noise cells found")
    return(seurat_obj)
  }

  # Subset to non-noise cells
  seurat_subset <- subset(seurat_obj, cells = non_noise_indices)
  cnv_matrix <- as.data.frame(seurat_subset@assays$RNA@layers$counts)

  # Replace zeros with NA for analysis
  cnv_matrix[cnv_matrix == 0] <- NA
  cnv_matrix$chr <- scdna_matrix_locs$chr

  # Calculate chromosome-wise averages
  chr_averages <- cnv_matrix %>%
    group_by(chr) %>%
    summarise_all(list(mean), na.rm = TRUE) %>%
    as.data.frame()

  chr_averages$chr <- NULL

  # Melt data for analysis
  chr_averages_melt <- melt(chr_averages, measure.vars = colnames(chr_averages))

  # Classify cells based on CNV threshold
  cell_classification <- chr_averages_melt %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(is_tumor = any(value > cnv_threshold, na.rm = TRUE))

  # Update celltype annotation
  seurat_obj$celltype[non_noise_indices] <- ifelse(
    cell_classification$is_tumor, 'tumor', 'normal'
  )

  return(seurat_obj)
}

#' Identify replication cells (S phase vs G0G1) using mixture model analysis
#'
#' This function identifies S phase and G0G1 cells from tumor cells using a mixture model
#' approach based on Euclidean distance to normal CNV profile (CN=2).
#' @importFrom magrittr %>%
#' @importFrom stats dist
#' @importFrom grDevices pdf dev.off
#' @importFrom ggplot2 ggplot
#' @importFrom mixtools normalmixEM
#' @param seurat_obj Seurat object with celltype annotation, or path to RDS file
#' @param scdna_matrix_locs Data frame with chromosome information for segments
#' @param output_plot Whether to save mixture model plot, default TRUE
#' @param plot_filename Output plot filename, default 'mixture_cells.pdf'
#' @param refine_replicates Whether to further refine S phase/G0G1 classification using reidentify_replicates, default FALSE
#' @param mean_threshold Threshold for mean CNV value in reidentify_replicates, default 2.5
#' @param auto_threshold Whether to automatically determine threshold in reidentify_replicates, default 'manual'
#' @param output_dir Directory to save output files, default is current working directory
#' @return Seurat object with updated celltype annotation (S phase or G0G1)
#' @export
identify_replication_cells <- function(seurat_obj = 'seurat_scDNA.rds',
                                     scdna_matrix_locs,
                                     output_plot = TRUE,
                                     plot_filename = 'mixture_cells.pdf',
                                     refine_replicates = TRUE,
                                     mean_threshold = 2.5,
                                     auto_threshold = 'manual',
                                     output_dir = ".") {


  # Handle seurat_obj input (can be file path or Seurat object)
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }

  # Check if celltype column exists, if not, initialize as 0
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- rep('tumor', ncol(seurat_obj))
  }

  # Get tumor cells
  tumor_indices <- which(seurat_obj$celltype == 'tumor')

  if (length(tumor_indices) == 0) {
    warning("No tumor cells found in the Seurat object")
    return(seurat_obj)
  }

  # Subset to tumor cells
  seurat_tumor <- subset(seurat_obj, cells = tumor_indices)

  # Filter out chromosome X and Y segments
  locs_filtered <- scdna_matrix_locs %>% dplyr::filter(chr != 'chrX' & chr != 'chrY')

  # Extract CNV matrix for filtered segments
  cnv_matrix <- as.matrix(seurat_tumor@assays$RNA@layers$counts[1:nrow(locs_filtered), ])

  # Create normal CNV profile (CN=2) for comparison
  normal_cnv <- matrix(2, nrow = 1, ncol = nrow(cnv_matrix))

  # Calculate Euclidean distances from each cell to normal profile
  distances <- as.matrix(dist(t(cbind(t(as.matrix(normal_cnv)), cnv_matrix))))[1, ]
  distances <- distances[1:ncol(cnv_matrix) + 1]

  # Fit mixture model using EM algorithm
  mixture_model <- normalmixEM(distances)

  # Create mixture model plot if requested
  if (output_plot) {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    pdf(file.path(output_dir, plot_filename))
    mixtools::plot.mixEM(
      mixture_model,
      whichplots = 2,
      main2 = 'S phase cells',
      xlab2 = 'The Euclidean distance to normal cells',
      lwd2 = 3,
      marginal = TRUE
    )
    dev.off()
  }

  # Determine S phase component based on mixture model parameters
  mu <- mixture_model$mu
  lambda <- mixture_model$lambda
  sigma <- mixture_model$sigma

  # Get posterior probabilities for each cell
  posterior_labels <- apply(mixture_model$posterior, 1, function(x) which(x == max(x)))

  # Identify S phase component: either largest mu or largest sigma
  if (max(mu) > 2 * min(mu)) {
    s_phase_component <- which(mu == max(mu))
  } else {
    s_phase_component <- which(sigma == max(sigma))
  }

  # Identify S phase cells
  s_phase_cells <- which(posterior_labels == s_phase_component)

  # Update labels
  new_labels <- rep('G0G1', length(posterior_labels))
  new_labels[s_phase_cells] <- 'S phase'

  # Update celltype annotation in original Seurat object
  seurat_obj$celltype[tumor_indices] <- new_labels

  # Further refine S phase/G0G1 classification if requested
  if (refine_replicates) {
    message("Refining S phase/G0G1 classification using reidentify_replicates...")
    seurat_obj <- reidentify_replicates(
      seurat_obj = seurat_obj,
      mean_threshold = mean_threshold,
      auto_threshold = auto_threshold
    )
  }

  return(seurat_obj)
}

#' Re-identify replicates from cellranger noise and S phase cells
#' @param seurat_obj Seurat object with celltype annotation, or path to RDS file
#' @param mean_threshold Threshold for mean CNV value to distinguish S phase from noise, default 2.5
#' @param auto_threshold Whether to automatically determine threshold from data, default 'manual'
#' @return Seurat object with updated celltype annotation
#' @export
reidentify_replicates <- function(seurat_obj='seurat_scDNA.rds',
                                 mean_threshold = 2.5,
                                 auto_threshold = 'mannual') {

  # Handle seurat_obj input (can be file path or Seurat object)
  if (is.character(seurat_obj)) {
    seurat_obj <- readRDS(seurat_obj)
  }

  # Check if celltype column exists, if not, initialize as 0
  if (!"celltype" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj$celltype <- rep(0, ncol(seurat_obj))
  }

  # Get indices of cellranger noise and S phase cells
  noise_s_indices <- which(seurat_obj$celltype %in% c('cellranger noise', 'S phase'))

  if (length(noise_s_indices) == 0) {
    warning("No cellranger noise or S phase cells found")
    return(seurat_obj)
  }

  # Subset Seurat object to noise and S phase cells
  noise_subset <- subset(seurat_obj, cells = noise_s_indices)

  noise_matrix <- noise_subset@assays$RNA@layers$counts

  # Calculate mean CNV for each cell
  noise_cell_means <- apply(noise_matrix, 2, mean)

  # Determine threshold based on auto_threshold parameter
  if (auto_threshold == 'auto') {
    # Use median of cell_means as automatic threshold
    non_noise_subset <- subset(seurat_obj, cells = which(seurat_obj$celltype == 'G0G1'))
    non_noise_matrix <- non_noise_subset@assays$RNA@layers$counts
    # Calculate mean CNV for each cell
    cell_means <- apply(non_noise_matrix, 2, mean)
    threshold <- median(cell_means, na.rm = TRUE)
    # message(paste("Using automatic threshold:", round(threshold, 3)))
  } else {
    # Use provided mean_threshold
    threshold <- mean_threshold
    message(paste("Using manual threshold:", threshold))
  }

  # Re-identify based on determined threshold
  new_labels <- rep('noise', length(noise_s_indices))
  new_labels[which(noise_cell_means > threshold)] <- 'S phase'

  # Update celltype annotation
  seurat_obj$celltype[noise_s_indices] <- new_labels

  return(seurat_obj)
}

