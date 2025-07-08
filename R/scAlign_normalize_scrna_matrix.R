#' Main scRNA normalization pipeline
#' @importFrom methods slotNames
#' @param seurat_epi Seurat RDS file
#' @param gene_locs_path Path to gene location file
#' @param output_rds Path to save normalized matrix
#' @param output_dir Directory to save output files, default is current working directory
#' @return Normalized matrix (invisible)
#' @export
run_scrna_normalization <- function(seurat_epi, gene_locs_path, output_rds = 'RNA_normalized_matrix.rds', output_dir = ".") {

  if ('condition' %in% names(seu_epi@meta.data)){
    seu_epi <- subset(seu_epi, cells = which(seu_epi$condition == 'tumor'))
  } else if ("Phase" %in% names(seu_epi@meta.data)){
    seu_epi <- subset(seu_epi, cells = which(seu_epi$Phase != 'S'))
  } else {
    seu_epi <- seu_epi
  }

  # Check if layers exist in Seurat object
  if ("layers" %in% slotNames(seu_epi@assays$RNA)) {
    raw_matrix <- seu_epi@assays$RNA@layers$counts
  } else {
    raw_matrix <- seu_epi@assays$RNA@counts
  }

  gene_locs <- fread(gene_locs_path)
  normalized_rna_matrix <- Normalize_RNA(raw_matrix, gene_locs)
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  saveRDS(normalized_rna_matrix, file = file.path(output_dir, output_rds))
  return(normalized_rna_matrix)
}


#' Filter genes below mean expression cutoff
#' @param data Expression matrix
#' @param cutoff Numeric, mean expression cutoff
#' @return Filtered matrix
#' @export
Filter_genes_below_mean_exp_cutoff <- function(data, cutoff = 0.1) {

  data <- as.matrix(data)
  average_gene <- rowMeans(data)
  remain_indices <- which(average_gene > cutoff)
  data[remain_indices, , drop = FALSE]
}

#' Filter genes expressed in fewer than min_num_cells
#' @param data Expression matrix
#' @param min_num_cells Integer, minimum number of cells
#' @return Filtered matrix
#' @export
Filter_genes <- function(data, min_num_cells = 3) {
  genes_exp_cn <- apply(data, 1, function(x) sum(x > 0 & !is.na(x)))
  genes_remain <- which(genes_exp_cn > min_num_cells)
  if (length(genes_remain) == 0) stop('All genes were removed.')
  data[genes_remain, , drop = FALSE]
}

#' Normalize counts by sequencing depth
#' @param data Expression matrix
#' @param normalize_factor Numeric, normalization factor (default: median of column sums)
#' @return Normalized matrix
#' @export
normalize_counts_by_seq_depth <- function(data, normalize_factor = NA) {
  cs <- colSums(data)
  data <- sweep(data, 2, cs, '/')
  if (is.na(normalize_factor)) normalize_factor <- median(cs)
  data * normalize_factor
}

#' Log2 transform
#' @param data Matrix
#' @return Log2 transformed matrix
#' @export
Logtransform <- function(data) log2(data + 1)

#' Invert log2 transform
#' @param data Matrix
#' @return Inverse log2 transformed matrix
#' @export
Invert_logtransform <- function(data) 2^data - 1

#' Subtract average reference expression
#' @param data Matrix
#' @return Centered matrix
#' @export
Subtract_ref_mean_exp <- function(data) {
  average_gene <- rowMeans(data)
  sweep(data, 1, average_gene, '-')
}

#' Apply max/min threshold bounds
#' @param data Matrix
#' @param threshold Numeric
#' @return Thresholded matrix
#' @export
Max_threshold_bounds <- function(data, threshold = 3) {
  data[data > threshold] <- threshold
  data[data < -threshold] <- -threshold
  data
}

#' Smoothing by chromosome using moving average
#' @importFrom magrittr %>%
#' @param data Matrix
#' @param window_length Integer
#' @param gene_locs Data.frame with columns: gene, chr, start, end
#' @return Smoothed matrix
#' @export
Smooth_by_chromosome <- function(data, window_length = 101, gene_locs) {
  names(gene_locs) <- c('gene', 'chr', 'start', 'end')
  gene_locs_remain <- gene_locs %>% dplyr::filter(gene_locs$gene %in% rownames(data))
  data <- data[gene_locs_remain$gene, , drop = FALSE]
  chrs <- unique(gene_locs_remain$chr)
  for (chr in chrs) {
    chr_genes_id <- which(gene_locs_remain$chr == chr)
    chr_data <- data[chr_genes_id, , drop = FALSE]
    if (nrow(chr_data) > 1) {
      smoothed_chr_data <- smooth_window(chr_data, window_length)
      data[chr_genes_id, ] <- smoothed_chr_data
    }
  }
  data
}

#' Center cells across chromosome
#' @importFrom stats median
#' @param data Matrix
#' @param method 'median' or 'mean'
#' @return Centered matrix
#' @export
Center_cells_across_chromosome <- function(data, method = 'median') {
  if (method == 'median') {
    row_median <- apply(data, 2, median, na.rm = TRUE)
    t(apply(data, 1, '-', row_median))
  } else {
    row_mean <- apply(data, 2, mean, na.rm = TRUE)
    t(apply(data, 1, '-', row_mean))
  }
}

#' Helper for smoothing window (internal)
#' @keywords internal
smooth_helper <- function(obs_data, window_length) {
  orig_obs_data = obs_data
  nas = is.na(obs_data)
  obs_data = obs_data[!nas]
  obs_length <- length(obs_data)
  end_data <- obs_data
  tail_length = (window_length - 1)/2
  if (obs_length >= window_length){
    nas_tmp = is.na(obs_data)
    vals = obs_data[! nas_tmp]
    custom_filter_denominator = ((window_length-1)/2)^2 + window_length
    custom_filter_numerator = c(seq_len((window_length-1)/2), ((window_length-1)/2)+1, c(((window_length-1)/2):1))
    custom_filter = custom_filter_numerator/rep(custom_filter_denominator, window_length)
    smoothed = stats::filter(vals, custom_filter, sides=2)
    ind = which(! is.na(smoothed))
    vals[ind] = smoothed[ind]
    obs_data[! nas_tmp] = vals
    end_data <- obs_data
  }
  obs_count <- length(obs_data)
  numerator_counts_vector = c(seq_len(tail_length), tail_length + 1, c(tail_length:1))
  iteration_range = ifelse(obs_count > window_length, tail_length, ceiling(obs_count/2))
  for (tail_end in seq_len(iteration_range)) {
    end_tail = obs_count - tail_end + 1
    d_left = tail_end - 1
    d_right = obs_count - tail_end
    d_right = ifelse(d_right > tail_length, tail_length, d_right)
    r_left = tail_length - d_left
    r_right = tail_length - d_right
    denominator = (((window_length - 1)/2)^2 + window_length) - ((r_left * (r_left + 1))/2) - ((r_right * (r_right + 1))/2)
    left_input_vector_chunk = obs_data[seq_len(tail_end + d_right)]
    right_input_vector_chunk = obs_data[(end_tail - d_right):obs_length]
    numerator_range = numerator_counts_vector[(tail_length + 1 - d_left):(tail_length + 1 + d_right)]
    end_data[tail_end] = sum(left_input_vector_chunk * numerator_range)/denominator
    end_data[end_tail] = sum(right_input_vector_chunk * rev(numerator_range))/denominator
  }
  orig_obs_data[!nas] = end_data
  orig_obs_data
}

#' Smoothing window for matrix
#' @param data Matrix
#' @param window_length Integer
#' @return Smoothed matrix
#' @export
smooth_window <- function(data, window_length) {
  if (window_length < 2) {
    message('window length < 2, returning original unmodified data')
    return(data)
  }
  data_sm <- apply(data, 2, smooth_helper, window_length=window_length)
  row.names(data_sm) <- row.names(data)
  colnames(data_sm) <- colnames(data)
  data_sm
}

#' Subset Seurat object by condition or Phase
#' @param seu Seurat object
#' @param condition_value Character, e.g. 'tumor'
#' @return Subsetted Seurat object
#' @export
subset_seurat_by_condition_or_phase <- function(seu, condition_value = 'tumor') {
  if ('condition' %in% names(seu@meta.data)) {
    seu <- subset(seu, cells = which(seu$condition == condition_value))
  } else if ("Phase" %in% names(seu@meta.data)) {
    seu <- subset(seu, cells = which(seu$Phase != 'S'))
  }
  seu
}

#' Full RNA normalization pipeline
#' @param scRNA_raw_matrix Raw expression matrix
#' @param gene_locs Data.frame with gene locations
#' @return Normalized matrix
#' @export
Normalize_RNA <- function(scRNA_raw_matrix, gene_locs) {
  tmp_mat <- Filter_genes_below_mean_exp_cutoff(scRNA_raw_matrix, 0.1)
  tmp_mat <- Filter_genes(tmp_mat, 3)
  mat_normalize <- normalize_counts_by_seq_depth(tmp_mat)
  mat_log <- Logtransform(mat_normalize)
  mat_subtract <- Subtract_ref_mean_exp(mat_log)
  mat_bounds <- Max_threshold_bounds(mat_subtract, 3)
  mat_smooth <- Smooth_by_chromosome(mat_bounds, window_length = 101, gene_locs)
  mat_center <- Center_cells_across_chromosome(mat_smooth, method = 'median')
  mat_subtract_adj <- Subtract_ref_mean_exp(mat_center)
  mat_exp <- Invert_logtransform(mat_subtract_adj) + 1
  return(mat_exp)
}

