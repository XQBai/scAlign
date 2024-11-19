library(dplyr)
library(Seurat)
library(SeuratObject, lib.loc="/venvs2/anaconda3/envs/scrna-4/lib/R/library")

Normalize_RNA <- function(scRNA_raw_matrix, gene_locs){
  
  #Step1
  mat <- Filter_genes_below_mean_exp_cutoff(scRNA_raw_matrix, 0.1)
  #Step2
  mat <- Filter_genes(mat, 3)
  #Step3
  mat_normalize <- normalize_counts_by_seq_depth(mat)
  #Step4
  mat_log <- Logtransform(mat_normalize)
  # Step5
  mat_subtract <- Subtract_ref_mean_exp(mat_log)
  # Step6 
  mat_bounds <- Max_threshold_bounds(mat_subtract, 3)
  # Step7: smooth data along chromosome with gene windows
  mat_smooth <- Smooth_by_chromosome(mat_bounds, window_length = 101, gene_locs)
  # Step8: Center cells by median
  mat_center <- Center_cells_across_chromosome(mat_smooth, method = 'median')
  # Step9: Adjustment to subtract the relative normal cells
  mat_subtract_adj <- Subtract_ref_mean_exp(mat_center)
  # Step10: Revert the log transformation 
  mat_exp <- Invert_logtransform(mat_subtract_adj) + 1
  return(mat_exp)
}

## Filter genes that expressed fewer than certain number of cells
Filter_genes_below_mean_exp_cutoff <- function(mat, cutoff){
  
  average_gene <- rowMeans(mat)
  remain_indices <- which(average_gene > cutoff)
  mat <- mat[remain_indices, ]
  
  return(mat)
}

Filter_genes <- function(mat, min_num_cells){
  
  total_genes <- dim(mat)[1]
  genes_exp_cn <- apply(mat, 1, function(x){sum(x>0 & !is.na(x))})
  genes_remain = which(genes_exp_cn > min_num_cells)
  
  if (length(genes_remain) > 0){
    sprintf("%d genes were remained for downstream analysis", length(genes_remain))
    if (length(genes_remain) == 0){
      stop('All genes were removed.')
    }
  
  }
  mat <- mat[genes_remain, ]
  return(mat)
}

## Normalizes count data by total sum scaling
normalize_counts_by_seq_depth <- function(mat, normalize_factor = NA){
  
  # Sum of reads count per cell
  cs = colSums(mat)
  
  # make fraction of total counts 
  data <- sweep(mat, STATS = cs, MARGIN = 2, FUN='/')
  
  if (is.na(normalize_factor)){
    normalize_factor = median(cs)
  }
  data <- data * normalize_factor
  return(data)
}

## Log transformation 
Logtransform <- function(data){
  data <- log2(data + 1)
  return(data)
}

## Invert Log transformation
Invert_logtransform <- function(data){
  data <- 2^data - 1
  return(data)
}

## Subtract average reference
Subtract_ref_mean_exp <- function(data, inv_log = FALSE){
  
  average_gene <- apply(data, 1, mean)
  data <- sweep(data, STATS = average_gene, MARGIN = 1, FUN='-')
  
  return(data)
}

# Apply max/min threshold bounds 
Max_threshold_bounds <- function(data, threshold){
  
  data[data > threshold] <- threshold
  data[data < (-1 * threshold) ] <- -1 * threshold
  
  return(data)
}

### Smoothing by chromosome through moving average
smooth_helper <- function(obs_data, window_length){
  
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
  
  # defining the iteration range in cases where the window size is larger than the number of genes. In that case we only iterate to the half since the process is applied from both ends.
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
  # replace the original values by smoothed values
  
  orig_obs_data[!nas] = end_data
  
  return(orig_obs_data)
}

smooth_window <- function(data, window_length){
  
  if (window_length < 2){
    print('window length < 2, returning orginal unmodified data')
    return(data)
  }
  
  # Fix ends that couldn't be smoothed 
  data_sm <- apply(data, 2, smooth_helper, window_length=window_length)
  # Set back row and columns names
  row.names(data_sm) <- row.names(data)
  colnames(data_sm) <- colnames(data)
  
  return(data_sm)
}


Smooth_by_chromosome <- function(data, window_length, gene_locs, smooth_end = TRUE){
  
  names(gene_locs) <- c('gene', 'chr', 'start', 'end')
  ## Select the genes in matrix
  gene_locs_remain <- gene_locs %>% dplyr::filter(gene %in% rownames(data))
  ## reorder the genes order in matrix 
  data <- data[gene_locs_remain$gene, ]
  
  chrs = unique(gene_locs_remain$chr)
  for (chr in chrs){
    chr_genes_id = which(gene_locs_remain$chr == chr)
    sprintf(paste0('smooth_by_chromosome: chr: ', chr))
    chr_data = data[chr_genes_id, ]
    
    # consider the only one gene case
    if (is.null(nrow(chr_data))){
      chr_data <- as.matrix(t(chr_data))
    }
    
    if(nrow(chr_data) > 1){
      smoothed_chr_data = smooth_window(chr_data, window_length)
      data[chr_genes_id, ] <- smoothed_chr_data
    }
  }
  return(data)
}

## Centering cells after smoothing 
Center_cells_across_chromosome <- function(data, method = 'median'){
  
  if (method == 'median'){
    row_median <- apply(data, 2, function(x){median(x, na.rm = TRUE)})
    data <- t(apply(data, 1, '-', row_median))
  }else if (method == 'mean'){
    row_mean <- apply(data, 2, function(x){mean(x, na.rm = TRUE)})
    data <- t(apply(data, 1, '-', row_mean))
  }
  return(data)
}

########################################################
# options(future.globals.maxSize = 8000 * 1024^2)
seu_epi <- readRDS('seurat_epi.rds')
gene_locs <- read.table('../gene_locs.sorted.bed')

if ('condition' %in% names(seu_epi@meta.data)){
  seu_epi <- subset(seu_epi, cells = which(seu_epi$condition == 'tumor'))
  #seu_epi <- subset(seu_epi, subset = condition == 'tumor')
  
} else if ("Phase" %in% names(seu_epi@meta.data)){
  seu_epi <- subset(seu_epi, cells = which(seu_epi$Phase != 'S'))
  #seu_epi <- subset(seu_epi, subset = Phase != 'S')
} else {
  seu_epi <- seu_epi
}

raw_matrix <- seu_epi@assays$RNA@counts
normalized_rna_matrix <- Normalize_RNA(raw_matrix, gene_locs)

saveRDS(normalized_rna_matrix, file = 'RNA_normalized_matrix.rds')











