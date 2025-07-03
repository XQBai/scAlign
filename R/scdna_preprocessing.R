#' Preprocessing function for scDNA-seq data
#' @importFrom data.table fread
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by mutate filter
#' @param genome_reference_file Path to genome reference bed file
#' @param scdna_tsv Path to scDNA-seq matrix tsv file
#' @param barcode_txt Path to barcode txt file
#' @param bin_size Number of bins to merge, default 50
#' @param dims_reduce Number of dimensions for reduction, default 50
#' @return list(seurat_obj, matrix_rds, locs_rds, seurat_rds)
#' @export
scdna_preprocessing <- function(genome_reference_file, scdna_tsv, barcode_txt, bin_size=50, dims_reduce=50) {

  # Read genome reference
  genome_reference <- fread(genome_reference_file, header=F, data.table=F)
  # Read scDNA-seq matrix
  scdna_matrix <- fread(scdna_tsv)
  # Read barcode
  barcodes <- read.csv(barcode_txt, sep='\n', header = FALSE)
  colnames(scdna_matrix) <- as.character(barcodes$V1)

  # Merge bins
  genome_reference <- genome_reference %>% group_by(V1) %>% mutate(bin=floor(row_number(V1)/bin_size))
  genome_reference$bin_corrected <- cumsum(c(0,as.numeric(diff(genome_reference$bin))!=0))

  scdna_matrix$chr <- genome_reference$V1
  scdna_matrix$bin <- genome_reference$bin_corrected

  # Merge CNV by median
  scdna_matrix_merge <- scdna_matrix %>% group_by(chr,bin) %>% summarise_all(list(median))
  scdna_matrix_locs <- scdna_matrix_merge[,c("chr","bin")]
  scdna_matrix_locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY')
  scdna_matrix_merge <- scdna_matrix_merge %>% filter(chr != 'chrX' & chr != 'chrY')
  scdna_matrix_merge$chr <- NULL
  scdna_matrix_merge$bin <- NULL
  row.names(scdna_matrix_merge) <- paste0('segment', 1:dim(scdna_matrix_merge)[1])


  # Save results
  saveRDS(scdna_matrix_merge, file='scdna_matrix_all_barcodes.rds')
  saveRDS(scdna_matrix_locs, file='scdna_matrix_locs.rds')

  return(list(scdna_matrix_merge,
              scdna_matrix_locs))
}

