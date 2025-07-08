utils::globalVariables(c("chr", "bin"))

#' Plot subclonal CNV heatmap for scDNA-seq data
#' This function generates a subclonal CNV heatmap from a Seurat object and segment location information.
#' It supports plotting with or without noise/normal cells, and saves the heatmap as a PDF.
#' @param seurat_obj A Seurat object with subclones and celltype annotation
#' @param scdna_matrix_locs A data.frame of segment locations (e.g., read from scdna_matrix_locs.rds)
#' @param output_pdf Output PDF file name (default 'Subclones_heatmap.pdf')
#' @param celltype Character vector of celltypes from the plot (e.g., c('normal','noise')), default NULL
#' @param max_cnv_value Maximum CNV value to cap in the heatmap (default 6)
#' @param hc_rds Optional path to hierarchical clustering RDS file for subclone ordering (default NULL)
#' @return Invisibly returns the matrix used for plotting
#' @export
plot_subclonal_heatmap <- function(
  seurat_obj,
  scdna_matrix_locs,
  output_pdf = 'Subclones_heatmap.pdf',
  celltype = "G0G1",
  max_cnv_value = 6,
  hc_rds = NULL
) {
  # Data preparation
  scdna_matrix_merge <- as.data.frame(seurat_obj@assays$RNA@layers$counts)
  colnames(scdna_matrix_merge) <- colnames(seurat_obj)
  scdna_matrix_merge$chr <- scdna_matrix_locs$chr
  scdna_matrix_merge <- scdna_matrix_merge %>% dplyr::filter(chr != 'chrX' & chr != 'chrY')
  scdna_matrix_merge$chr <- NULL
  scdna_matrix_locs <- scdna_matrix_locs %>% dplyr::filter(chr != 'chrX' & chr != 'chrY')

  sum_row <- apply(scdna_matrix_merge, 1, mean)
  row_index <- which(sum_row != 0)
  scdna_matrix_merge <- scdna_matrix_merge[row_index, ]
  scdna_matrix_locs <- scdna_matrix_locs[row_index, ]

  scdna_matrix_merge$chr <- scdna_matrix_locs$chr
  scdna_matrix_merge$bin <- scdna_matrix_locs$bin
  scdna_matrix_merge <- scdna_matrix_merge %>% dplyr::arrange(bin)
  scdna_matrix_locs <- scdna_matrix_locs %>% dplyr::arrange(bin)
  scdna_matrix_merge$chr <- NULL
  scdna_matrix_merge$bin <- NULL

  scdna_matrix_locs$seg <- paste0('seg', 1:dim(scdna_matrix_locs)[1])
  scdna_matrix_locs[['Gene_ID']] = as.character(scdna_matrix_locs[["seg"]])
  scdna_matrix_locs[['Gene_Chr']] = as.character(scdna_matrix_locs[["chr"]])
  row.names(scdna_matrix_locs) = scdna_matrix_locs[['Gene_ID']]
  rownames(scdna_matrix_merge) <- scdna_matrix_locs$seg

  gene_chr = rle(scdna_matrix_locs[["Gene_Chr"]][match(rownames(scdna_matrix_merge), scdna_matrix_locs[["Gene_ID"]])])
  gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
  names(gene_chr_mark) = gene_chr$values

  # Optionally remove specified celltypes
  if (!is.null(celltype)) {
    keep_cells <- which(seurat_obj$celltype %in% celltype)
    seurat_obj <- subset(seurat_obj, cells = keep_cells)
  }

  # Prepare matrix and labels
  mat <- scdna_matrix_merge[, colnames(seurat_obj)]
  label <- as.factor(seurat_obj$subclones)
  levels(label) <- 1:(length(unique(label)) + 1)
  label <- as.numeric(label)

  # Recorder cells in each subclones based on the distance to normal CNV profile
  normal_CNV <- matrix(2, nrow = 1, ncol = dim(scdna_matrix_merge)[1])
  d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
  d <- d[1:dim(mat)[2]+1]

  # Sort cells and cap CNV values
  mat <- SortCells(mat, label, d)
  mat[mat > max_cnv_value] = max_cnv_value
  labels <- rep(names(table(seurat_obj$subclones)), table(seurat_obj$subclones))
  labels <- data.frame(cluster=labels)

  # Optionally reorder subclones using hierarchical clustering
  if (!is.null(hc_rds) && file.exists(hc_rds)) {
    hc <- readRDS(hc_rds)
    obj <- Sort_subclones(mat, labels, hc)
    mat <- obj[[1]]
    labels <- obj[[2]]
  }

  # Plot and save heatmap
  filename <- paste0(celltype, '_', output_pdf)
  Plot_CNV_heatmap(mat, labels, gene_chr_mark, filename)
  invisible(mat)
}

#' Sort subclones based on hierarchical clustering order
#' This function reorders subclones in a matrix and their corresponding labels
#' based on the order from hierarchical clustering results.
#' @param mat A matrix of CNV values (cells x features) to be reordered
#' @param labels A data.frame with a 'cluster' column containing subclone labels
#' @param hc A hierarchical clustering object (e.g., from hclust) with 'order' and 'labels' components
#' @return A list containing:
#'   \item{mat}{The reordered matrix with cells sorted by subclone order}
#'   \item{labels}{The reordered labels data.frame with cluster factor levels set to match the new order}
#' @export
Sort_subclones <- function(mat, labels, hc){
  # order <- paste0('C', hc$order)
  order <- hc$labels[hc$order]
  index <- c()
  for (i in order){
    index_sub <- which(labels$cluster == i)
    index <- c(index, index_sub)
  }
  tmp_cluster <- labels$cluster[index]
  labels$cluster <- factor(tmp_cluster, levels = unique(tmp_cluster))
  mat <- mat[index, ]
  return(list(mat, labels))
}

#' Sort cells within a subclone by correlation
#' This function sorts cells within a subclone based on a correlation vector.
#' @importFrom dplyr %>% arrange
#' @param scdna_matrix A matrix of CNV values (features x cells)
#' @param distance A numeric vector of correlation or distance values for sorting
#' @return A matrix with columns (cells) sorted by distance
#' @export
SortCells_in_subclone <- function(scdna_matrix, distance) {
  scdna_matrix <- as.data.frame(t(scdna_matrix))
  scdna_matrix$distance <- distance
  scdna_matrix <- scdna_matrix %>% arrange(distance)
  scdna_matrix$distance <- NULL
  scdna_matrix <- t(scdna_matrix)
  return(scdna_matrix)
}

#' Sort cells by subclone and correlation
#'
#' This function sorts cells first by subclone label, then within each subclone by correlation.
#'
#' @param scdna_matrix A matrix of CNV values (features x cells)
#' @param label A vector of subclone labels (numeric or factor)
#' @param distance A numeric vector of correlation or distance values for sorting
#' @return A matrix with columns (cells) sorted by subclone and distance
#' @export
SortCells <- function(scdna_matrix, label, distance) {
  clusterlabel <- c()
  barcodes <- c()
  tmp_matrix <- as.data.frame(matrix(0, dim(scdna_matrix)[1], dim(scdna_matrix)[2]))
  for (j in 1:length(unique(label))) {
    clusterlabel <- c(clusterlabel, length(which(label == j)))
    barcodes <- c(barcodes, colnames(scdna_matrix)[which(label == j)])
    if(j == 1){
      tmp_matrix[, 1:clusterlabel] = scdna_matrix[, which(label == j)]
      tmp_matrix[, 1:clusterlabel] <- SortCells_in_subclone(tmp_matrix[, 1:clusterlabel], distance[which(label == j)])
    }else{
      a = sum(clusterlabel[1:j-1]) + 1
      b = sum(clusterlabel[1:j])
      tmp_matrix[, a:b]=scdna_matrix[, which(label == j)]
      tmp_matrix[, a:b] <- SortCells_in_subclone(tmp_matrix[, a:b], distance[which(label == j)])
    }
  }
  tmp_matrix <- as.matrix(tmp_matrix)
  rownames(tmp_matrix) <- rownames(scdna_matrix)
  colnames(tmp_matrix) <- barcodes
  tmp_matrix <- t(tmp_matrix)
  return(tmp_matrix)
}

#' Plot subclonal CNV heatmap
#'
#' This function plots a CNV heatmap for subclones, with chromosome and subclone annotations.
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation decorate_heatmap_body
#' @importFrom grid gpar grid.lines unit
#' @importFrom grDevices pdf dev.off
#' @param mat A matrix of CNV values (cells x features)
#' @param labels A data.frame with a 'cluster' column for subclone labels
#' @param gene_chr_mark A named vector marking chromosome boundaries
#' @param filename Output PDF file name
#' @return Invisibly returns NULL
#' @export
Plot_CNV_heatmap <- function(mat, labels, gene_chr_mark, filename) {
  ht_font = 50
  grid_height = 1
  col_fun = colorRamp2(c(0,1,2,3,4,5,6), c('blue', 'light blue', 'white', 'light salmon', 'coral', 'red', 'dark red'))
  sort_cluster = sort(as.character(unique(labels$cluster)))
  if(length(unique(labels$cluster)) > 9){
    mycol = brewer.pal(9, 'Set1')
    Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666", brewer.pal(6, 'Set1'), 'black')
    names(Label_color) <- sort_cluster
  }else if(length(unique(labels$cluster)) >= 3 & length(unique(labels$cluster)) <= 9){
    Label_color <- brewer.pal(length(unique(labels$cluster)), "Set1")
    names(Label_color) <- sort_cluster
  }else if(length(unique(labels$cluster)) == 2){
    Label_color <- c("#e41a1c", "#377eb8")
    names(Label_color) <- sort_cluster
  }else if(length(unique(labels$cluster)) == 1){
    Label_color <- c('#e41a1c')
    names(Label_color) <- sort_cluster
  }

  column_ha = HeatmapAnnotation(Chr = ComplexHeatmap::anno_mark(at=gene_chr_mark[1:22],
                                                                side="bottom",
                                                                labels=names(gene_chr_mark)[1:22],
                                                                labels_gp = grid::gpar(fontsize = ht_font)))

  ht1 = Heatmap(mat,name = 'CNV', na_col = 'white', show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                col = col_fun,
                border = TRUE,
                row_split = labels$cluster,
                row_title = NULL,
                bottom_annotation = column_ha,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_gap = grid::unit(0.3, "in"),
                column_gap = grid::unit(1, 'mm'),
                height = grid::unit(18, 'in'),
                width = grid::unit(36, 'in'),
                heatmap_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                            direction = "horizontal",
                                            grid_width = grid::unit(grid_height, 'inch'),
                                            grid_height = grid::unit(grid_height,"inch" ),
                                            labels_gp = gpar(fontsize = 0.8 * ht_font)))

  ht2 = Heatmap(labels$cluster,name = 'Subclone', col = Label_color, show_row_dend = FALSE,
                show_column_dend = FALSE,
                show_row_names = FALSE,
                show_column_names = FALSE,
                row_split = labels$cluster,
                row_title = NULL,
                border = TRUE,
                cluster_columns = FALSE,
                cluster_rows = FALSE,
                row_gap = grid::unit(0.3, "in"),
                width = grid::unit(1, 'in'),
                heatmap_legend_param = list(title_gp = gpar(fontsize = ht_font, fontface = "bold"),
                                            direction = "horizontal",
                                            grid_width = grid::unit(grid_height, 'inch'),
                                            grid_height = grid::unit(grid_height,"inch" ),
                                            labels_gp = gpar(fontsize = 1 * ht_font)))
  ht_list = ht1 + ht2

  pdf(file=filename, width=44, height=20)
  ComplexHeatmap::draw(ht_list,
                       ht_gap = unit(1, 'cm'),
                       heatmap_legend_side = "left",
                       annotation_legend_side = 'left',
                       merge_legend = T
  )
  chr_line <- gene_chr_mark/dim(mat)[2]
  chr_line[1] <- 0
  chr_line <- c(chr_line, 1)
  for (i in 1:length(unique(labels$cluster))){
    lapply(chr_line,
           FUN = function(p){
             decorate_heatmap_body("CNV", {
               grid.lines(c(p, p), c(0, 1), gp = gpar(lty = 1, lwd = 5))
             }, slice = i)
           }
    )
  }
  dev.off()
  invisible(NULL)
}

