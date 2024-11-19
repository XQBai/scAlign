library(readxl)
library(tibble)
library(dplyr)
library(plyr)
library(stringr)
library(gplots)
library(shiny)
# library(heatmaply)
library(RColorBrewer)
library(Hmisc) #for translate 
#library(rjags)
#library(devtools)
library(ComplexHeatmap)
library(circlize)

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

SortCells_in_subclone <- function(tmp_scdna_matrix, Corr){
  
  tmp_scdna_matrix <- as.data.frame(t(tmp_scdna_matrix))
  tmp_scdna_matrix$corr <- Corr
  tmp_scdna_matrix <- tmp_scdna_matrix %>% arrange(corr)
  tmp_scdna_matrix$corr <- NULL
  tmp_scdna_matrix <- t(tmp_scdna_matrix)
  return(tmp_scdna_matrix)
}


# b <- SortCells(CNVsample)
SortCells <- function(tmp_scdna_matrix_arm, S1, Corr){
  ## Re-order cells according to the cluster labels
  clusterlabel <- c()
  barcodes <- c()
  
  tmp_matrix <- as.data.frame(matrix(0, dim(tmp_scdna_matrix_arm)[1], dim(tmp_scdna_matrix_arm)[2]))
  for (j in 1:length(unique(S1))){
    clusterlabel <- c(clusterlabel, length(which(S1 == j)))
    barcodes <- c(barcodes, colnames(tmp_scdna_matrix_arm)[which(S1 == j)])
    if(j == 1){
      tmp_matrix[, 1:clusterlabel] = tmp_scdna_matrix_arm[, which(S1 == j)]
      tmp_matrix[, 1:clusterlabel] <- SortCells_in_subclone(tmp_matrix[, 1:clusterlabel], Corr[which(S1 == j)])
    }else{
      a = sum(clusterlabel[1:j-1]) + 1
      b = sum(clusterlabel[1:j])
      tmp_matrix[, a:b]=tmp_scdna_matrix_arm[, which(S1 == j)]
      tmp_matrix[, a:b] <- SortCells_in_subclone(tmp_matrix[, a:b], Corr[which(S1 == j)])
    }
  }
  tmp_matrix <- as.matrix(tmp_matrix)
  rownames(tmp_matrix) <- rownames(tmp_scdna_matrix_arm)
  colnames(tmp_matrix) <- barcodes
  tmp_matrix <- t(tmp_matrix)
  return(tmp_matrix)
}


Plot_CNV_heatmap <- function(mat, labels, gene_chr_mark, filename){
  
  ht_font = 50
  grid_height = 1
  col_fun = colorRamp2(c(0,1,2,3,4,5,6), c('blue', 'light blue', 'white', 'light salmon', 'coral', 'red', 'dark red'))
  # col_fun = colorRamp2(c(0, 1, 2, 3, 4, 5, 6), c(brewer.pal(3, 'Blues')[3], brewer.pal(3, 'Blues')[2], 'white', brewer.pal(4, 'Reds'))) 
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
                # row_km = 4, 
                #column_km = 22, 
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
                # row_km = 4, 
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

}

