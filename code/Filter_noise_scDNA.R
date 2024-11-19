library(IRanges)
library(GenomicRanges)
library(data.table)
library(dplyr)
library(Seurat)
library(optparse)
library(ggplot2)
library(dplyr)
library(future)
library(pheatmap)
# library(mixdist)
library(mixtools)
library(GenomicRanges)


Iden_signal_segements <- function(scdna_matrix){
  
  arm_sd <- apply(scdna_matrix, 1, sd)
  arm_mean <- apply(scdna_matrix-2, 1, mean)
  arm_sd <- arm_sd[union(which(arm_sd != 0), which(arm_mean != 0))]
  # calculate the distance in each rows
  # arm_sd <- apply(scdna_matrix-2, 1, function(x){sqrt(sum(x^2))})
  
  # pdf('./Segements_sd_density.pdf')
  # plot(density(arm_sd), main='Density of sd', ylab = 'Density')
  # dev.off()
  
  pdf('./Segements_density_curve.pdf')
  mixtools::plot.mixEM(normalmixEM(arm_sd), whichplots = 2, main2 = 'Selection of informative segments', xlab2 = 'Standard deviation of CNVs on segments' ,  
                       lwd2 = 3, marginal = TRUE)
  dev.off()
  
  estimate_mix_model <- normalmixEM(arm_sd)
  index_seg_cluster <- which(estimate_mix_model$lambda == max(estimate_mix_model$lambda))
  index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
  index_segement <- which(index_label != index_seg_cluster)
  
  index_segement_tmp <- index_segement
  index_seg_cluster_tmp <- index_seg_cluster
  
  ## Stop select when length of selected segements is less than 20% of total
  while(length(index_segement) <= round(0.2 * dim(scdna_matrix)[1])){
    
    if (length(index_segement) >= round(0.21 * dim(scdna_matrix)[1])){
      break 
    }
    segement_index_remain <- which(index_label == index_seg_cluster_tmp)
    arm_sd_tmp <- arm_sd[segement_index_remain]
    estimate_mix_model_tmp <- normalmixEM(arm_sd_tmp)
    index_seg_cluster_tmp <- which(estimate_mix_model_tmp$lambda == max(estimate_mix_model_tmp$lambda))
    index_label_tmp <- apply(estimate_mix_model_tmp$posterior, 1, function(x){which(x == max(x))})
    index_segement_tmp <- segement_index_remain[which(index_label_tmp != index_seg_cluster_tmp)]
    index_segement <- c(index_segement, index_segement_tmp)
  }
  return(index_segement)
}

