## Identify cellular components 

library(ggplot2)
library(optparse)
library(future)
library(Seurat)
library(IRanges)
library(data.table)
library(dplyr)
library(mixtools)
library(GenomicRanges)
library(RColorBrewer)


# options(future.globals.maxSize = 24000 * 1024^2)
# plan('multiprocess', workers=24)
# 
# # project_path <- '../example/P5931/scDNA'
# # scdna_matrix_locs <- readRDS('scdna_matrix_locs.rds')
# # seurat_scDNA <- readRDS('seurat_scDNA.rds')
# 
# ## read the cell metrics file conducted by cellranger-dna
# per_cell_metrics <- read.csv(file.path(project_path, 'per_cell_summary_metrics.csv'), sep=",", header=T)

# Identify the cellranger noise cells by cellranger metrics 
Identify_cellranger_noise <- function(per_cell_metrics, seurat_scDNA){
  
  celltype <- matrix(0, nrow = dim(seurat_scDNA)[2])
  celltype[which(per_cell_metrics$is_noisy == 1)] <- 'cellranger noise'
  seurat_scDNA$celltype <- celltype
  return(seurat_scDNA)
}

# Identify the technical noise cells 
Identify_technoise <- function(seurat_scDNA){
  
  scdna_matrix_merge_allbarcodes <- seurat_scDNA@assays$RNA@counts
  num_NA <- apply(scdna_matrix_merge_allbarcodes, 2, function(x){length(which(x==0))})
  prop_NA <- num_NA/dim(scdna_matrix_merge_allbarcodes)[1]
  prop_format <- data.frame(prop = prop_NA, index = seq(1, length(prop_NA)))
  
  p<- ggplot(data=prop_format, aes(x= index, y=prop)) + 
    geom_point(size = 1, alpha = 0.5) + 
    xlab("Cell index") +
    ylab("CNV dropout proportion") +
    scale_y_continuous(breaks= seq(0, max(prop_format$prop), 0.1)) + 
    geom_hline(aes(yintercept = 0.1, colour = '#990000', linetype = 'dashed'), show.legend = FALSE) + 
    theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
    theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
    theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
    theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
    theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) 
  ggsave(p, file = 'Miss_val_proportion.pdf')
  
  technoise_index <- which(prop_NA > 0.1)
  
  celltype <- seurat_scDNA$celltype
  celltype[technoise_index] <- 'technical noise'
  seurat_scDNA$celltype <- celltype
  
  return(seurat_scDNA)
}

### Distinguish the noise cells and S phase cells from the union of cellranger noise and S phase
ReIdentify_replicates <- function(seurat_scDNA){
  
  celltype <- seurat_scDNA$celltype
  index <- union(which(celltype == 'cellranger noise'), which(celltype == 'S phase'))
  seurat.noise <- subset(seurat_scDNA, cells = index)
  noise.mat <- seurat.noise@assays$RNA@counts
  noise.mean <- apply(noise.mat, 2, mean)
  
  identifylabel <- rep(0, length(index))
  identifylabel[which(noise.mean >= 2.5)] <- 'S phase'
  identifylabel[which(noise.mean < 2.5)] <- 'cellranger noise'
  celltype[index] <- identifylabel
  seurat_scDNA$celltype <- celltype
  return(seurat_scDNA)
}

### Identify the normal cells
Identify_normal <- function(seurat_scDNA, scdna_matrix_locs){
  
  # Filter the noise cells
  celltype <- seurat_scDNA$celltype
  seurat_scDNA_sub <- subset(seurat_scDNA, cells = which(celltype == 0))

  tmp_matrix <- as.data.frame(seurat_scDNA_sub@assays$RNA@counts)
  tmp_matrix[tmp_matrix == 0] <- NA
  tmp_matrix$chr <- scdna_matrix_locs$chr
  
  cell_averages <- tmp_matrix %>% group_by(chr) %>% summarise_all(list(mean), na.rm=T) %>% as.data.frame()
  cell_averages$chr <- NULL
  cell_averages_melt <- data.table::melt(cell_averages, measure.vars=colnames(cell_averages))
  # detach(package:plyr)
  cell_averages_classification <- cell_averages_melt %>% dplyr::group_by(variable) %>% 
    dplyr::summarise(c=any(value > 2.5, na.rm=T))
  
  celltype[which(celltype == 0)] <- cell_averages_classification$c
  celltype[which(celltype == "FALSE")] <- 'normal'
  celltype[which(celltype == 'TRUE')] <- 'tumor'
  seurat_scDNA$celltype <- celltype
  return(seurat_scDNA)
}

## Identify the replication cells
Identify_replicates <- function(seurat_scDNA, scdna_matrix_locs){
  
  celltype <- seurat_scDNA$celltype
  # subset the object with tumor cells
  seurat_scDNA_tumor <- subset(seurat_scDNA, cells = which(celltype == 'tumor'))
  # Filter the genome segments on chromosome X and Y
  locs <- scdna_matrix_locs %>% dplyr::filter(chr != 'chrX' & chr != 'chrY')
  
  mat <- seurat_scDNA_tumor@assays$RNA@counts[1:length(locs$chr), ]
  
  # set the normal CNV profile as CN=2
  normal_CNV <- matrix(2, nrow = 1, ncol = dim(mat)[1])
  d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
  d <- d[1:dim(mat)[2]+1]

  # Fit the mixture models for the distances from tumor cell to normal by EM algorithm
  estimate_mix_model <- normalmixEM(d)
  est <- data.frame(lambda = estimate_mix_model$lambda, mu = estimate_mix_model$mu, sigma = estimate_mix_model$sigma)
  
  pdf('mixture_cells.pdf')
  mixtools::plot.mixEM(estimate_mix_model, whichplots = 2, main2 = 'S phase cells', xlab2 = 'The Euclidean distance to normal cells' , lwd2 = 3, marginal = TRUE)
  dev.off()
  
  #S phase cells from two categories, either from the distribution with big mu and big sigma, or
  #the distribution from small mu with big sigma
  mu <- estimate_mix_model$mu
  lambda <- estimate_mix_model$lambda
  sigma <- estimate_mix_model$sigma
  index_label <- apply(estimate_mix_model$posterior, 1, function(x){which(x == max(x))})
  if(max(mu) > 2*min(mu)){
    index_Sphase <- which(mu == max(mu))
  }else{
    index_Sphase <- which(sigma == max(sigma))
  }
  
  index_Scells <- which(index_label == index_Sphase)
  
  index_label[index_Scells] <- 'S phase'
  index_label[which(index_label != 'S phase')] <- 'G0G1'
  celltype[which(celltype == 'tumor')] <- index_label
  seurat_scDNA$celltype <- celltype
  return(seurat_scDNA)
}








