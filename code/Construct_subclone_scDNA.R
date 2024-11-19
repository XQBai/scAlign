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

source('Plot_all_chr_heatmap_latest.R', echo=TRUE)
source('Filter_noise_scDNA.R', echo=TRUE)

options(future.globals.maxSize = 24000 * 1024^2)
plan('multiprocess', workers=24)


### Run seurat cluster
run_seurat <- function(mat){
  
  row.names(mat) <- paste0('segement', 1:dim(mat)[1])
  dims.reduce <- 50
  ncells <- dim(mat)[2]
  
  # Set the resolution parameter in Seurat 
  if ( ncells <= 100 ){
    r <- 0.1
  }else if(100 < ncells & ncells < 1200){
    r <- 0.2
  }else if (1200 <= ncells & ncells  < 1800){
    r <- 0.3
  }else if (1800 <= ncells & ncells < 2400){
    r <- 0.4
  }else{
    r <- 0.5
  }
  
  ## Feature selection 
  y <- tryCatch({
    Iden_signal_segements(mat)
  },warning= function(w){
    # index_segement <- Iden_signal_segements(mat)
    # print('warning')
    return(0)
  }, error = function(e){
    # index_segement <- seq(1, dim(mat)[1])
    # print('error')
    return(0)
  })
  
  if(y == 0){
    index_segement <- seq(1, dim(mat)[1])
  }else{
    index_segement <- Iden_signal_segements(mat)
  }
  
  seurat_scDNA <- CreateSeuratObject(counts = as.matrix(mat), 
                                     project='seurat-v3', min.cells = -Inf, 
                                     min.features = -Inf)
  
  signal_segement <- rownames(mat)[index_segement]
  seurat_scDNA<- subset(seurat_scDNA, features = signal_segement)
  seurat_scDNA <- seurat_scDNA %>% ScaleData(do.center=F, do.scale = F) %>% 
    RunPCA(verbose=T, features=rownames(seurat_scDNA), npcs = min(50, ncol(seurat_scDNA)-1)) %>% 
    RunUMAP(umap.method = 'uwot', n.neighbors = min(20, ncol(seurat_scDNA)), dims = 1:min(20, ncol(seurat_scDNA)-1), verbose=F) %>% 
    FindNeighbors(reduction='umap', dims=1:2)%>% 
    FindClusters(reduction = 'umap', resolution = r)
  
  Idents(seurat_scDNA) <- paste0('C', as.numeric(Idents(seurat_scDNA)))
  return(seurat_scDNA)
}

## Run Seurat to detect main clusters in G0G1 tumor cells
Construct_subclones <- function(seurat_scDNA){
  
  celltype <- seurat_scDNA$celltype
  seurat_scDNA_subset <- subset(seurat_scDNA, cells = which(celltype == 'G0G1'))
  seurat_scDNA_subset <- run_seurat(seurat_scDNA_subset@assays$RNA@counts)
  seurat_scDNA_subset$subclones <- paste0('C', as.numeric(seurat_scDNA_subset$seurat_clusters))
  
  celltype <- seurat_scDNA$celltype
  celltype[which(celltype == 'G0G1')] <- paste0('C', as.numeric(seurat_scDNA_subset$seurat_clusters))
  seurat_scDNA$subclones <- celltype
  saveRDS(seurat_scDNA_subset, file = 'seurat_scDNA_tumor.rds')
  return(seurat_scDNA)
}


