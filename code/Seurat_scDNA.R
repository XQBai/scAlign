library(dplyr)
library(mixtools)
library(GenomicRanges)
library(RColorBrewer)

source('/mnt/ix1/Projects/M070_200622_GI_multiomics/GithubCode_202411/code/Plot_all_chr_heatmap_latest.R', echo=TRUE)
source('/mnt/ix1/Projects/M070_200622_GI_multiomics/GithubCode_202411/code/Filter_noise_scDNA.R', echo=TRUE)

### Run seurat cluster
run_seurat <- function(mat){
  
  row.names(mat) <- paste0('segement', 1:dim(mat)[1])
  dims.reduce <- 50
  ncells <- dim(mat)[2]
  
  r <- 0.5
  # # Set the resolution parameter in Seurat 
  # if ( ncells <= 100 ){
  #   r <- 0.1
  # }else if(100 < ncells & ncells < 1200){
  #   r <- 0.2
  # }else if (1200 <= ncells & ncells  < 1800){
  #   r <- 0.3
  # }else if (1800 <= ncells & ncells < 2400){
  #   r <- 0.4
  # }else{
  #   r <- 0.5
  # }
  
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
                                     project='seurat-v4', min.cells = -Inf, 
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