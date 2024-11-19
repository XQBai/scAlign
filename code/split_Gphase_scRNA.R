library(Seurat)

seu_obj <- readRDS('seurat_epi.rds')

seu_obj <- subset(seu_obj, cells = which(seu_obj$Phase != 'S'))

saveRDS(seu_obj, file = 'seurat_epi.rds')