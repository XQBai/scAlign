## Author: Sue Grimes <sgrimes@stanford.edu> 
## Merged normal samples' objects with tumor samples' objects 
library(Seurat)
library(Matrix)
library(dplyr)

modifyBarcodeSuffix <- function(barcodes, cell.suffix) {
  bc_has_suffix = (length(grep("-", barcodes)) > 0)
  
  if (cell.suffix == 'N' && bc_has_suffix) {
    new_barcodes = gsub("-\\d", '', barcodes)  #Delete barcode suffix
    
  } else if (cell.suffix %in% c('K','N')) {
    new_barcodes = barcodes
    
  } else if (bc_has_suffix) {
    new_bc = paste0('-',cell.suffix)
    new_barcodes = gsub("-\\d", new_bc, barcodes)  #Replace barcode suffix
    
  } else {
    new_barcodes = paste(barcodes, paste0('-',cell.suffix), sep="")
  }
  return(new_barcodes)
}  

robjs = list.files('.', pattern="*seurat_dfx.rds", full.names=F)
nr_seu = length(robjs)
sample_fn = gsub('.seurat_dfx.rds','',robjs)

seu <- vector("list", nr_seu)
for (i in 1:nr_seu) {
  sprintf("Loading Seurat object: %s", robjs[i])
  seu[[i]] <- readRDS(robjs[i])
  
  barcodes <- colnames(seu[[i]])
  seu[[i]] <- RenameCells(seu[[i]], new.names=modifyBarcodeSuffix(barcodes, i))
  
  fn_parts <- unlist(strsplit(sample_fn[i], '_', fixed=T))
  seu[[i]]$patientID = fn_parts[1]
  seu[[i]]$sample    = fn_parts[2]
  seu[[i]]$condition = fn_parts[3]
  
  sprintf("%i cells loaded", length(Idents(seu[[i]])))
}

seu_merged <- merge(x=seu[[1]], y=seu[-1])

sprintf("Aggregate Seurat object has %i cells", length(Idents(seu_merged)))
save(seu_merged, file="seurat_merged.Robj")

#checking
head(seu_merged@meta.data)
tail(seu_merged@meta.data)
table(seu_merged@meta.data$orig.ident, seu_merged@meta.data$condition)

