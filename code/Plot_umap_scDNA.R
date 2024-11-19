library(ape)
library(Seurat)
library(dplyr)
library(RColorBrewer)
library(ggplot2)

seurat_scDNA_subset <- readRDS('seurat_scDNA_tumor.rds')
labels <- seurat_scDNA_subset$subclones

if(length(unique(labels)) > 9){
  mycol = brewer.pal(9, 'Set1')
  Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666", brewer.pal(5, 'Set1'))
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) >= 3 & length(unique(labels)) <= 9){
  Label_color <- brewer.pal(length(unique(labels)), "Set1")
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) == 2){
  Label_color <- c("#e41a1c", "#377eb8")
  names(Label_color) <- sort(unique(labels))
}else if(length(unique(labels)) == 1){
  Label_color <- c('#e41a1c')
  names(Label_color) <- sort(unique(labels))
}

pdf('scDNA_umap_subclones.pdf')
DimPlot(seurat_scDNA_subset, group.by = 'subclones', cols = Label_color, reduction = "umap") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  labs(title = '') + 
  guides(color=guide_legend("Subclone")) +
  theme(plot.title = element_text(size = 11, color = 'black', face = 'bold', hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.text.y = element_text(size = 15, color = 'black', face = 'bold')) +
  theme(axis.title.x = element_text(size = 15, color = 'black', face = 'bold')) + 
  theme(axis.title.y = element_text(size = 15, color = 'black', face = 'bold')) + 
  theme(legend.text= element_text(size=15,color="black", face= "bold", vjust=0.5, hjust=0.5)) +
  theme(legend.title = element_text(size=15,color="black", face= "bold"))
dev.off()



