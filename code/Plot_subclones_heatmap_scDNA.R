library(future)
library(optparse)
library(Seurat)
library(IRanges)
library(data.table)
library(dplyr)
library(ggplot2)
# library(mixdist)
library(mixtools)
library(GenomicRanges)


options(future.globals.maxSize = 6000*1024^2)
plan('multiprocess', workers = 12)

source('Plot_all_chr_heatmap_latest.R', echo=TRUE)

# load('scdna_matrix_locs.Robj')
# load('seurat_scDNA.Robj')
# hc <- readRDS('hc.rds')
# hc_no_normal <- readRDS('hc_no_normal.rds')

scdna_matrix_merge <- as.data.frame(seurat_scDNA@assays$RNA@counts)
scdna_matrix_merge$chr <- scdna_matrix_locs$chr

scdna_matrix_merge <- scdna_matrix_merge %>% filter(chr != 'chrX' & chr != 'chrY') 
scdna_matrix_merge$chr <- NULL
scdna_matrix_locs <- scdna_matrix_locs %>% filter(chr != 'chrX' & chr != 'chrY') 

sum_row <- apply(scdna_matrix_merge, 1, mean)
row_index <- which(sum_row != 0)
scdna_matrix_merge <- scdna_matrix_merge[row_index, ]
scdna_matrix_locs <- scdna_matrix_locs[row_index, ]

scdna_matrix_merge$chr <- scdna_matrix_locs$chr
scdna_matrix_merge$bin <- scdna_matrix_locs$bin

scdna_matrix_merge <- scdna_matrix_merge %>% arrange(bin)
scdna_matrix_locs <- scdna_matrix_locs %>% arrange(bin)

scdna_matrix_merge$chr <- NULL
scdna_matrix_merge$bin <- NULL

scdna_matrix_locs$seg <- paste0('seg', 1:dim(scdna_matrix_locs)[1])
scdna_matrix_locs[['Gene_ID']] = as.character(scdna_matrix_locs[["seg"]])
scdna_matrix_locs[['Gene_Chr']] = as.character(scdna_matrix_locs[["chr"]])
row.names(scdna_matrix_locs) = scdna_matrix_locs[['Gene_ID']]

rownames(scdna_matrix_merge) <- scdna_matrix_locs$seg

#mat = t(CNVsample) #reduce size for debugging
gene_chr = rle(scdna_matrix_locs[["Gene_Chr"]][match(rownames(scdna_matrix_merge), scdna_matrix_locs[["Gene_ID"]])])
gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
names(gene_chr_mark) = gene_chr$values
# View(gene_chr_mark)


#####################################################################
## Plot figure with noise cells
mat = scdna_matrix_merge
label <- as.factor(seurat_scDNA$subclones)
levels(label) <- 1:(length(unique(label)) + 1)
label <- as.numeric(label)

# Recorder cells in each subclones based on the distance to normal CNV profile
normal_CNV <- matrix(2, nrow = 1, ncol = dim(scdna_matrix_merge)[1])
d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), mat))))[1, ]
d <- d[1:dim(mat)[2]+1]

mat <- SortCells(mat, label, d)
mat[mat > 6] = 6
labels <- rep(names(table(seurat_scDNA$subclones)), table(seurat_scDNA$subclones))
labels <- data.frame(cluster=labels)
# Do not need to reorder the subclones
Plot_CNV_heatmap(mat, labels, gene_chr_mark, 'Subclones_with_noise.pdf')

## Plot figure when removing normal cells
seurat.no.normal <- subset(seurat_scDNA, cells = which(seurat_scDNA$celltype == 'G0G1'))
labels.no.normal <- seurat.no.normal$subclones
labels.no.normal <- data.frame(cluster=labels.no.normal)
mat.no.normal <- scdna_matrix_merge[, rownames(labels.no.normal)]
mat.no.normal <- SortCells(mat.no.normal, label[which(seurat_scDNA$celltype == 'G0G1')], d[which(seurat_scDNA$celltype == 'G0G1')])
mat.no.normal[mat.no.normal > 6] = 6

label.no.normal <- rep(names(table(seurat.no.normal$subclones)), table(seurat.no.normal$subclones))
label.no.normal <- data.frame(cluster=label.no.normal)

if(length(unique(label.no.normal$cluster)) > 1){
  obj <- Sort_subclones(mat.no.normal, label.no.normal, hc_no_normal)
  mat.no.normal <- obj[[1]]
  label.no.normal <- obj[[2]]
}

Plot_CNV_heatmap(mat.no.normal, label.no.normal, gene_chr_mark, 'Subclones_no_normal.pdf')

## Plot figure with normal cells
index <- union(which(seurat_scDNA$celltype == 'normal'), which(seurat_scDNA$celltype == 'tumor'))
seurat.normal <- subset(seurat_scDNA, cells = index)
labels.normal <- seurat.normal$subclones
labels.normal <- data.frame(cluster=labels.normal)

labelnormal <- as.factor(seurat.normal$subclones)
levels(labelnormal) <- 1:(length(unique(labelnormal)) + 1)
labelnormal <- as.numeric(labelnormal)

mat.normal <- scdna_matrix_merge[, colnames(seurat.normal)]
mat.normal <- SortCells(mat.normal, labelnormal , d[index])
mat.normal[mat.normal > 6] = 6

label.normal <- rep(names(table(seurat.normal$subclones)), table(seurat.normal$subclones))
label.normal <- data.frame(cluster=label.normal)
obj <- Sort_subclones(mat.normal, label.normal, hc)
mat.normal <- obj[[1]]
labels.normal <- obj[[2]]
Plot_CNV_heatmap(mat.normal, labels.normal, gene_chr_mark, 'Subclones_no_noise.pdf')


## Plot figure for technical noise cells
seurat.tech <- subset(seurat_scDNA, cells = which(seurat_scDNA$celltype == 'technical noise'))
labels.tech <- seurat.tech$subclones
labels.tech <- data.frame(cluster=labels.tech)
mat.tech <- scdna_matrix_merge[, colnames(seurat.tech)]
mat.tech <- SortCells(mat.tech, matrix(1, nrow = dim(seurat.tech)[2]) , d[which(seurat_scDNA$celltype == 'technical noise')])
mat.tech[mat.tech > 6] = 6
Plot_CNV_heatmap(mat.tech, labels.tech, gene_chr_mark, 'Subclones_technical_noise.pdf')


## Plot figure for cellranger noise and s phase cells
index <- union(which(seurat_scDNA$celltype == 'cellranger noise'), which(seurat_scDNA$celltype == 'S phase'))
seurat.noise <- subset(seurat_scDNA, cells = index)
labels.noise <- seurat.noise$subclones
labels.noise <- data.frame(cluster=labels.noise)

labelnoise <- as.factor(seurat.noise$subclones)
levels(labelnoise) <- 1:(length(unique(labelnoise)) + 1)
labelnoise <- as.numeric(labelnoise)

mat.noise <- scdna_matrix_merge[, colnames(seurat.noise)]
mat.noise <- SortCells(mat.noise, labelnoise , d[index])
mat.noise[mat.noise > 6] = 6

Plot_CNV_heatmap(mat.noise, labels.noise, gene_chr_mark, 'Subclones_S_phase.pdf')



