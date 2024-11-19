library(dplyr)
library(GenomicRanges)
library(Seurat)
library(ggplot2)
library(future)

options(future.globals.maxSize = 6000 * 1024^2)
plan('multisession', workers= 18)

source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/code_v1/Plot_all_chr_heatmap_latest.R', echo=TRUE)

expr.infercnv.17 <- readRDS('./test/preliminary.infercnv_obj')
pred_cnv_genes <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.genes_used.dat', header = T)
pred_cnv_regions <- read.table('./test/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat', header = T)

identical(rownames(expr.infercnv.17), rownames(pred_cnv_genes))

### Only keep 1:22 chromosomes
pred_cnv_genes <- pred_cnv_genes %>% filter(chr %in% seq(1:22))
expr.infercnv.17 <- expr.infercnv.17[rownames(pred_cnv_genes), ]


idents <- read.table('idents.txt')
cellbarcode <- gsub('-', '.', colnames(expr.infercnv.17))
identical(cellbarcode, colnames(expr.infercnv.17))
##############################

###############
scdna_matrix_merge <- expr.infercnv.17
scdna_matrix_locs <- pred_cnv_genes
scdna_matrix_locs$chr <- paste0('chr', scdna_matrix_locs$chr)
scdna_matrix_merge$chr <- scdna_matrix_locs$chr
chrs <- paste0('chr', 1:22)

### Remove the normal cells
if(length(which(idents$V2 == 'normal')) == 0){
  scdna_matrix_merge = scdna_matrix_merge
  idents = idents
}else{
  scdna_matrix_merge <- scdna_matrix_merge[, which(idents$V2 != 'normal')]
  idents <- idents %>% filter(V2 != 'normal')
}

scdna_matrix_merge$chr <- NULL

sum_row <- apply(scdna_matrix_merge, 1, mean)
row_index <- which(sum_row != 0)
scdna_matrix_merge <- scdna_matrix_merge[row_index, ]
scdna_matrix_locs <- scdna_matrix_locs[row_index, ]

scdna_matrix_locs$gene <- rownames(scdna_matrix_locs)
# scdna_matrix_locs$seg <- paste0('seg', 1:dim(scdna_matrix_locs)[1])
scdna_matrix_locs[['Gene_ID']] = as.character(scdna_matrix_locs[["gene"]])
scdna_matrix_locs[['Gene_Chr']] = as.character(scdna_matrix_locs[["chr"]])
row.names(scdna_matrix_locs) = scdna_matrix_locs[['Gene_ID']]
rownames(scdna_matrix_merge) <- scdna_matrix_locs$gene

#mat = t(CNVsample) #reduce size for debugging
gene_chr = rle(scdna_matrix_locs[["Gene_Chr"]][match(rownames(scdna_matrix_merge), scdna_matrix_locs[["Gene_ID"]])])
gene_chr_mark = c(0, cumsum(gene_chr$lengths))[seq_along(gene_chr$lengths)]+1
names(gene_chr_mark) = gene_chr$values
# View(gene_chr_mark)
scdna_matrix_merge[which(scdna_matrix_merge > 6)] = 6
# mat.normal = scdna_matrix_merge[, which(idents$V2 == 'normal')]
# mat.tumor = scdna_matrix_merge[, which(idents$V2 != 'normal')]

# Recorder cells in each subclones based on the distance to normal CNV profile
normal_CNV <- matrix(2, nrow = 1, ncol = dim(scdna_matrix_merge)[1])
d <- as.matrix(dist(t(cbind(t(as.matrix(normal_CNV)), scdna_matrix_merge))))[1, ]
d <- d[1:dim(scdna_matrix_merge)[2]+1]

label = as.factor(idents$V2)
levels(label) <- 1:(length(unique(label)) + 1)

mat = scdna_matrix_merge
mat <- SortCells(mat, label, d)
mat[mat > 6] = 6
mat <- mat-1
# mat <- mat - median(as.matrix(mat))
labels <- rep(names(table(idents$V2)), table(idents$V2))
labels <- data.frame(cluster=labels)

df <- mat
#### Plot RNA infercnv heatmap
# loading data containing vector chr
chr_lengths <- gene_chr$lengths
chr_binary <- rep(c(2, 1), 11)
chr <- data.frame(chr = rep.int(x = chr_binary, times = chr_lengths))

## Geting lengths for numbers annotation
chr_rl_c <- c(1, cumsum(chr_lengths))

# Creating lengths for chr numbers annotation
chr_df <- data.frame(a = chr_rl_c[1:length(chr_rl_c) - 1],
                     b = chr_rl_c[2:length(chr_rl_c)])

chr_1_means <- round(rowMeans(chr_df))

chrom.names <- c(1:22)

## creating the vectors for chr number annotations
v <- vector(length = sum(chr_lengths), mode = 'character')
v[chr_1_means] <- chrom.names
v[is.na(v)] <- ''

# chr bar with the chr names
chr_bar <- HeatmapAnnotation(
  chr_text = anno_text(v[1:ncol(df)],
                       location = 0.5,
                       #rot = 0,
                       just = 'center',
                       gp = gpar(fontsize = 16)),
  df = as.character(chr[1:dim(mat)[2], ]),
  show_legend = FALSE,
  show_annotation_name = FALSE,
  which = 'column',
  col = list(df=c('1' = 'grey88', '2' = 'black'))
)

## colors
ploidy_trunc = 6
ploidy_VAL = 2
color_heat <- structure(pals::ocean.balance(length(0 : ploidy_trunc)),
                        names = 0:ploidy_trunc) 

# special ploidy colors if ground state rounds to 2
if (round(ploidy_VAL) == 2) {
  color_heat <- structure(
    c(
      "#3787BA",
      "#95B8C5",
      "#F0ECEB",
      "#D7A290",
      "#BF583B",
      "#8D1128",
      "#3C0912"
    ),
    names = c("0", "1", "2", "3", "4", "5", "6")
  )
}


# set color for subclones
if(length(unique(labels$cluster)) > 9){
  mycol = brewer.pal(9, 'Set1')
  Label_color <- c(mycol, brewer.pal(12, 'Set3')[10:12], "#666666")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) >= 3 & length(unique(labels$cluster)) <= 9){
  Label_color <- brewer.pal(length(unique(labels$cluster)), "Set1")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) == 2){
  Label_color <- c("red", "green")
  names(Label_color) <- sort(unique(labels$cluster))
}else if(length(unique(labels$cluster)) == 1){
  Label_color <- c('red')
  names(Label_color) <- sort(unique(labels$cluster))
}

## Row annotation for subclones
ann <- rowAnnotation(
  df = labels$cluster,
  col = list(df = Label_color),
  border = FALSE,
  show_annotation_name = FALSE,
  gap = unit(1, "in"),
  annotation_name_gp = gpar(fontsize = 17),
  simple_anno_size = unit(0.5, 'cm'),
  annotation_legend_param = list(
    title = 'Subclones',
    subclones = list(
      labels = gtools::mixedsort(unique(as.character(labels$cluster))),
      at = gtools::mixedsort(unique(as.character(labels$cluster)))
    )
  ),
  show_legend = TRUE
)

popseg_heatmap <- df

pdf('infercnv_heatmap.pdf', width = 10, height = 5)
Heatmap(popseg_heatmap, name = 'infercnv',
        use_raster = TRUE,
        column_title = 'Genomic coordinates',
        column_title_side = 'bottom',
        top_annotation = chr_bar, 
        cluster_rows = FALSE,
        row_gap = unit(0.1, "in"),
        height = grid::unit(4, 'in'),
        width = grid::unit(8.5, 'in'),
        border = TRUE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = FALSE,
        col = color_heat,
        heatmap_legend_param = list(title = 'Inferred CNV'),
        show_heatmap_legend = TRUE,
        left_annotation = ann) 

chr_line <- gene_chr_mark/dim(mat)[2]
chr_line[1] <- 0
chr_line <- c(chr_line, 1)
lapply(chr_line,
       FUN = function(p){
         decorate_heatmap_body("infercnv", {
           grid.lines(c(p, p), c(0, 1), gp = gpar(lty = 1, lwd = 1))
         }, slice = 1)
       }
)
dev.off()

