library(dplyr)
library(proxy)
library(Seurat)
library(data.table)


source('/mnt/ix1/Projects/M070_200622_GI_multiomics/scDNA/00_Code/code_v1/Filter_noise.R', echo=TRUE)

########### Only common genes #############################
Preproc <- function(cnv.gene, rna.norm, cnv.clone){
  
  if (identical(rownames(cnv.gene), rownames(rna.norm)) && identical(rownames(cnv.gene), colnames(cnv.clone))){
    print("The genes are consistent between scRNA and scDNA")
  } else{
    d <- apply(rna.norm, 1, sum)
    rna.norm <- rna.norm[which(d != 0), ]
    common <- intersect(rownames(cnv.gene), rownames(rna.norm))
    rna.norm <- rna.norm[common, ]
    cnv.gene <- cnv.gene[common, ]
    cnv.clone <- as.matrix(cnv.clone[, common])
    }
  return(list(rna.norm, cnv.gene, cnv.clone))
}


############# Projection function #############################
### Y represents the DNA matrix
### X represents the RNA matrix
# X <- L[[1]]
# Y <- L[[2]]
# Y_proj <- Projection(X, Y)
Projection <- function(X, Y){

  beta <- matrix(0, ncol = 1, nrow = dim(X)[1])
  X <- t(X)
  Y <- t(Y)
  
  for (i in 1:dim(beta)[1]){
    beta[i] <- (sum(X[, i])*sum(Y[, i]))/dim(Y)[1]/(norm(X[, i], type = '2')^2)
  }
  
  ### Projection space
  Y_proj <- matrix(NA, nrow = dim(X)[1], ncol = dim(X)[2])
  
  for (i in 1:dim(Y_proj)[1]){
    Y_proj[i, ] <- X[i, ] * beta
  }
  return(Y_proj)
}


#### alignment by the euclidean distance between rna projection & dna
Cluster_Alignment <- function(cnv.clone, Y_proj){
  
  dis.mat <- apply(cnv.clone, 1, function(col1) {
    apply(Y_proj, 1, function(col2) {
      dist(rbind(col1,col2), method = "Euclidean")
      #cor(col1, col2, method = "pearson")
    })
  })

  label <- unlist(apply(dis.mat, 1, function(x){which(x == min(x))}))

  # dis.mat <- apply(cnv.clone, 1, function(col1) {
  #   apply(Y_proj, 1, function(col2) {
  #     cor(col1, col2, method = "pearson")
  #   })
  # })
  # 
  # label <- unlist(apply(dis.mat, 1, function(x){which(x == max(x))}))

  lookup_table <- sort(rownames(cnv.clone))
  names(lookup_table) <- as.character(seq(1, length(lookup_table), 1))
  
  # Replace numeric values with characters using the lookup table
  label <- lookup_table[as.character(label)]
  names(label) <- colnames(rna.norm)
  
  return(label)
}

Assignment_convergence <- function(dna.label, cnv.gene, rna.norm, rna.label, cnv.clone){
  
  item = 0
  mismatch_rate = 1
  start_time <- Sys.time()
  
  while(mismatch_rate > 0.05 && item <= 2){
    item = item + 1
    clone.ids <- unique(sort(dna.label))
    assign.ids <- unique(sort(rna.label))
    all_rna.dis = matrix(0, nrow = length(clone.ids),
                         ncol = dim(rna.norm)[2])
    
    rownames(all_rna.dis) = clone.ids
    colnames(all_rna.dis) = colnames(rna.norm)
    
    for (j in 1:length(assign.ids)){
      index1 = which(rna.label == assign.ids[j])
      tmp_rna = rna.norm[, index1]
      tmp_rna.dis = matrix(0, nrow = length(clone.ids), ncol = length(index1))
      rownames(tmp_rna.dis) = clone.ids
      colnames(tmp_rna.dis) = colnames(tmp_rna)
      
      for (i in 1:length(clone.ids)){
        index2 = which(dna.label == clone.ids[i])
        tmp_cnv = cnv.gene[, index2]
        tmp_Y_proj = Projection(tmp_rna, tmp_cnv)
        tmp_dis = apply(tmp_Y_proj, 1, function(col2) {
          dist(rbind(cnv.clone[i, ],col2), method = "Euclidean")
          #cor(cnv.clone[i, ], col2, method = "pearson")
          
        })
        tmp_rna.dis[i, ] = tmp_dis
      }
      all_rna.dis[, index1] = tmp_rna.dis
    }
    
    label <- apply(all_rna.dis, 2, function(x){which(x == max(x))[1]})
    lookup_table <- clone.ids
    names(lookup_table) <- as.character(seq(1, length(lookup_table), 1))
    
    # Replace numeric values with characters using the lookup table
    label <- lookup_table[as.character(label)]
    names(label) <- colnames(rna.norm)
    
    mismatch_count <- sum(label != rna.label)
    mismatch_rate <- mismatch_count/length(rna.label)
    ## Update the RNA assignment
    rna.label <- label
  }
  
  print(paste0("Iteative:", item))
  
  end_time <- Sys.time()
  execution_time <- end_time - start_time
  print(execution_time)
  print(mismatch_rate)
  
  return(label)
}

############# input #############################
rna.norm <- readRDS('RNA_normalized_matrix.rds')
cnv.gene <- readRDS('CNV_tumor_matrix_gene.rds')
cnv.clone <- readRDS('CNV_subclones.rds')
load('seurat_scDNA.Robj')
seu_epi <- readRDS('seurat_epi.rds')
dna.label <- seurat_scDNA$subclone[which(seurat_scDNA$celltype == 'tumor')]

#######################remove cells #######################
#### calculate the correlation matrix between rna & dna
###############Choose the signal chromosomes#################
gene_locs <- read.table('../gene_locs.sorted.bed')
signal <- fread('par_chr.txt')
gene_locs_sub <- gene_locs %>% filter(V2 %in% signal$V1)
common <- intersect(rownames(rna.norm), gene_locs_sub$V1)
common <- intersect(common, colnames(cnv.clone))

# ###############################################################
L <- Preproc(cnv.gene, rna.norm[common, ], cnv.clone)
rna.common <- L[[1]]
dna.common <- L[[2]]
cnv.clone.common <- L[[3]]


Y_proj <- Projection(rna.common, dna.common)

###############Choose the signal chromosomes#################
gene_locs <- read.table('../gene_locs.sorted.bed')
signal <- fread('par_chr.txt')
gene_locs_sub <- gene_locs %>% filter(V2 %in% signal$V1)
common <- intersect(colnames(cnv.clone.common),
                    gene_locs_sub$V1)
colnames(Y_proj) <- colnames(cnv.clone.common)
rownames(cnv.clone.common) <- rownames(cnv.clone)

label <- Cluster_Alignment(cnv.clone.common[, common], Y_proj[, common])

#new label
label.new <- seu_epi$condition
label.new[which(label.new == 'tumor')] <- label

seu_epi$pro_alignment <- label.new
saveRDS(seu_epi, file = 'seurat_epi.rds')

df <- as.data.frame(table(label.new)) 
write.csv(df, file = 'project_alignment.csv')

