#Custom R code for Azcorra & Gaertner et al 2023
#R scripts by Zack Gaertner 2022
#All other analyses described in the paper utilize standard seurat functions; an example of the standard seurat integration/clustering workflow with used settings is presented below


library(Seurat)
library(dplyr)
library(SeuratWrappers)
library(SeuratData)
library(tidyverse)

#STANDARD SEURAT WORKFLOW WITH SETTINGS USED IN MANUSCRIPT
#
#
#build seurat objects
f2 <- Read10X(data.dir =  "your-file-here")
f2.seurat <- CreateSeuratObject(counts = f2, project = "Female Mice", min.cells = 3, min.features = 100)
m2 <- Read10X(data.dir =  "your-file-here")
m2.seurat <- CreateSeuratObject(counts = m2, project = "Male Mice", min.cells = 3, min.features = 100)
f2.seurat[["source"]] <- "Sample F2"
m2.seurat[["source"]] <- "Sample M2"
f2.seurat[["percent.mt"]] <- PercentageFeatureSet(object = f2.seurat, pattern = "^mt-")
m2.seurat[["percent.mt"]] <- PercentageFeatureSet(object = m2.seurat, pattern = "^mt-")
#make qc plots and subset based on them
VlnPlot(f2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(m2.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plotf2 <- FeatureScatter(f2.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotf2
f2.seurat <- subset(f2.seurat, subset = nFeature_RNA > 500)
plotm2 <- FeatureScatter(m2.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plotm2
m2.seurat <- subset(m2.seurat, subset = nFeature_RNA > 500)
#integrate the datasets
f2.seurat <- SCTransform(object = f2.seurat, vst.flavor = "v2", method = "glmGamPoi", return.only.var.genes = FALSE, vars.to.regress = "percent.mt")
m2.seurat <- SCTransform(object = m2.seurat, vst.flavor = "v2", method = "glmGamPoi", return.only.var.genes = FALSE, vars.to.regress = "percent.mt")
#integrate objects
hq.Da.list <- list(f2.seurat, m2.seurat)
features <- SelectIntegrationFeatures(object.list = hq.Da.list, nfeatures = 3000)
hq.Da.list <- PrepSCTIntegration(object.list = hq.Da.list, anchor.features = features)
Da.anchors <- FindIntegrationAnchors(object.list = hq.Da.list, normalization.method = "SCT", anchor.features = features)
combined.sct <- IntegrateData(anchorset = Da.anchors, normalization.method = "SCT")
combined.sct <- RunUMAP(combined.sct, reduction = "pca", dims = 1:30)
combined.sct <- FindNeighbors(combined.sct, dims = 1:30)
combined.sct <- FindClusters(combined.sct, resolution = 0.5)


#ADDITIONAL CUSTOM R SCRIPTS FOR ANALYSES
#
#
#Bootstrapping the difference in correlations of gene expression between two pairs of clusters
bootstrapCorDiff <- function(object, ident1, ident2, ident3, ident4, reps = 1000, histo = FALSE){
  message("Note: this is only implemented for 2 pairs of unique clusters (e.g. 1&2 vs 3&4, NOT 1&2 vs 1&3. DO NOT REPEAT CLUSTERS")
  set.seed(1)
  DefaultAssay(object) <- "RNA"
  subset <- subset(object, idents = c(ident1, ident2, ident3, ident4), seed = NULL)
  levels(subset) <- c(ident1, ident2, ident3, ident4)
  avg.subset <- as.data.frame(log1p(AverageExpression(subset, verbose = FALSE)$RNA))
  avg.subset$gene <- rownames(avg.subset)
  colnames(avg.subset) <- c("Ident1", "Ident2", "Ident3", "Ident4")
  cor1 <- rep(0, reps)
  cor2 <- rep(0, reps)
  for(i in 1:reps){
    samp <- sample(1:nrow(avg.subset), nrow(avg.subset), TRUE)
    cor1[i] <- cor(avg.subset$Ident1[samp], avg.subset$Ident2[samp])
    cor2[i] <- cor(avg.subset$Ident3[samp], avg.subset$Ident4[samp])
    corDiff <- cor1 - cor2
  }
  if(histo == TRUE){
    hist(corDiff)
    abline(v = quantile(corDiff,probs=c(0.025)), col="red")
    abline(v = quantile(corDiff,probs=c(0.975)), col="red")
    abline(v = mean(corDiff), col="blue")
  }
  return(corDiff)
}

#making scatter plots for correlation of gene expression for a pair of clusters
library(ggplot2)
library(cowplot)
makeCorScatter <- function(object, ident1, ident2, title = NULL, label.genes = NULL, remove.zero = FALSE){
  theme_set(theme_cowplot())
  DefaultAssay(object) <- "RNA"
  subset <- subset(object, idents = c(ident1, ident2))
  levels(subset) <- c(ident1, ident2)
  avg.subset <- as.data.frame(log1p(AverageExpression(subset, verbose = FALSE)$RNA))
  avg.subset$gene <- rownames(avg.subset)
  colnames(avg.subset) <- c("Ident1", "Ident2")
  message("Number of genes")
  message(nrow(avg.subset))
  if (remove.zero == TRUE) {
    avg.subset <- subset(avg.subset, (avg.subset$Ident1 + avg.subset$Ident2 > 0.15))
    message("Number of genes after removing zeros")
    message(nrow(avg.subset))
  }
  genes.to.label = label.genes
  cor <- cor(avg.subset$Ident1, avg.subset$Ident2)
  p1 <- ggplot(avg.subset, aes(x = Ident1, y = Ident2)) + geom_point() + ggtitle(as.character(title), subtitle = paste0("Correlation coefficient = ", cor)) + xlab(paste0("Cluster #", ident1)) + ylab(paste0("Cluster #", ident2))
  if (is.null(label.genes)) {
    return(p1)
  } else {
    p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE, size = 5)
    p1
  }
}


#Cluster stability calculations
#Several helper functions are needed to execute code at bottom
#
#set downsample size
sample.size <- as.integer(0.8 * ncol(combined.sct))
#helper function for a single run of downsampling and calculating max jaccard for each old cluster
ClusterStability.singlerun <- function(object){
  object.downsample <- subset(object, cells = Clusterstability.cellnames(object), seed = NULL)
  message(head(cellnames))
  object.downsample <- RunPCA(object.downsample, verbose = F, seed.use = NULL)
  object.downsample <- RunUMAP(object.downsample, reduction = "pca", dims = 1:30, verbose = F, seed.use = NULL)
  object.downsample <- FindNeighbors(object.downsample, dims = 1:30, verbose = F)
  object.downsample <- FindClusters(object.downsample, resolution = 0.5, verbose = F)
  nclusters.old <- length(levels(object))
  nclusters.new <- length(levels(object.downsample))
  jac.max <- NULL
  jac.max <- vector(mode = "numeric", length = nclusters.old)
  for (i in 1:nclusters.old) {
    for (z in 1:nclusters.new) {
      cellNames.old <- WhichCells(object, idents = (i-1), seed = NULL)
      cellNames.New <- WhichCells(object.downsample, idents = (z-1), seed = NULL)
      jac.max[i] <- max(jac.max[i], calculateJac(cellNames.old, cellNames.New))
    }
  }
  jac.max
}
#helper function for calculating jaccard
calculateJac <- function(set1, set2){
  in_length = length(intersect(set1, set2))
  un_length = length(union(set1, set2))
  jaccard = in_length/un_length
  return(jaccard)
}
#helper function for getting cell names that are used in downsampled dataset
Clusterstability.cellnames <- function(object){
  allCells <- Cells(object, seed = NULL)
  cellnames <<- sample(allCells, size = sample.size, replace = F)
  cellnames
}
#wrapper function for calculating jaccard n times and saving it as a dataframe
Clusterstability.replicated <- function(object, reps){
  message("Running replication number 1")
  outputs <- data.frame(ClusterStability.singlerun(object))
  rownames(outputs) <- levels(object)
  for (i in 1:(reps - 1)) {
    message("Running replication number ", i+1)
    outputs <- cbind(outputs, ClusterStability.singlerun(object))
  }
  colnames(outputs) <- 1:reps
  outputs
}
#execute code for our dataset (here, it was named "combined.sct")
#can replace the dataset here as needed
#clusters are reassigned numbers to match those listed in the manuscript
stability <- Clusterstability.replicated(combined.sct, 100)
stability.norm <- (stability / 0.8)
medians <- rowMedians(as.matrix(stability))
medians.norm <- (medians / 0.8)
cluster_1 <- as.numeric(stability.norm[4,])
cluster_2 <- as.numeric(stability.norm[5,])
cluster_3 <- as.numeric(stability.norm[1,])
cluster_4 <- as.numeric(stability.norm[3,])
cluster_5 <- as.numeric(stability.norm[12,])
cluster_6 <- as.numeric(stability.norm[2,])
cluster_7 <- as.numeric(stability.norm[10,])
cluster_8 <- as.numeric(stability.norm[11,])
cluster_9 <- as.numeric(stability.norm[9,])
cluster_10 <- as.numeric(stability.norm[6,])
cluster_11 <- as.numeric(stability.norm[7,])
cluster_12 <- as.numeric(stability.norm[8,])
cluster_13 <- as.numeric(stability.norm[13,])
cluster_14 <- as.numeric(stability.norm[14,])
cluster_15 <- as.numeric(stability.norm[15,])
list_stabilities <- list(cluster_1, cluster_2, cluster_3, cluster_4, cluster_5, cluster_6, cluster_7, cluster_8, cluster_9, cluster_10, cluster_11, cluster_12, cluster_13, cluster_14, cluster_15)
boxplot(list_stabilities, main = "Cluster Stability", xlab ="Cluster Number", ylab = "Stability (Normalized Jaccard similarity index)", cex.main = 2, cex.lab = 1.5)

