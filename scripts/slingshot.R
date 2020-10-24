library(gam)
library(slingshot)
library(SingleCellExperiment)
library(RColorBrewer)
library(Seurat)
library(scater)
library(rhdf5)
library(clusterExperiment)
#so
#print(snakemake@config[["general"]][["amount_of_hvgs"]])

#function to convert anndata object in a SingleCellExperiment
ReadH5AD_2 <- function(h5.path) {
  data <- Seurat::ReadH5AD(h5.path)
  uns <- rhdf5::h5read(h5.path, "uns")

  names(uns) <- gsub("_", ".", names(uns)) # Normalize UNS names
  metadata <- lapply(colnames(data@meta.data), function(meta.name) {
    uns.name <- paste0(meta.name, ".categories")
    if (uns.name %in% names(uns)) {
      uns.array <- as.character(uns[[uns.name]])
      meta.index <- as.numeric(data@meta.data[[meta.name]]) + 1
      return(uns.array[meta.index])
    } else {
      return(data@meta.data[[meta.name]])
    }
  })
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  colnames(metadata) <- colnames(data@meta.data)
  rownames(metadata) <- colnames(data)
  data@meta.data <- metadata
  return(data)
}

#converting
sce <- ReadH5AD_2(unlist(snakemake@input[1]))
saveRDS(sce, "adata_ent.rds")
temp <- readRDS(file="adata_ent.rds")
sce <- as.SingleCellExperiment(temp)

#selecting subset containing only hvgs
sce_sub <- sce[rowData(sce)$highly.variable == TRUE]

colData(sce_sub)[,ncol(colData(sce_sub))-1] <- as.factor(colData(sce_sub)[,ncol(colData(sce_sub))-1])

#find position of initial clustering ------- hier raus wahrscheinlich
#pos = match("louvain.r1",names(colData(sce_sub)))
#initial_clustering = names(colData(sce_sub))[pos+1]
#print(initial_clustering)
#--------

#select latest clustering
current_clustering <- names(colData(sce_sub)[ncol(colData(sce_sub))-1])

#plot visualization of the chosen clusters
colour_map = brewer.pal(8, "Set1")
par(xpd=TRUE)
par(mar=c(4.5,5.5,2,11))
jpeg("figures/slingshot_data_vis.jpg")
plot(reducedDims(sce_sub)$PCA[,1], reducedDims(sce_sub)$PCA[,2], col=colour_map[colData(sce_sub)[,ncol(colData(sce_sub))-1]], bty='L', xlab='PC1', ylab='PC2')
legend(x=12, y=12, legend=unique(colData(sce_sub)[,ncol(colData(sce_sub))-1]), fill=colour_map[as.integer(unique(colData(sce_sub)[,ncol(colData(sce_sub))-1]))])
dev.off()

#infer trajectory
start = snakemake@config[["trajectory_inference"]][["trajectory_start"]]
end = snakemake@config[["trajectory_inference"]][["trajectory_ending"]]

if (start != "" & end != "") { 
sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
} else if (start == "" & end  == "") {
start = toString(unique(colData(sce_sub)[current_clustering])[[1]][1])
end = toString(unique(colData(sce_sub)[current_clustering])[[1]][2])
sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
} else if (start == "" & end != "") {
  if (end == toString(unique(colData(sce_sub)[current_clustering])[[1]][1])) {
  start = toString(unique(colData(sce_sub)[current_clustering])[[1]][2])
  sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
  } else {
  start = toString(unique(colData(sce_sub)[current_clustering])[[1]][1])
  sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
  }
} else if (start != "" & end == "") {
  if (start == toString(unique(colData(sce_sub)[current_clustering])[[1]][1])) {
  end = toString(unique(colData(sce_sub)[current_clustering])[[1]][2])
  sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
  } else {
  end = toString(unique(colData(sce_sub)[current_clustering])[[1]][1])
  sce_startend <- slingshot(sce_sub, clusterLabels = current_clustering, reducedDim = 'PCA', start.clus = start, end.clus = end)
  }
}

#plot of lineage 1 (pseudotime1)
jpeg("figures/slingshot_pst_1.jpg")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = colors[cut(sce_startend$slingPseudotime_1,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(sce_startend)$curve1, lwd=2)
dev.off()

#plot of lineage 2 (pseudotime2)
jpeg("figures/slingshot_pst_2.jpg")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = colors[cut(sce_startend$slingPseudotime_2,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(sce_startend)$curve2, lwd=2)
dev.off()

#plot annotated trajectories
jpeg("figures/slingshot_annotated.jpg")
par(xpd=TRUE)	#sce$louvain_r0.53
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = brewer.pal(11,'Set1')[colData(sce_sub)[,ncol(colData(sce_sub))-1]], pch=16, asp = 1, bty='L', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(sce_startend), lwd=2, type='lineages')
legend(x=10, y=20, legend=unique(colData(sce_sub)[,ncol(colData(sce_sub))-1]), fill=brewer.pal(11,'Set1')[as.integer(unique(colData(sce_sub)[,ncol(colData(sce_sub))-1]))])
dev.off()


#Gene Expression Dynamics
#set the pseudotime variable
t <- sce_startend$slingPseudotime_1
#extract the gene expression matrix
Y <- assay(sce_startend)

#fit a GAM with a loess term for pseudotime
#note: This takes a while
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

#select the top 100 most significant genes that change over pseudotime
topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:100]
#topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:as.numeric(snakemake@config[["amount_of_topgenes"]])]
heatdata <- assay(sce_startend)[rownames(assay(sce_startend)) %in% topgenes, 
                        order(t, na.last = NA)]

#scale the data per gene for visualization
heatdata <- t(scale(t(heatdata)))

#trimm z-score scale
heatdata[heatdata > 3] = 3
heatdata[heatdata < -3] = -3

#get cluster assignment
heatclus <- colData(sce_startend)[,ncol(colData(sce_startend))-4][order(t, na.last = NA)]

#set up a clusterExperiment data structure for heatmap visualization
ce <- ClusterExperiment(heatdata, heatclus, transformation = function(x){x})

#plot the heatmap
#pdf("figures/Heatmap.pdf", height = as.numeric(snakemake@config[["height"]]), width = as.numeric(snakemake@config[["width"]]))
pdf("figures/Heatmap.pdf", height = 15, width = 8)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', fontsize=10)
dev.off()


#non-batch-corrected trajectory inference
#read in files
data_mat <- read.csv("file_dir/data_mat3.csv", header=FALSE)
data_mat <- sapply(data_mat, as.numeric)

obs <- read.csv("file_dir/obs_nbc.csv", header=TRUE)
var <- read.csv("file_dir/var_nbc.csv", header=TRUE)
pca <- read.csv("file_dir/pca_nbc.csv", header=FALSE)
pca_mat <- as.matrix(pca)

#set up a SingleCellExperiment
sce <- SingleCellExperiment(
        assays= data_mat,
        colData = obs,
        rowData = var,
        reducedDims = list(PCA = pca_mat)
        )

colData(sce)[,ncol(colData(sce))] <- as.factor(colData(sce)[,ncol(colData(sce))])

#plot visualization of the chosen clusters
colour_map = brewer.pal(8, "Set1")
par(xpd=TRUE)
par(mar=c(4.5,5.5,2,11))
jpeg("figures/slingshot_data_vis_nbc.jpg")
plot(reducedDims(sce)$PCA[,1], reducedDims(sce)$PCA[,2], col=colour_map[colData(sce)[,ncol(colData(sce))]], bty='L', xlab='PC1', ylab='PC2')
legend(x=12, y=12, legend=unique(colData(sce)[,ncol(colData(sce))]), fill=colour_map[as.integer(unique(colData(sce)[,ncol(colData(sce))]))])
dev.off()

#infer trajectory
sce_startend <- slingshot(sce, clusterLabels = colnames(obs)[ncol(obs)], reducedDim = 'PCA', start.clus= start, end.clus= end) 

#plot of lineage 1 (pseudotime1)
jpeg("figures/slingshot_pst_1_nbc.jpg")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = colors[cut(sce_startend$slingPseudotime_1,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(sce_startend)$curve1, lwd=2)
dev.off()

#plot of lineage 2 (pseudotime2)
jpeg("figures/slingshot_pst_2_nbc.jpg")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = colors[cut(sce_startend$slingPseudotime_2,breaks=100)], pch=16, asp = 1, xlab='PC1', ylab='PC2')
lines(slingCurves(sce_startend)$curve2, lwd=2)
dev.off()

#plot annotated trajectories
jpeg("figures/slingshot_annotated_nbc.jpg")
par(xpd=TRUE)   #sce$louvain_r0.53
plot(reducedDims(sce_startend)$PCA[,c(1,2)], col = brewer.pal(11,'Set1')[colData(sce)[,ncol(colData(sce))]], pch=16, asp = 1, bty='L', xlab='PC1', ylab='PC2')
lines(SlingshotDataSet(sce_startend), lwd=2, type='lineages')
legend(x=10, y=20, legend=unique(colData(sce)[,ncol(colData(sce))]), fill=brewer.pal(11,'Set1')[as.integer(unique(colData(sce)[,ncol(colData(sce))]))])
dev.off()












