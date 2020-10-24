library(RColorBrewer)
library(monocle)

#read in files
data_mat_mon <- read.csv("file_dir/data_mat_mon.csv", header=FALSE)
data_mat_mon <- sapply(data_mat_mon, as.numeric)

obs_mon <- read.csv("file_dir/obs_mon.csv", header=TRUE)
var_mon <- read.csv("file_dir/var_mon.csv", header=TRUE)

#latest clustering
current = colnames(obs_mon)[ncol(obs_mon)]

#set up the CellDataSet data structure
pd <- AnnotatedDataFrame(data = obs_mon)
fd <- AnnotatedDataFrame(data = var_mon)
colnames(data_mat_mon) <- rownames(pd)
rownames(data_mat_mon) <- rownames(fd)
ie_regions_cds <- newCellDataSet(cellData=data_mat_mon, phenoData=pd, featureData=fd, expressionFamily=negbinomial.size())

#normalize the count data
ie_regions_cds <- estimateSizeFactors(ie_regions_cds)

#calculate dispersions to filter for highly variable genes
ie_regions_cds <- estimateDispersions(ie_regions_cds)

#filter specified clusters
cell_types = as.character(pData(ie_regions_cds)[,ncol(pData(ie_regions_cds))-1])

if (snakemake@config[["trajectory_inference"]][["clusters_to_include"]] != "") {
cell_mask = rep(FALSE, length(cell_types))
cells_to_keep = unlist(strsplit(snakemake@config[["trajectory_inference"]][["clusters_to_include"]],","))

for (item in cells_to_keep) {cell_mask = cell_mask | startsWith(cell_types, item)}
print("Number of cells after filtering:")
print(sum(cell_mask))
print("")
} else {
cell_mask = rep(TRUE, length(cell_types))
}

#filter highly variable genes from our analysis
hvg_mask = fData(ie_regions_cds)$highly_variable
hvg_mask <- sapply(hvg_mask, as.logical)
ie_regions_cds <- ie_regions_cds[hvg_mask, cell_mask]

#do dimensionality reduction
ie_regions_cds <- reduceDimension(ie_regions_cds, norm_method = 'vstExprs', reduction_method='DDRTree', verbose = F, max_components = 7)

#run for the first time to get the ordering
ie_regions_cds <- orderCells(ie_regions_cds)

#find the correct root state the corresponds to the specified cluster
tab1 <- table(pData(ie_regions_cds)$State, pData(ie_regions_cds)[,ncol(pData(ie_regions_cds))-3])

if (snakemake@config[["trajectory_inference"]][["trajectory_start"]] != "") {
id = which(colnames(tab1) == snakemake@config[["trajectory_inference"]][["trajectory_start"]])
} else {
id = which(colnames(tab1) == unique(pData(ie_regions_cds)[,ncol(pData(ie_regions_cds))-3])[1])
}

root_name = names(which.max(tab1[,id]))

#run a second time to get the correct root state that overlaps with the specified cells
ie_regions_cds <- orderCells(ie_regions_cds, root_state=root_name)

#get a nice colour map
custom_colour_map = brewer.pal(20,'Set1')

#visualize trajectory
options(repr.plot.width=6, repr.plot.height=3)
temp <- plot_complex_cell_trajectory(ie_regions_cds, color_by = current, show_branch_points = T, 
                             cell_size = 2, cell_link_size = 1, root_states = c(root_name)) +
scale_size(range = c(0.2, 0.2)) +
theme(legend.position="left", legend.title=element_blank(), legend.text=element_text(size=rel(1.5))) +
guides(colour = guide_legend(override.aes = list(size=6))) + 
scale_color_manual(values = custom_colour_map)

ggsave("monocle_trajectory.jpg", plot = temp, path = "figures", width = 13.33, height = 10.66)

#visualize pseudotime found
options(repr.plot.width=5, repr.plot.height=4)
temp2 <- plot_cell_trajectory(ie_regions_cds, color_by="Pseudotime")

ggsave("monocle_pseudotime.jpg", plot = temp2, path ="figures", width= 8, height = 10.66)
