library(scran)

#read in data
data_mat_df <- read.csv("file_dir/data_mat.csv", header=FALSE)
data_mat <- data.matrix(data_mat_df)
input_groups <- read.csv("file_dir/input_groups.csv", header=TRUE)

#calculate sizefactors
size_factors <- computeSumFactors(data_mat, clusters=input_groups$groups, min.mean=0.1)

#export sizefactors
write.table(size_factors, file= "file_dir/size_factors.csv", sep= ",", row.names = FALSE, col.names="size_factors")
