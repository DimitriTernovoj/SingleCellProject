library(SingleCellExperiment)
library(MAST)
library(Seurat)
library(scater)
library(rhdf5)

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

#function for differential testing
Differential_Testing <- function(sca,name,pos) {

#subsetting the object
sca_ent <- subset(sca, with(colData(sca), startsWith(colData(sca)[,ncol(colData(sca))-2],name)))

print("Dimensions before subsetting:")
print(dim(sca_ent))

sca_ent_filt = sca_ent[rowSums(assay(sca_ent)) != 0, ]

print("Dimensions after subsetting:")
print(dim(sca_ent_filt))

#choosing correct formula; checking if user specified a variable column in the units-file; define & run a hurdle model
if (pos != "") {
zlmCond_ent <- zlm(formula = as.formula(paste("~contrast +",colnames(colData(sce))[pos],"+ nFeatures_RNA")), sca=sca_ent_filt)
} else {
zlmCond_ent <- zlm(formula = ~contrast + nFeatures_RNA, sca=sca_ent_filt)
print("second case")
}

summaryCond_ent <- summary(zlmCond_ent, doLRT='contrastB')

summaryDt_ent <- summaryCond_ent$datatable


result_ent <- merge(summaryDt_ent[contrast=='contrastB' & component=='H',.(primerid, `Pr(>Chisq)`)], #P-vals
                 summaryDt_ent[contrast=='contrastB' & component=='logFC', .(primerid, coef)],
                 by='primerid') #logFC coefficients

#correct for multiple testing (FDR correction) and filtering
result_ent[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
ent_de = result_ent[result_ent$FDR<as.numeric(snakemake@config[["differential_testing"]][["DE_threshold"]]),, drop=F]
ent_de = ent_de[order(ent_de$FDR),]

print("hier kommt DE_threshold")
print(as.numeric(snakemake@config[["differential_testing"]][["DE_threshold"]]))
#ent_all = result_ent[order(result_ent$FDR),]

name_of_csv = paste("results/",name,"_diff_testing.csv", sep="")
write.csv(ent_de, name_of_csv, row.names=FALSE)

#name_of_csv2 = paste("results/",name,"_diff_testing.csv", sep="")
#write.csv(ent_all,name_of_csv2, row.names=FALSE)

print("done mit einem Durchlauf")
}



#read in celltypes used for dt
celltypes <- gsub("-",",",unlist(strsplit(snakemake@config[["differential_testing"]][["clusters"]],",")))

#convert anndata object to SingleCellExperiment
sce <- ReadH5AD_2(unlist(snakemake@input[1]))
saveRDS(sce, "adata_dt.rds")
temp <- readRDS(file="adata_dt.rds")
sce <- as.SingleCellExperiment(temp)

#convert SingleCellExperiment to SingleCellAssay type as required by MAST
sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")

#scale Gene detection rate
colData(sca)$nFeatures_RNA = scale(colData(sca)$nFeatures_RNA)

#finding the position of the user variable if available
if (names(colData(sca))[4] == "nCount_RNA") {
  pos <- ""
} else if (names(colData(sca))[1] != "contrast") {
  pos <- 1
} else if (names(colData(sca))[2] != "region") {
  pos <- 2
} else if (names(colData(sca))[3] != "sample") {
  pos <- 3
} else {
  pos <- 4
}

#dt for all specified clusters (if no clusters are specifed, for all clusters)
if (length(celltypes)!= 0){
   for (i in celltypes){
   print(i)
   Differential_Testing(sca,i,pos)
  }
} else {
    for (i in unique(colData(sca)[,ncol(colData(sca))-2])){
    print(i)
    Differential_Testing(sca,i, pos)
  }
}
