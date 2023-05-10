seurat_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    inputmatrix<-read.csv(inputmatrix,header=TRUE,row.names=1)
    seuratout <- NULL
    cso <- CreateSeuratObject(counts =inputmatrix, project = "cso")
    seuratout<- NormalizeData(cso, normalization.method = "LogNormalize", scale.factor = 10000,verbose=FALSE)
    seuratout <- FindVariableFeatures(seuratout, selection.method = "vst", nfeatures = 2000,verbose=FALSE)
    all.genes <- rownames(seuratout)
    seuratout <- ScaleData(seuratout, features = all.genes,verbose=FALSE)
    seuratout<- RunPCA(seuratout, features = VariableFeatures(object = seuratout),verbose=FALSE)
    seuratout <- FindNeighbors(seuratout, dims = 1:10,verbose=FALSE)
    for (i in 1:100){
    seuratout <- FindClusters(seuratout,resolution=i*0.02,verbose=FALSE)
    if (length(unique(seuratout@active.ident))==c_nums){
        break
        }
    if (i==100){
        seuratout <- FindClusters(seuratout,verbose=FALSE)
        }
    }
    # seuratout <- FindClusters(seuratout)
    if(save.results == TRUE){
        saveRDS(seuratout, file = "seurat_SCluster.rds")
    }
    seuratout<-as.data.frame(seuratout@active.ident)
    return(seuratout)
}

suppressPackageStartupMessages(library(Seurat))

args <- commandArgs(trailingOnly = TRUE)

clusterpath<-args[1]
c_nums<-args[2]
inputmatrix<-args[3]
save.results<-args[4]

if(is.na(save.results)){
    save.results<-FALSE
}

seuratout<-seurat_SCluster(as.integer(c_nums),inputmatrix,save.results)
write.csv(seuratout, file =clusterpath,quote = TRUE)