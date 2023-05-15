###################################################
###################################################
## sc3 ##
sc3_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){

    sc3out<-NULL

    sc3out <- SingleCellExperiment(
    assays = list(
        counts = as.matrix(inputmatrix),
        # logcounts = log2(as.matrix(inputmatrix)+1)
        logcounts=as.matrix(inputmatrix)
    )
    )
    rowData(sc3out)$feature_symbol <-rownames(sc3out)
    sc3out <- sc3out[!duplicated(rowData(sc3out)$feature_symbol), ]

    # sce <- sc3_estimate_k(sce)
    # k <- metadata(sce)$sc3$k_estimation

    sc3out <- sc3(sc3out, ks =c_nums, biology = FALSE,n_cores=120,svm_max =30000)
    
    if(save.results == TRUE){
        saveRDS(sc3out, file = "sc3_SCluster.rds")
    }
    sc3out<-as.data.frame(sc3out@colData)

    return(sc3out)
}

## SOUP ##
soup_SCluster<- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    soupout <- NULL

    soup_counts<-t(as.matrix(inputmatrix))
    # log.expr = log2(scaleRowSums(soup_counts)*(10^6) + 1)
    soupout = SOUP(soup_counts,Ks=c_nums, type="log")

    if(save.results == TRUE){
        saveRDS(soupout, file = "soup_SCluster.rds")
    }

    soupout <- as.data.frame(soupout$major.labels)
    return(soupout)
}

## SIMLR ##
simlr_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){

    simlrout <- NULL

    counts_matrix <- as.matrix(inputmatrix)

    ##Estimate_Number_of_Clusters##
    # NUMC=2:30
    # res_example = SIMLR_Estimate_Number_of_Clusters(counts_matrix,NUMC = NUMC,cores.ratio = 0)
    # NUMC[which.min(res_example$K2)]
    # K=NUMC[which.min(res_example$K1)]
    
    if (dim(counts_matrix)[2] < 500) {
    simlrout <- SIMLR(counts_matrix, c = c_nums, cores.ratio = 0)
        } else {
    simlrout <- SIMLR_Large_Scale(counts_matrix, c = c_nums,kk=10)
    }

    if(save.results == TRUE){
        saveRDS(simlrout, file = "soup_SCluster.rds")
    }

    simlrout <- as.data.frame(simlrout$y$cluster)
    return(simlrout)
}

##CIDR##
cidr_SCluster<-function(c_nums=NULL,inputmatrix,save.results=FALSE){

    cidrout <- NULL
    cidrout <- scDataConstructor(as.matrix(inputmatrix))
    cidrout <- determineDropoutCandidates(cidrout)
    cidrout <- wThreshold(cidrout)
    cidrout <- scDissim(cidrout)
    cidrout <- scPCA(cidrout)
    cidrout <- nPC(cidrout)
    cidrout <- scCluster(cidrout,nCluster=c_nums)


    if(save.results == TRUE){
        saveRDS(cidrout, file = "cidr_SCluster.rds")
    }
    cidrout<- as.data.frame(cidrout@clusters)
    return(cidrout)
}

##SINCERA##
sincera_SCluster<- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    sinceraout <- NULL
    counts=as.data.frame(inputmatrix)
    sinceraout <- construct(exprmatrix=counts, samplevector=colnames(counts))

    sinceraout <- expr.minimum(sinceraout, value=0.01)
    # sinceraout <- batch.analysis(sinceraout, analysis=c("distribution"), min.expression=1)
    sinceraout <- prefilterGenes(sinceraout, pergroup=FALSE, min.expression=1, min.cells=1, min.samples=1)
    sinceraout <- normalization.zscore(sinceraout, pergroup=FALSE)
    sinceraout <- doPCA(sinceraout, genes=NULL, use.fast = T)
    sinceraout <- doTSNE(sinceraout, genes=NULL, dims = 1:5, use.fast = T)
    obj <- cluster.geneSelection(sinceraout, method="specificity")
    sinceraout <- setGenesForClustering(sinceraout,value=getGenesForClustering(obj))
    # use gap statistics to determine the number of clusters
    if (FALSE) {
    x <- getExpression(sinceraout, scaled=T, genes=getGenesForClustering(sinceraout))
    cordist <- function(y) as.dist((1-cor(t(y)))/2)
    hclustForGap <- function(y, k) list(cluster=cutree(hclust(cordist(y), method = "average"),k=k))
    gapstats <- clusGap(t(x),
                        FUN = hclustForGap,
                        K.max = 10,
                        B = 100)
    }

    sinceraout <- cluster.assignment(sinceraout,k=c_nums)

    if(save.results == TRUE){
        saveRDS(sinceraout, file = "sincera_SCluster.rds")
    }
    sinceraout <- as.data.frame(getCellMeta(sinceraout, name="GROUP"))
    return(sinceraout) 
}

##SHARAP##
sharp_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){

    sharpout <- NULL

    counts<-as.data.frame(inputmatrix)
    sharpout<-try(SHARP(counts,N.cluster = c_nums))
    if('try-error' %in% class(sharpout)){
        sharpout=SHARP(counts)
    }
    
    if(save.results == TRUE){
        saveRDS(sharpout, file = "sharp_SCluster.rds")
    }

    sharpout<-as.data.frame(sharpout$pred_clusters)
    return(sharpout)
}

##Seurat##

seurat_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    seuratout <- NULL
    cso <- CreateSeuratObject(counts =inputmatrix, project = "cso")
    seuratout<- NormalizeData(cso, normalization.method = "LogNormalize", scale.factor = 10000)
    seuratout <- FindVariableFeatures(seuratout, selection.method = "vst", nfeatures = 2000)
    all.genes <- rownames(seuratout)
    seuratout <- ScaleData(seuratout, features = all.genes)
    seuratout<- RunPCA(seuratout, features = VariableFeatures(object = seuratout))
    seuratout <- FindNeighbors(seuratout, dims = 1:10)
    for (i in 1:100){
    seuratout <- FindClusters(seuratout,resolution=i*0.02)
    if (length(unique(seuratout@active.ident))==c_nums){
        break
        }
    if (i==100){
        seuratout <- FindClusters(seuratout)
        }
    }
    # seuratout <- FindClusters(seuratout)
    if(save.results == TRUE){
        saveRDS(seuratout, file = "seurat_SCluster.rds")
    }
    seuratout<-as.data.frame(seuratout@active.ident)
    return(seuratout)
}

##RACEID3##
RaceID3_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){

    RaceID3out<- NULL
    RaceID3out <- SCseq(inputmatrix)
    RaceID3out <- filterdata(RaceID3out,mintotal=1)
    RaceID3out <- compdist(RaceID3out,metric="pearson")
    RaceID3out <- clustexp(RaceID3out,cln=c_nums,sat=FALSE)

    if(save.results == TRUE){
        saveRDS(RaceID3out, file = "RaceID3_SCluster.rds")
    }

    RaceID3out<-as.data.frame(RaceID3out@cluster$kpart)
    return(RaceID3out)
}
