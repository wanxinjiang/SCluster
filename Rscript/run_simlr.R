simlr_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    inputmatrix<-read.csv(inputmatrix,header=TRUE,row.names=1)
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
    print('asdasdasdasdasda')
    if(save.results == TRUE){
        saveRDS(simlrout, file = "soup_SCluster.rds")
    }

    simlrout <- as.data.frame(simlrout$y$cluster)
    return(simlrout)
}

library(SIMLR)

args <- commandArgs(trailingOnly = TRUE)

clusterpath<-args[1]
c_nums<-args[2]
inputmatrix<-args[3]
save.results<-args[4]

if(is.na(save.results)){
    save.results<-FALSE
}

write.csv('',file =clusterpath,quote = TRUE)
simlrout<-simlr_SCluster(as.integer(c_nums),inputmatrix,save.results)

write.csv(simlrout, file =clusterpath,quote = TRUE)