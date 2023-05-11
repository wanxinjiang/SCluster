RaceID3_SCluster <- function(c_nums=NULL,inputmatrix,save.results=FALSE){
    inputmatrix<-read.csv(inputmatrix,header=TRUE,row.names=1)
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

library(RaceID)
args <- commandArgs(trailingOnly = TRUE)
clusterpath<-args[1]
c_nums<-args[2]
inputmatrix<-args[3]
save.results<-args[4]

if(is.na(save.results)){
    save.results<-FALSE
}

RaceID3out<-try(RaceID3_SCluster(as.integer(c_nums),inputmatrix,save.results))

if('try-error' %in% class(RaceID3out)){
        RaceID3out<-''
    }
    
write.csv(RaceID3out, file =clusterpath,quote = TRUE)