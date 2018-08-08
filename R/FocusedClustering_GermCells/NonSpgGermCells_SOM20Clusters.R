### R script to generate 20 ordered non-SPG germ cell clusters using SOM and markers for the 20 clusters in Jan 2018
### Related to Figure S2F and Table S3B 

######### Generate 20 ordered non-SPG germ cell clusters using SOM by Qianyi on 1/3/2018
###### load data for non-SPG germ cells with >1k genes (n=18450 cells)
load(file="data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj")
### use normalized exp with highly-variable genes (N=1879 HVG)
data=as.matrix(dge@data)[dge@var.genes,]
data=t(data) # in order to run SOM on cells, I should put cells as rows
dim(data) # [1]  18450 1879 
###### generate SOM 20 clusters for non-SPG germ cells
library(som)
date() # [1] "Wed Jan  3 16:47:58 2018" [1] "Mon Jan  8 13:40:24 2018"
som201<-som(data,20,1)
date() # [1] "Wed Jan  3 16:56:39 2018" 
### save cluster results
save(som201,file=paste0("data_DGE/SOM-20-1_NonSpgGermCells1kgenes_HVG.Robj"))
som=som201
tmp=cbind(rownames(som$data),som$visual)
write.table(tmp,paste0("data_DGE/SOM-20-1_NonSpgGermCells1kgenes_HVG.txt"),row.names=T,col.names=T,quote=FALSE,sep="\t")
### number of cells per cluster
which(rownames(som$data)!=rownames(dge@data.info)) # integer(0)
som.id<-som$visual
table(som.id$x)
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
#1742  835  821 1412 1033  645  753  775  876  886  862  853  828  884  671  747
#  16   17   18   19
# 677  777  639 1734
###### calculate 20 SOM cluster centroids 
tmp=cbind(rownames(som$data),som$visual)
dataall=t(as.matrix(dge@data))
### for gene expression matrix, order cells by cluster id 1, 2, 3, 4, ...
mousecluster=cbind(tmp$x,dataall)
colnames(mousecluster)[1]="ident"
mousecluster=mousecluster[order(mousecluster[,1]),]
head(mousecluster)[,1:8]
### for each cluster, calculate average expression of each gene
genecounts=matrix(,length(unique(mousecluster[,1])),dim(mousecluster)[2])
colnames(genecounts)=colnames(mousecluster)
genecounts[,1]=unique(mousecluster[,1])
for(i in 1:length(unique(mousecluster[,1]))){
    id=unique(mousecluster[,1])[i]
    genecounts[i,-1]=apply(mousecluster[mousecluster[,1]==id,-1],2,function(x) ExpMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}
### for each cluster, calculate total expression of each gene
genetotalcounts=matrix(,length(unique(mousecluster[,1])),dim(mousecluster)[2])
colnames(genetotalcounts)=colnames(mousecluster)
genetotalcounts[,1]=unique(mousecluster[,1])
for(i in 1:length(unique(mousecluster[,1]))){
    id=unique(mousecluster[,1])[i]
    genetotalcounts[i,-1]=apply(mousecluster[mousecluster[,1]==id,-1],2,function(x) log(sum(exp(as.numeric(x))-1)+1)) # function(x) sum(exp(as.numeric(x))-1)
    print(i)
}
### for each cluster, calculate the fraction of cells positively expressing each gene
# the counts of Non-0 values for each gene
genepct=matrix(,length(unique(mousecluster[,1])),dim(mousecluster)[2])
colnames(genepct)=colnames(mousecluster)
genepct[,1]=unique(mousecluster[,1])
for(i in 1:length(unique(mousecluster[,1]))){
    id=unique(mousecluster[,1])[i]
    genepct[i,-1]=apply(mousecluster[mousecluster[,1]==id,-1],2,function(x) sum(x!=0))/table(mousecluster[,1])[i]
    print(i)
}
### download, paste transpose in Excel
write.table(t(genecounts),paste0(dgefile,"Cells1kgenes_SOM-20-1_genesexpavg.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")
write.table(t(genetotalcounts),paste0(dgefile,"Cells1kgenes_SOM-20-1_genesexptotal.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")
write.table(t(genepct),paste0(dgefile,"Cells1kgenes_SOM-20-1_genesfractionposcells.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")

write.table(t(genecounts),paste0(dgefile,"Cells2kgenes_SOM-20-1_genesexpavg.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")
write.table(t(genetotalcounts),paste0(dgefile,"Cells2kgenes_SOM-20-1_genesexptotal.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")
write.table(t(genepct),paste0(dgefile,"Cells2kgenes_SOM-20-1_genesfractionposcells.txt"),quote=FALSE,row.names=T,col.names=F,sep="\t")

