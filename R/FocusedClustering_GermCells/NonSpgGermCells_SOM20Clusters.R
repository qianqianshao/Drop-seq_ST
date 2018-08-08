### R script to generate 20 ordered non-SPG germ cell clusters using SOM and their markers in Jan 2018
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
write.table(t(genecounts),"Cen20.txt",quote=FALSE,row.names=T,col.names=F,sep="\t")
write.table(t(genetotalcounts),"Cen20_tot.txt",quote=FALSE,row.names=T,col.names=F,sep="\t")

                        
######### Identify markers for SOM-20-1 non-SPG germ cells clusters by Jun in Jan 2018
cen20<-as.matrix(read.delim("Cen20.txt",header=T,row.names=1,sep="\t"))
tot20<-as.matrix(read.delim("Cen20_tot.txt",header=T,row.names=1,sep="\t"))

tot.mean20<-apply(tot20,1,mean)
tot.var20<-apply(tot20,1,var)

plot(tot.mean20,tot.var20/(tot.mean20),cex=0.5,ylim=c(0,3))
plot(g.mean20,g.var20/(g.mean20),cex=0.5,ylim=c(0,3))

filter20<-(tot.mean20>2)&( tot.var20/tot.mean20>0.5)
sum(filter20) #3384

cen20.sorted<-cen20[order(tot.mean20),]
g.mean20<-apply(cen20.sorted,1,mean)
g.var20<-apply(cen20.sorted,1,var)
tot20.sorted<-tot20[order(tot.mean20),]
tot.mean20<-apply(tot20.sorted,1,mean)
tot.var20<-apply(tot20.sorted,1,var)

lowe<-lowess(tot.mean20,tot.var20/tot.mean20, f = 1/10, delta = 0.01)
plot(tot.mean20,tot.var20/tot.mean20,cex=0.5,ylim=c(0,3))
lines(lowe$x,lowe$y,col=3,lwd=2)
#lines(lowe1$x,lowe1$y+0.3,col=4,lwd=2)

filter20<-(tot.mean20>1)&( tot.var20/tot.mean20> lowe$y)#8742
points(tot.mean20[filter20],tot.var20[filter20]/tot.mean20[filter20],cex=0.5,col=2)

### kmeans
cen20.norm<-normalize(cen20.sorted[filter20,])
> max(cen20.norm)
[1] 4.18464
> min(cen20.norm)
[1] -3.912198
> exp(4.2)
[1] 66.68633
> sum(cen20.norm>4)
[1] 24

k10<-kmeans(cen20.norm,10)$cluster

heatmap(cen20.norm[k10==1,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==2,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==3,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==4,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==5,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==6,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==7,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==8,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==9,],col=redblue100, Colv = NA,zlim=c(-4,4))
heatmap(cen20.norm[k10==10,],col=redblue100, Colv = NA,zlim=c(-4,4))
 
### Find specific markers for each of 20 SOM clusters
uni20<-cbind(c(1,rep(0,19)), c(0, 1,rep(0,18)), c(rep(0,2), 1,rep(0,17)), c(rep(0,3), 1,rep(0,16)), c(rep(0,4), 1,rep(0,15)), c(rep(0,5), 1,rep(0,14)), c(rep(0,6), 1,rep(0,13)), c(rep(0,7), 1,rep(0,12)), c(rep(0,8), 1,rep(0,11)), c(rep(0,9), 1,rep(0,10)), c(rep(0,10), 1,rep(0,9)), c(rep(0,11), 1,rep(0,8)), c(rep(0,12), 1,rep(0,7)), c(rep(0,13), 1,rep(0,6)), c(rep(0,14), 1,rep(0,5)), c(rep(0,15), 1,rep(0,4)), c(rep(0,16), 1,rep(0,3)), c(rep(0,17), 1,rep(0,2)), c(rep(0,18), 1,0), c(rep(0,19), 1))

uni20.cor<-matrix(0,8742,20)

date()
for (i in 1:20) {
for (j in 1:8742) {
uni20.cor[j,i]<-cor(cen20.norm[j,],uni20[i,],method="spearman")
}
}
date() #10s
plot(apply(uni20.cor,2,function(x) sum(x>0.45)))
 
heatmap(cen20.norm[uni20.cor[,1]>0.5,],col=redblue100, Colv = NA,zlim=c(-4,4))

tmp<-uni20.cor
for (i in 1:8742) {
tmp[i,]<-rank(uni20.cor[i,])
}

> select<-c(1:8742)[(uni20.cor[,1]>0.43)&(tmp[,1]==20)]
> length(select)
[1] 92
> select<-append(select, c(1:8742)[(uni20.cor[,2]>0.42)&(tmp[,2]==20)])
> length(select)
[1] 113
> select<-append(select, c(1:8742)[(uni20.cor[,3]>0.40)&(tmp[,3]==20)])
> length(select)
[1] 137
> select<-append(select, c(1:8742)[(uni20.cor[,4]>0.40)&(tmp[,4]==20)])
> length(select)
[1] 157
> select<-append(select, c(1:8742)[(uni20.cor[,5]>0.40)&(tmp[,5]==20)])
> length(select)
[1] 171
> select<-append(select, c(1:8742)[(uni20.cor[,6]>0.40)&(tmp[,6]==20)])
> length(select)
[1] 191
> select<-append(select, c(1:8742)[(uni20.cor[,7]>0.43)&(tmp[,7]==20)])
> length(select)
[1] 212
> select<-append(select, c(1:8742)[(uni20.cor[,8]>0.43)&(tmp[,8]==20)])
> length(select)
[1] 237
> select<-append(select, c(1:8742)[(uni20.cor[,9]>0.43)&(tmp[,9]==20)])
> length(select)
[1] 253
> select<-append(select, c(1:8742)[(uni20.cor[,10]>0.40)&(tmp[,10]==20)])
> length(select)
[1] 282
> select<-append(select, c(1:8742)[(uni20.cor[,11]>0.43)&(tmp[,11]==20)])
> length(select)
[1] 300
> select<-append(select, c(1:8742)[(uni20.cor[,12]>0.43)&(tmp[,12]==20)])
> length(select)
[1] 317
> select<-append(select, c(1:8742)[(uni20.cor[,13]>0.43)&(tmp[,13]==20)])
> length(select)
[1] 334
> select<-append(select, c(1:8742)[(uni20.cor[,14]>0.43)&(tmp[,14]==20)])
> length(select)
[1] 358
> select<-append(select, c(1:8742)[(uni20.cor[,15]>0.43)&(tmp[,15]==20)])
> length(select)
[1] 397
> select<-append(select, c(1:8742)[(uni20.cor[,16]>0.43)&(tmp[,16]==20)])
> length(select)
[1] 428
> select<-append(select, c(1:8742)[(uni20.cor[,17]>0.43)&(tmp[,17]==20)])
> length(select)
[1] 512
> select<-append(select, c(1:8742)[(uni20.cor[,18]>0.43)&(tmp[,18]==20)])
> length(select)
[1] 625
> select<-append(select, c(1:8742)[(uni20.cor[,19]>0.43)&(tmp[,19]==20)])
> length(select)
[1] 684
> select<-append(select, c(1:8742)[(uni20.cor[,20]>0.43)&(tmp[,20]==20)])
> length(select)
[1] 725
image(t(cen20.norm[select,]),col=redblue100, zlim=c(-4.2,4.2),axes=F)
write.table(cen20.norm[select,], "SOM20_725genes",sep="\t")

select<-c(1:8742)[(uni20.cor[,1]>0.47)&(tmp[,1]==20)]
> length(select)
[1] 44
> select<-append(select, c(1:8742)[(uni20.cor[,2]>0.42)&(tmp[,2]==20)])
> length(select)
[1] 65
> select<-append(select, c(1:8742)[(uni20.cor[,3]>0.40)&(tmp[,3]==20)])
> length(select)
[1] 89
> select<-append(select, c(1:8742)[(uni20.cor[,4]>0.40)&(tmp[,4]==20)])
> length(select)
[1] 109
> select<-append(select, c(1:8742)[(uni20.cor[,5]>0.40)&(tmp[,5]==20)])
> length(select)
[1] 123
> select<-append(select, c(1:8742)[(uni20.cor[,6]>0.40)&(tmp[,6]==20)])
> length(select)
[1] 143
> select<-append(select, c(1:8742)[(uni20.cor[,7]>0.43)&(tmp[,7]==20)])
> length(select)
[1] 164
> select<-append(select, c(1:8742)[(uni20.cor[,8]>0.43)&(tmp[,8]==20)])
> length(select)
[1] 189
> select<-append(select, c(1:8742)[(uni20.cor[,9]>0.43)&(tmp[,9]==20)])
> length(select)
[1] 205
> select<-append(select, c(1:8742)[(uni20.cor[,10]>0.40)&(tmp[,10]==20)])
> length(select)#186
[1] 234
> select<-append(select, c(1:8742)[(uni20.cor[,11]>0.43)&(tmp[,11]==20)])
> length(select)
[1] 252
> select<-append(select, c(1:8742)[(uni20.cor[,12]>0.43)&(tmp[,12]==20)])
> length(select)
[1] 269
> select<-append(select, c(1:8742)[(uni20.cor[,13]>0.43)&(tmp[,13]==20)])
> length(select)
[1] 286
> select<-append(select, c(1:8742)[(uni20.cor[,14]>0.43)&(tmp[,14]==20)])
> length(select)
[1] 310
> select<-append(select, c(1:8742)[(uni20.cor[,15]>0.43)&(tmp[,15]==20)])
> length(select)
[1] 349
> select<-append(select, c(1:8742)[(uni20.cor[,16]>0.43)&(tmp[,16]==20)])
> length(select)
[1] 380
> select<-append(select, c(1:8742)[(uni20.cor[,17]>0.47)&(tmp[,17]==20)])
> length(select)
[1] 406
> select<-append(select, c(1:8742)[(uni20.cor[,18]>0.47)&(tmp[,18]==20)])
> length(select)
[1] 447
> select<-append(select, c(1:8742)[(uni20.cor[,19]>0.46)&(tmp[,19]==20)])
> length(select)
[1] 482
> select<-append(select, c(1:8742)[(uni20.cor[,20]>0.46)&(tmp[,20]==20)])
> length(select)
[1] 508
image(t(cen20.norm[select,]),col=redblue100, zlim=c(-4.2,4.2),axes=F)
### save as Figure S2F right panel
write.table(cen20.norm[select,], "SOM20_508genes",sep="\t")
### save as Table S3B
