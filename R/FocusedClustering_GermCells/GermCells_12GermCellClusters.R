### R script to generate 12 ordered germ cell clusters in Fall 2017 by Qianyi
### Related to Figure 2 and S2A, and Table S3A: germ cells (N=20,646) with >1k detected genes 

### load file path and libraries
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
library(Seurat)
library(dplyr)
library(Matrix)

######### Extract all germ cells (24 germ cell clusters) with >1k Genes
### load Seurat object with 33k filtered cells from 24 batches (without INT6) and genes with >20 total UMIs and expressed in >15 cells
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge20=dge   # Gene Filter 1 (N=24482): genes >20totalUMIsOverAllCells & >15Cells +31 genes

table(dge@data.info[,52:53]) 
# Cluster31SeriationOLO_GeneFilter1 is the reordered cluster ID 1-31
# Cluster 8-31 are germ cells

### Extract cells in 24 germ cell clusters
cells.use=rownames(dge@data.info)[which(dge@data.info$Cluster31SeriationOLO_GeneFilter1>=8)]

anno <- data.frame(V1=rownames(dge@data.info),V2=dge@data.info$Cluster31SeriationOLO_GeneFilter1,dge@data.info)[which(dge@data.info$Cluster31SeriationOLO_GeneFilter1>=8),]
table(anno$V2)
summary(anno$nGene)
summary(anno$nUMI)
summary(anno$nUMIperCell2)
anno=anno[,-c(10:48,65:66)]

### change GermClusterID 8-31 to 1-24
for(i in 1:24){
    anno$V2[which(anno$Cluster31SeriationOLO_GeneFilter1 == i+7)] <- i
}
anno$V2=factor(anno$V2,levels=1:24)
table(anno$V2)

### extract gene expression matrix for 24 germ cell clusters 
library(Seurat)
dge=SetAllIdent(dge20,id="Cluster31SeriationOLO_GeneFilter1")
table(dge@ident)
dge=SubsetData(dge,ident.use=8:31)
dim(dge@data)       # [1] 24482 29552
dim(dge@data.info)  # [1] 29552    64

dge=SubsetData(dge,cells.use=rownames(anno))
table(dge@ident)

germclusters=anno$V2
names(germclusters)=anno$V1
dge=AddMetaData(dge,germclusters,"GermClusters")
dge=SetAllIdent(dge,id="GermClusters")
table(dge@ident)

###### Remove cells with <=1000 Genes, because germ cells tend to have more transcripts
# note: these are already filtered by keeping cells with >1k genes
# so I do not need to filter by >1k UMIs
anno=anno[anno$nGene > 1000,]
table(anno$V2)
dim(anno)  # 20646 cells

dge=SubsetData(dge24,cells.use=rownames(anno))
table(dge@ident)

### genes expressed in the germ cells with >1k genes subset 
nCellperGene <- rowSums(as.matrix(dge@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # 24475
dge@data=dge@data[genes.use,]
print(table(dge@ident))
print(dim(dge@data))

### Re-Scale data!!! important because PCA uses scale.data
dge=ScaleData(dge)

### PCA for all germ cells with >1k genes using all genes
print(Sys.time())
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())

dgefile="figAug2017_MouseAdultST24_ReCluster24GermClustersLargeCells1kgenes/"
save(dge,file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled.Robj"))
dge24=dge

dge=SetAllIdent(dge,id="GermClusters")
myBrewerPalette=rep("grey60",24)
plot1=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=T)
plot2=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=T)
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster24GermClustersLargeCells_grey.pdf",sep=""),height=7,width=16)
MultiPlotList(list(plot1,plot2),cols = 2)
dev.off()
   

###### 09.23.2017 Highly Variable Genes and PCA
pdf(paste(dgefile,"dge_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
points(data.x[pass.cutoff],data.norm.y[pass.cutoff],col="green3",pch=20,cex=0.6)
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
 
length(dge@var.genes)
# 2047 y.cutoff=0.2

###### PCA for all germ cells with >1k genes using Highly Variable Genes
dgefile="figAug2017_MouseAdultST24_ReCluster24GermClustersLargeCells1kgenes/HVG_"
dge=SetAllIdent(dge,id="GermClusters")

print(Sys.time())
dge <- PCA(dge, pc.genes = dge@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dge,file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))
dge24HVG=dge

plot1=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=T)
plot2=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=T)
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster24GermClustersLargeCells.pdf",sep=""),height=7,width=16)
MultiPlotList(list(plot1,plot2),cols = 2)
dev.off()
   

### Scree Plot for PCA for all germ cells with >1k genes using HVG
numPCs=24 # HVG
i=1
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:200],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,5),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
     
###### Louvain-Jaccard clustering and tSNE using top PCs
j=1
print(j)
dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,2,0.1), print.output = 0, save.SNN = T)

dge24HVG=dge
save(dge,file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))


######### Subset reclustering for Non-SPG germ cells with >1k genes
dgefile="figOct2017_MouseAdultST24_ReCluster19GermClustersLargeCells1kgenes_RemovedSPG1-5/HVG_"
### 19 Germ Clusters: Remove SPG Cluster 1-5 (SPG123+Transition45)
dge=dge24
dge=SetAllIdent(dge,id="GermClusters")
table(dge@ident)
dge=SubsetData(dge,ident.use=6:24)
dge19=dge

###### PCA for non-SPG Germ Cell Cluster with >1k genes Subset
### genes expressed in the subset 
nCellperGene <- rowSums(as.matrix(dge@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
length(genes.use) # 24460
dge@data=dge@data[genes.use,]
print(table(dge@ident))
print(dim(dge@data)) # [1] 24460 18450
dge19=dge

### Re-Scale data!!! important because PCA uses scale.data
dge=ScaleData(dge)
save(dge,file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled.Robj"))
dge19=dge

### Highly Variable Genes and PCA
pdf(paste(dgefile,"dge_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
points(data.x[pass.cutoff],data.norm.y[pass.cutoff],col="green3",pch=20,cex=0.6)
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"HVG"))
dev.off()
 
length(dge@var.genes)
# 1879 y.cutoff=0.2

## PCA for non-SPG germ cells with >1k genes using Highly Variable Genes
print(Sys.time())
dge <- PCA(dge, pc.genes = dge@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dge,file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))
dge19HVG=dge

dge=dge19HVG
myBrewerPalette=rep("gray60",24)
plot1=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=T)
plot2=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=T)
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster19GermClustersLargeCells_RemovedSPG1-5_HVG_grey.pdf",sep=""),height=7,width=16)
MultiPlotList(list(plot1,plot2),cols = 2)
dev.off()
   

### Scree Plot for single PCA
numPCs=20
i=1
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:200],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,5),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.3,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
     


###### Louvain-Jaccard clustering and tSNE using top PCs
j=1
print(j)
dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,2,0.1), print.output = 0, save.SNN = T)
### used 9 clusters with resolution = 0.3
dge=SetAllIdent(dge,id="res.0.3")
PCAPlot(dge,pt.size=1,do.label=T)
PCAPlot(dge,1,3,pt.size=1,do.label=T)
   
dge19HVG=dge
save(dge,file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))

###### order the 9 non-SPG germ cell clusters using optimal leaf reordering by Seriation
### randomly shuffling cells within each cluster
res="res.0.3"
resi=1
dge=SetAllIdent(dge,id=res[resi])
ident=dge@ident
levels=levels(ident)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroid using dissimilarity matrix
library(seriation)
n=ncluster=length(levels)
pdf(file=paste0(dgename,"Centroid_norm_Seriation_Dissimilarity_",res[resi],".pdf"))
tmp=genecountsall[,levels][,1:(sum(ncluster))]
tmp=tmp
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
# note: dist calculate distance between each row
length(da) 
da
 # plot with seriation
 hmap(da) # default method="OLO"
dev.off()
levelss=get_order(seriate(da,method="OLO"))
 

### Reordered clusters for all cells
levels=levelss
cells.use=colnames(dge@data)
# random shuffling cells within ordered clusters
ident=factor(ident,levels=levels)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}

### save the reordered ident to dge file
which(unique(sidecol2[,1])!=rev(get_order(do)-1)) # integer(0)
cells.ident.ordered=factor(as.numeric(sidecol2[,2]),ordered=TRUE)
names(cells.ident.ordered)=names(cells.ident)

ordered="res.0.3order"

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@data.info[,ordered]=factor(dge@data.info[,ordered])
dge=SetAllIdent(dge,ordered)
dge19HVG=dge
save(dge, file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))

###### visualize 9 ordered non-SPG germ cell clusters load(file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))
dge19HVG=dge
dgefile=dgename="figOct2017_MouseAdultST24_ReCluster19GermClustersLargeCells1kgenes_RemovedSPG1-5/HVG"
library(RColorBrewer)
myBrewerPalette=brewer.pal(9,"Set1")
dge=SetAllIdent(dge,id="res.0.3order")
dge@ident=dge@ident+3
### flip PC3
dge@pca.obj[[1]]$x[,3]=-dge@pca.obj[[1]]$x[,3]
dge@pca.rot[,3]=-dge@pca.rot[,3]
### use Cluster 1:3 for SPG-Transition1-Transition2, and use cluster 4-12 for new 9 clusters
tmp=dge@ident
tmp=as.numeric(tmp)+3
tmp=as.factor(tmp)
names(tmp)=names(dge@ident)
dge=AddMetaData(dge,tmp,"res.0.3orderReNum")
dge=SetAllIdent(dge,id="res.0.3orderReNum")
dge@pca.obj[[1]]$x[,3]=-dge@pca.obj[[1]]$x[,3]
dge@pca.rot[,3]=-dge@pca.rot[,3]
pdf(paste(dgefile,"dge_tSNE_res.0.3orderReNum.pdf",sep=""),height=6,width=6.5)
TSNEPlot(dge,pt.size=1,do.label=T,label.size=8)
PCAPlot(dge,pt.size=1,do.label=T)
PCAPlot(dge,1,3,pt.size=1,do.label=T)
PCAPlot(dge,1,4,pt.size=1,do.label=T)
PCAPlot(dge,1,5,pt.size=1,do.label=T)
PCAPlot(dge,1,6,pt.size=1,do.label=T)
PCAPlot(dge,1,7,pt.size=1,do.label=T)
dev.off()
save(dge,file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))
dge19HVG=dge
       
     

######### Plot Figure 2A on 10.7.2017
### load PCA for all germ cells with >1k genes using HVG (dge24HVG)
load(file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))
dge24HVG=dge
### Cluster IDs: SPG-Transition1-Transition2-new 9 clusters
# Transtion1&2 are old CLuster 4 and 5 from 24 germ cell clusters (31 clusters solution)
# new 9 clusters are ordered from reclustering of 19 clusters (24 germ cell clusters removing SPG1-5)
dge=dge24HVG
### save cluster ID for SPG-Transition1-Transition2 as cluster 1-2-3
dge=SetAllIdent(dge,id="GermClusters")
tmp1=dge@ident[which(dge@ident %in% c(1:5))]
table(tmp1)
tmp11=as.numeric(tmp1)
names(tmp11)=names(tmp1)
tmp11[which(tmp1 %in% c(1:3))]<-1
tmp11[which(tmp1 == 4)]<-2
tmp11[which(tmp1 == 5)]<-3
table(tmp11)
### load new 9 clusters for non-SPG germ cells with >1k genes (dge19HVG)
load(file=paste0(home,"data_DGE/19GermClusters1kgenes_RemovedSPG1-5_ReScaled_HVG.Robj"))
dge19HVG=dge
dge=dge19HVG
tmp2=dge@ident
tmp22=as.numeric(tmp2)+3
names(tmp22)=names(tmp2)
table(tmp22)
### New 12 Cluster IDs=SPG-Transition1-Transition2-new 9 non-SPG germ cell clusters
tmp=c(tmp11,tmp22)
table(tmp)
tmp=factor(tmp,levels=1:12,ordered=T)

### save new 12 cluster IDs for germ cells
dge=dge24HVG
dge=AddMetaData(dge,tmp,"SPG45new9clusters")
dge=SetAllIdent(dge,id="SPG45new9clusters")
table(dge@ident)
save(dge,file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))
dge24HVG=dge

### color palette
library(RColorBrewer)
myBrewerPalette=brewer.pal(12,"Paired")

### flip PC1 and PC2
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.obj[[1]]$x[,1]=-dge@pca.obj[[1]]$x[,1]
dge@pca.rot[,2]=-dge@pca.rot[,2]
dge@pca.obj[[1]]$x[,2]=-dge@pca.obj[[1]]$x[,2]

### plot each cluster labeled with cluster ID
pdf(paste(dgefile,"dge_tSNE_SPG45new9clusters.pdf",sep=""),height=6,width=6.5)
TSNEPlot(dge,pt.size=1,do.label=T,label.size=8)
PCAPlot(dge,pt.size=1,do.label=T)
PCAPlot(dge,1,3,pt.size=1,do.label=T) # this is Figure 2A
PCAPlot(dge,1,4,pt.size=1,do.label=T)
PCAPlot(dge,1,5,pt.size=1,do.label=T)
PCAPlot(dge,1,6,pt.size=1,do.label=T)
PCAPlot(dge,1,7,pt.size=1,do.label=T)
dev.off()
# save as Figure 2A
  
                            
######## Select markers for 12 germ cell clusters
### 10.10.2017 analyze markers for 12 germ clusters using dge24HVG (SPG-4-5-new9clusters)
dge=dge24HVG
dge=SetAllIdent(dge,id="SPG45new9clusters")

markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)
table(markersall$cluster)[levels(dge@ident)]
write.table(markersall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_12GermClusterTrajectory_pct0.2_diffpct0.2_thresh2fold_10.10.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")
### save as Table S3A
                            

######### Visualize representative known markers for germ cell types in Figure 2C
### 11.28.2017 known germ cell markers 
gene1=c("Zbtb16","Sall4","Sohlh1","Id4","Gfra1","Uchl1","Kit","Stra8","Prdm9")
gene2=c("H2afx","Spo11","Hormad1","Piwil1","Sycp3","Spag6","Mns1","Tbpl1","Acrv1","Tssk1","Tnp1","Prm1","Pgk2","Hspa1l")
genes=c(gene1,gene2)
setname="Markers"

### color scheme
redblue100.alpha<-rgb(read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)

PC=3;object=dge;

for(j in 1:length(genes)){
feature=features.plot=genes[j]
features.plot; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
pch.use = 16; reduction.use = "tsne";
use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE

cells.use <- set.ifnull(cells.use, colnames(object@data))
dim.code <- translate.dim.code(reduction.use)
dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
data.plot <- FetchData(object, dim.codes, cells.use = cells.use)

x1 <- paste(dim.code, dim.1, sep = "")
x2 <- paste(dim.code, dim.2, sep = "")

data.plot$x <- data.plot[, x1]
data.plot$y <- data.plot[, x2]
data.plot$pt.size <- pt.size
data.use <- data.frame(t(FetchData(object, features.plot, cells.use = cells.use,
use.imputed = use.imputed)))
data.gene0 <- na.omit(data.frame(data.use[features.plot, ]))
data.plot$gene <- t(data.gene0)
st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

### plot heatmap for expression of known markers for germ cells
jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

### To maximize red color, use the same color for top 1% of expression for each marker
z<-st6[,6]
top=0.99
top1percent=round(dim(st6)[1]*top,0):dim(st6)[1]
top1n=length(top1percent)
top1=quantile(z,top)
max=diff(range(z))
z=z[which(z<top1)]
jpeg(paste0(dgename,setname,"_",feature,"_top1percent.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
zcolor <- redblue100.alpha[c((z - min(z))/diff(range(z))*100 + 1, rep(100,top1n))]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()
}
# save as Figure 2C


###### mark transition populations in 12 germ cell clusters
dge=dge24HVG
dge=SetAllIdent(dge,id="SPG45new9clusters")
germ=levels(dge@ident)
### transtion (cluster 2,3,7,11) from 12 germ cell clusters with large cells 
transition=c(2,3,7,11)
nonts=germ[-transition]
dge@data.info$SPG45new9Transition=rep(NA,length(dge@data.info$SPG45new9clusters))
levels(dge@data.info$SPG45new9Transition)=c("GermNonTs","Transition")
dge@data.info$SPG45new9Transition[which(dge@data.info$SPG45new9clusters %in% nonts)] <- "GermNonTs"
dge@data.info$SPG45new9Transition[which(dge@data.info$SPG45new9clusters %in% transition)] <- "Transition"
save(dge,file=paste0(home,"data_DGE/24GermClusters1kgenes_ReScaled_HVG.Robj"))
dge24HVG=dge

dge=SetAllIdent(dge,id="SPG45new9Transition")
table(dge@ident) 
# GermNonTs Transition 
#     16903       3743

myBrewerPalette=c("grey70","forestgreen")
pdf(paste0(dgefile,"PCA_GermTransition.pdf"),width=7.5,height=6)
PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F) # this is insert for Figure 2C
TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=F)
dev.off()
# save as insert for Figure 2C


######### Rank correlation for germ cells ordered by 12 germ cell clusters
### order cells by 12 germ cell clusters and randomly shuffle cells within each cluster
ident=dge@ident
levels=levels(dge@ident)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   set.seed(i)
   tmp=cells[which(cells == levels[i])]
   if(length(tmp)>0){
      tmpname=sample(names(tmp),length(tmp),replace=FALSE)
      cells.use=c(cells.use,tmpname)
   }
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

###### Rank correlation for 12 germ cell cluster centroids - Figure 2B
### for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Calculate correlation for each normalized centroid
cc=cor(as.matrix(genecountsall),method="spearman")
dim(cc)
min(cc)

### labeling for axis 
data.use=cc[levels,levels]
colsep.use=cumsum(table(gsub("_.*","",levels)))
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels)))-table(gsub("_.*","",levels))/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cluster")
colnames(clab)=c("Cluster","")

### color scheme
col.use=redblue100

### heatmap of rank correlation for 12 germ cell cluster centroids
pdf(file=paste(dgename,"Centroid_RankedCorrelation_",res[resi],".pdf",sep=""),height=5.5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),RowSideColors=rlab,ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=1,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
# save as Figure 2B

                            
###### Rank correlation for germ cells with >1k genes ordered by 12 germ cell clusters - Figure S2A left panel
### Caculate rank correlation for germ cells with >1k genes
length(dge@var.genes) # 2047
nrho=cor(as.matrix(dge@data[dge@var.genes,]),method="spearman")
dim(nrho) # 20646

### rank correlation for germ cells ordered by 12 germ cell clusters
            data.use2=rho[cells.use,cells.use]
            data.use2=minmax(data.use2,min=disp.min,max=disp.max)

### labeling for axis
            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+15]=levels(cells.ident)

            row.lab2=gsub(".*_","",lab2)

            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+15]=levels(orig.ident)

            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) # draw a line between datasets

sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}
rlab2=rbind(c("white")[as.numeric(sidecol2[,1])],myBrewerPalette[as.numeric(sidecol2[,2])])
clab2=cbind(rlab2[2,],rlab2[1,])
rownames(rlab2)=c("","Cluster")
colnames(clab2)=c("Cluster","")

### modify color scheme
midrange=median(testcor[which(testcor!=1)])
maxrange=max(testcor[which(testcor!=1)])
maxrange
midrange2=maxrange/2
maxrange2=maxrange
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2)

### heatmap of rank correlation for germ cells ordered by 12 germ cell clusters
jpeg(file=paste(dgename,"RankedCorrelation.jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=2,cexRow=2,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
# save as Figure S2A left panel

                            
###### Jaccard distance for germ cells with >1k genes ordered by 12 germ cell clusters - Figure S2A right panel
### Jaccard distance for germ cells ordered by 12 germ cell clusters
            data.use2=dge@snn.dense[cells.use,cells.use]
            data.use2=minmax(data.use2,min=disp.min,max=disp.max)

### labeling for axis
            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+15]=levels(cells.ident)
            row.lab2=gsub(".*_","",lab2)
            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+15]=levels(orig.ident)
            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) # draw a line between datasets
sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}
rlab2=rbind(c("white")[as.numeric(sidecol2[,1])],myBrewerPalette[as.numeric(sidecol2[,2])])
clab2=cbind(rlab2[2,],rlab2[1,])
rownames(rlab2)=c("","Cluster")
colnames(clab2)=c("Cluster","")

### modify color scheme in order to visualize Jaccard distance
col.use2=redblue100[c(rep(c(41:100),each=10),rep(100,100000))]
length(col.use2)

### heatmap for Jaccard distance of germ cells ordered by 12 germ cell clusters
jpeg(file=paste(dgename,"SNN_",res[resi],".jpeg",sep=""),height=3000,width=3000,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=2,cexRow=2,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(7,5))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
# save as Figure S2A right panel

     
