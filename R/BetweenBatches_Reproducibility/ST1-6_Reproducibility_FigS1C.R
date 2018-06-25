### R script to compare 6 ST batches on 03/22/2017 by Qianyi, related to Figure S1C
### Unbiased Representation of Adult Mouse Seminferous Tubule (ST) - 6 ST batches

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

### 6 ST batches
topcells=read.table("data_DGE/topcellST", stringsAsFactors=F)
topcells # 3 columns: SeqID, BatchName, OrigBatchName
#73587   ST1 T-7wk-A
#73588   ST2 T-7wk-B
#76287   ST3 T-9wk-1
#76288   ST4 T-9wk-2
#73589   ST5 T-9wk-A
#73590   ST6 T-9wk-B
infile=topcells[,1] 
dataset=topcells[,3]


######### PCA and Louvain-Jaccard Clustering for each individual batch of 6 ST batches
dgelist=list()
for(i in 1:6){
name2=dataset[i]
dgename=dgefile=paste("figMarch2017_MouseAdultTestis6/",dataset[i],"_",sep="")

###### load raw gene expression matrix, filter for cells and genes
dgedata=read.table(paste("data_DGE/",infile[i],"/mouse_gene_exon_tagged_cleaned.dge.txt.gz",sep=""),header=T, row.names=1,stringsAsFactors=F)
names(dgedata)=paste(dataset[i],names(dgedata),sep="_")

### Filter for cells (>500 detected genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dgedata.tmp), value = T) # mouse mt-
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
print(length(which(percent.mito>0.1)))
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1]
print(dim(dgedata.tmp)[2])

### Filter for genes 
nUMIperGene <- rowSums(dgedata.tmp)
nCellperGene <- rowSums(dgedata.tmp>0)
thres.nUMIperGene = dim(dgedata.tmp)[2]*0.0017
dgedata2=dgedata.tmp[nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nUMIperGene/2,]
print(dim(dgedata2))

###### plot distribution ordered by number of unique UMI per gene or per cell for each ST batch
countincell0<-apply(dgedata2,2,sum)
countingene0<-apply(dgedata2,1,sum)
dgeorder <- dgedata2[order(countingene0,decreasing=TRUE),order(countincell0,decreasing=TRUE)]
dgeorder=as.matrix(dgeorder)
logdge=log(dgeorder+1)
dim(logdge)
countincell=apply(exp(logdge)-1,2,sum)
countingene=apply(exp(logdge)-1,1,sum)
dgesmooth=as.matrix(logdge)
by=20
dgesmooth=matrix(NA,nrow=floor(dim(logdge)[1]/by),ncol=dim(logdge)[2])
dim(dgesmooth)
c=1
for(i in 1:floor(dim(logdge)[1]/by)){
ooo=logdge[c:(c+by-1),]
dgesmooth[i,]=apply(ooo,2,mean)
c=c+by
print(i)
}
redblue100<-rgb(read.table('data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))
col=redblue100[c(10:60,rep(61:80,each=5),rep(81:100,each=20))]
jpeg(file=paste(dgename,"heatmap2.jpeg",sep=""),res=300,height=3200,width=3200)
par(fig=c(0.15,1,0,0.6),mar=c(2,2,0.5,1.5),mgp=c(2, 0.5, 0))
image(dgesmooth,xaxt="n",yaxt="n",col=col,xlab="",ylab="")
par(fig=c(0.15,1,0.65,1),mar=c(2.8,2,2.5,1.5),mgp=c(1.8, 0.5, 0),new=T)
plot(1:length(countingene),countingene,cex.lab=1.5,ylab="",xlab=paste("Unique UMIs in", length(countingene),"Genes (Heatmap log-transformed & bin 20 genes)"),cex.main=1.2)
mtext(name2, side=3, line=0.8, at=length(countingene)/2,cex=1.8)
mtext(paste("z =",ceiling(max(dgesmooth))), side=1, line=2.8, at=length(countingene),cex=1.3)
par(fig=c(0.15,1,0.6,0.65),mar=c(0,2,1,1.5),mgp=c(2, 0.5, 0),new=T)
zz=seq(0,ceiling(max(dgesmooth)),0.001)
image(cbind(zz,zz),col=col,xaxt="n",yaxt="n",cex.main=1) #
par(fig=c(0,0.15,0,0.6),mar=c(2,3,0.5,0),mgp=c(1.5, 0.5, 0),new=T)
plot(countincell,1:length(countincell),cex.lab=1.5,xlab="",ylab=paste("Unique UMIs in", length(countincell),"Cells"))
dev.off()


###### setup Seurat object, normalize for each cell, log-transform
dge <- new("seurat", raw.data = dgedata2)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jan2017-MouseAdultTestis6",names.field = 1,names.delim = "_")
dge

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

### Violin plot for distribution of nGene, nUMI and %MT for each ST batch
pdf(file=paste(dgefile,"NGeneUMImt.pdf",sep=""),height=5,width=10)
VlnPlot(dge, c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

### Regress out unwanted sources of variation, standardize for each gene
#note that this overwrites dge@scale.data. Therefore, if you intend to use RegressOut, you can set do.scale=F and do.center=F in the original object to save some time.
dge <- RegressOut(dge, latent.vars = c("nUMI", "percent.mito"))


###### Run PCA using all genes 
Sys.time()
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
Sys.time()
save(dge, file = paste0("data_DGE/MouseAdultTestis6_",dataset[i],".Robj"))

### visualize PC1-5
pdf(paste(dgefile,"dge_PCA_all.pdf",sep=""),height=7.5,width=16)
plot1=PCAPlot(dge,1,2,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
plot2=PCAPlot(dge,1,3,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
MultiPlotList(list(plot1,plot2),cols = 2)
plot1=PCAPlot(dge,1,4,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
plot2=PCAPlot(dge,1,5,do.return = TRUE,pt.size = dge@data.info$nUMIperCell2)
MultiPlotList(list(plot1,plot2),cols = 2)
dev.off()


###### Choose top PCs
### Scree Plot
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:120],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
legend("topright",legend=dataset[i],cex=1.5)
dev.off()
  
### density plot of Eigenvalue
numPCs=c(13,13,13,13,13,13)
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue,0.4,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=dataset[i],cex=1.5)
dev.off()
jpeg(paste(dgefile,"dge_PCA_perm100_topPCs",numPCs[i],".jpeg",sep=""),height=3000,width=3000,res=300)
JackStrawPlot( dge,PCs = 1:(numPCs[i]+1),nCol=ceiling(sqrt(numPCs[i]) ) )
dev.off()
  
### Jackstraw Permutation of 1k rounds
library(parallel) # mclapply belongs to package parallel
dge=JackStrawMC(dge,num.pc=40,num.replicate = 1000,do.print = TRUE,num.cores=4)
print(Sys.time())
# save jackstraw permutation to a new dge
save(dge, file = paste("/home/qzm/data_DGE/MouseAdultTestis6_",dataset,"_2.Robj",sep=""))
# plot p-value of top PCs from Jackstraw permuation
jpeg(paste(dgefile,"dge_PCA_perm1k_topPCs13.jpeg",sep=""),height=3000,width=3000,res=300)
JackStrawPlot(dge,PCs = 1:14,nCol=4)
dev.off()
  
### decided to use top 13 PCs for each of the 6 ST batches based on scree plot, density plot of eigenvalues, and Jackstraw permutation


###### Run tSNE using top 13 PCs
numPCs=c(13,13,13,13,13,13)
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)    # max_iter=2000
pdf(paste(dgefile,"tSNE.pdf"),height=7.5,width=9)
TSNEPlot(dge,pt.size = dge@data.info$nUMIperCell2)
PCAPlot(dge,pt.size = dge@data.info$nUMIperCell2)
dev.off()

###### Run Louvain-Jaccard Clustering using top 13 PCs
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.3,3,0.1), print.output = 0, save.SNN = T)
# resolution parameter: sets the granularity of the downstream clustering, with increased values leading to a greater number of clusters. setting this parameter between 0.6-1.2 typically returns good results for single cell datasets of around 3K cells. Optimal resolution often increases for larger datasets.
save(dge,file= paste("/home/qzm/data_DGE/MouseAdultTestis6_",dataset[i],".Robj",sep=""))

### use 14 cluster solutions for each of the 6 ST batches 
res=paste0("res.",c(1.2,1.8,1.2,1.2,1.8,1.2))
resi=i
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])
dge@data.info[,res[i]]=dge@ident
save(dge,file= paste("/home/qzm/data_DGE/MouseAdultTestis6_",dataset[i],".Robj",sep=""))
dgelist[[i]]=dge

###### Order the 14 clusters for each of the 6 ST batches
### for each cluster, calculate cluster centroids of each gene
tmpdge=data.frame(t(as.matrix(dge@data)))
mouseclustersall=data.frame(ident=dge@ident,t(as.matrix(dge@data)))
genecountsall=matrix(,dim(dge@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dge@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Reordering cluster centroids by seriation
library(seriation)
n=ncluster=length(levels)
tmp=genecountsall[,levels]
da <- dist(t(as.matrix(tmp)), method = "euclidean")
pdf(file=paste0(dgename,"Centroid_norm_Seriation_Dissimilarity_",res[resi],".pdf"))
 dissplot(da, method="OLO",options = list(main = paste("24ST Dissimilarity with seriation OLO")))
dev.off()
levelss=get_order(seriate(da,method="OLO"))
levels=levelss

### Order cells with re-ordered cluster IDs
ident=factor(dge@ident,levels=levelss)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   tmpname=names(cells[which(cells == levels[i])])
   cells.use=c(cells.use,tmpname)
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels
sidecol2=cells.ident
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}
cells.ident.ordered=factor(sidecol2[,2],levels=1:ncluster,ordered=TRUE)
names(cells.ident.ordered)=names(cells.ident)

### save the re-ordered cluster IDs
ordered="Clusters14seriation"
dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge=SetAllIdent(dge,ordered)
dgelist[[i]]=dge
save(dge,file= paste("/home/qzm/data_DGE/MouseAdultTestis6_",dataset[i],".Robj",sep=""))
}


###### Visualize 14 ordered clusters in PCA and tSNE spaces
plott=plot2=plot3=plot4=plot5=list()
for(i in 1:6){
dge=dgelist[[i]]
### flip PCs to make the directions consistent among the 6 batches
if(i==5){
dge@pca.rot$PC2=-dge@pca.rot$PC2
dge@pca.rot$PC3=-dge@pca.rot$PC3
}
if(i==2 | i==3 | i==5){
dge@pca.rot$PC4=-dge@pca.rot$PC4
}
if(i==4 | i==5| i==6){
dge@pca.rot$PC5=-dge@pca.rot$PC5
}
### re-define coordinates for DimPlot
plot2[[i]]=PCAPlot(dge,do.return=TRUE,pt.size = dge@data.info$nUMI,do.label=TRUE)
plot3[[i]]=PCAPlot(dge,1,3,do.return=TRUE,pt.size = dge@data.info$nUMI,do.label=TRUE)
plot4[[i]]=PCAPlot(dge,1,4,do.return=TRUE,pt.size = dge@data.info$nUMI,do.label=TRUE)
plot5[[i]]=PCAPlot(dge,1,5,do.return=TRUE,pt.size = dge@data.info$nUMI,do.label=TRUE)
plott[[i]]=TSNEPlot(dge,do.return=TRUE,pt.size = dge@data.info$nUMI,do.label=TRUE,label.size=8)
}
pdf(paste0(dgefile,"cluster_14clusters_seriation.pdf"),width=20,height=12)
MultiPlotList(plot2,cols = 3)
MultiPlotList(plot3,cols = 3)
MultiPlotList(plot4,cols = 3)
MultiPlotList(plot5,cols = 3)
MultiPlotList(plott,cols = 3)
dev.off()


######### Rank Correlation of 14 Cluster Centroids Across Each of the 6 ST Batches
###### Merged Gene Expression Data for 6 ST batches
### read in files in mouse cells and genes
mouse1=read.table(paste(paste0(home,"data_DGE/",infile[1],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
mouse2=read.table(paste(paste0(home,"data_DGE/",infile[2],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
mouse3=read.table(paste(paste0(home,"data_DGE/",infile[3],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
mouse4=read.table(paste(paste0(home,"data_DGE/",infile[4],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
mouse5=read.table(paste(paste0(home,"data_DGE/",infile[5],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
mouse6=read.table(paste(paste0(home,"data_DGE/",infile[6],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz",sep=""),header=T, stringsAsFactors=F)
# add batch name to Cell barcode names
names(mouse1)[-1]=paste(topcells[1,3],names(mouse1)[-1],sep="_")
names(mouse2)[-1]=paste(topcells[2,3],names(mouse2)[-1],sep="_")
names(mouse3)[-1]=paste(topcells[3,3],names(mouse3)[-1],sep="_")
names(mouse4)[-1]=paste(topcells[4,3],names(mouse4)[-1],sep="_")
names(mouse5)[-1]=paste(topcells[5,3],names(mouse5)[-1],sep="_")
names(mouse6)[-1]=paste(topcells[6,3],names(mouse6)[-1],sep="_")
# merge all files 
mousetmp=Reduce(function(x,y) merge(x,y,all=TRUE), list(mouse1,mouse2,mouse3,mouse4,mouse5,mouse6) )
# rm data to save memory
rm(mouse1,mouse2,mouse3,mouse4,mouse5,mouse6)
gc()
# convert back to data.frame, then change NA to 0
Sys.time()
mousetmp[is.na(mousetmp)] <- 0
Sys.time()                   # 14 min passed
# set gene names as row.names
row.names(mousetmp)=mousetmp[,1]
mousetmp[1:5,1:2]
testis6=mousetmp[,-1]

### Filter for cells (>500 detected genes and <10% MT)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dgedata.tmp), value = T) # mouse mt-
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
print(length(which(percent.mito>0.1)))
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1]
print(dim(dgedata.tmp)[2])

### Filter for genes 
nUMIperGene <- rowSums(dgedata.tmp)
nCellperGene <- rowSums(dgedata.tmp>0)
thres.nUMIperGene = 24
dgedata2=dgedata.tmp[nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nUMIperGene/2,]
print(dim(dgedata2)) # 12991 19243


###### setup Seurat object for merged 6 ST batches, normalize for each cell, log-transform
dge <- new("seurat", raw.data = dgedata2)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jan2017-MouseAdultTestis6",names.field = 1,names.delim = "_")
dge

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

### Run PCA using all genes 
Sys.time()
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
Sys.time()
save(dge,file = "data_DGE/MouseAdultTestis6.Robj")


###### Rank correlation of cluster centroids for each of 6 ST batches
load(file = "data_DGE/MouseAdultTestis6.Robj")
dgeall=dge
dgename=dgefile="figMarch2017_MouseAdultTestis6/6datasets_"

### cells ordered first by 6 ST batches, then by 14 cluster IDs of each batch
blockident=NULL
batch=levels(dgeall@data.info$orig.ident)
nbatch=length(batch)
for(bb in 1:nbatch){
  dge=dgelist[[i]]
  blockident=c(blockident,paste(dge@data.info$orig.ident,dge@ident,sep="_")
}

clusters=1:14
levels2=NULL
for(bb in 1:nbatch){
  for(cluster in clusters){
    levels2=c(levels2,paste(batch[bb],cluster,sep="_"))
  }
}
levels2=levels2[which(levels2 %in% unique(blockident))]
levels=levels2

ident=factor(blockident,levels=levels)

cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   tmpname=names(cells[which(cells == levels[i])])
   cells.use=c(cells.use,tmpname)
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### Calculate 14 cluster centroids for each gene in each of the 6 ST batches
tmpdge=data.frame(t(as.matrix(dgeall@data[,cells.use])))
## make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dgeall@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dgeall@data[,cells.use])))
### for each cluster, calculate average normalized expression of each gene
genecountsall=matrix(,dim(dgeall@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dgeall@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in 1:length(unique(mouseclustersall$ident))){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==as.character(unique(mouseclustersall$ident))[i],-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Calculate rank correlation for each normalized centroid
cc=cor(as.matrix(genecountsall),method="spearman")

### labeling and color scheme
data.use=cc[levels,levels]
col.use=redblue100

### Plot Heatmap for rank correlation of 14 cluster centroids across each of the 6 ST batches
# this is Figure S1C

pdf(file=paste(dgename,"Centroid_RankedCorrelation_indiv_14clusters_image.pdf",sep=""),height=5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
image(as.matrix(data.use),col=col.use,cex=1.2,cex.lab=1.3)
dev.off()

# add scale bar for color scale
pdf(file=paste(dgename,"RankedCorrelation_norm_scalebar_centroid.pdf"))
image(cbind(seq(0,1,0.01),seq(0,1,0.01)),col=col.use,cex.axis=2)
dev.off()

