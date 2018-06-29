### R script for focused analysis of all Somatic Cells in merged 25 ST batches in Dec 2017 by Qianyi
### Related to Figure 1B right panel, Figure S1E and Figure 5A-B

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)
### path
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
dgefile=dgename="figDec2017_MouseAdultST25/Somatic_"

######### Extract all Somatic cells in 25 ST datasets
###### 1. Extract filtered cells and all detected genes for somatic group of 24 ST batches
### load Seurat object for 24 ST batches
load(file =paste0(home,"data_DGE/MouseAdultST24genesall.Robj"))
dgeall # 36979 genes across 33180 samples
dge=dgeall
dge=SetAllIdent(dge,id="celltype2")
dge@ident=factor(dge@ident,levels=levels(dge@data.info$celltype2))
table(dge@ident)
#         Somatic    Spermatogonia     Spermatocyte   RoundSpermatid     Elongating
#            3628             4239             8792             9923            6598 
### extract somatic cells 
sets="Somatic"
setsname="Somatic"
setslabel="Somatic"
dgetmp=SubsetData(dge,ident.use=sets)
print(table(dgetmp@ident))
# Somatic 
#    3628
cells.use=colnames(dgetmp@data)
### extract detected genes for somatic group
nCellperGene <- rowSums(as.matrix(dgetmp@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
### extract gene expression matrix for somatic group
oldsomatic=dgetmp@raw.data[genes.use,cells.use]
oldsomatic=as.matrix(oldsomatic)
print(dim(oldsomatic)) # [1] 27504  3628
### save old 24 ST datasets Somatic cells
write.table(as.matrix(oldsomatic),file=paste0(home,"data_DGE/MouseAdultST24/MouseAdultST24Somatic_3628cells_27504genes.dge.txt"), quote=FALSE, sep='\t',row.names=T, col.names=T)
### read old 24 ST somatic cells
oldsomatic=read.table(file=paste0(home,"data_DGE/MouseAdultST24/MouseAdultST24Somatic_3628cells_27504genes.dge.txt"), row.names=1, header=T)
oldsomatic=data.frame(GENE=rownames(oldsomatic),oldsomatic)

### old 24 datasets somatic cell types
dgetmp=SetAllIdent(dgetmp,id="celltype")
oldsomaticID=factor(dgetmp@ident,levels=levels(dgetmp@data.info$celltype)[1:6])

###### 2. Extract filtered cells and all detected genes for INT6 (somatic cells of new Sca1 dataset)
### load new Sca1 dataset
dgedata <- Read10X(paste0(home,"data_DGE/92783"))
### Filter for cells (keep cells with >500 genes and <10% MT)
nUMIperCell <- colSums(dgedata)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dgedata.tmp), value = T) # mouse
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1]
dgedata = dgedata.tmp
dim(dgedata) # [1] 29786  1990
### load clusters for new Sca1 data
load(file=paste0(home,"data_DGE/NewSca1_Nov2017.Robj"))
dgeSca1=dge
### extract somatic only for new Sca1 dataset (remove germ clusters 1,2,6 from 9 clusters)
dge=dgeSca1
dge=SetAllIdent(dge,id="res.0.5")
dge=SubsetData(dge,ident.use=c(3:5,7:9))
dgeSca1Somatic=dge
cells.use=colnames(dgeSca1Somatic@data)
length(cells.use) # 1453 somatic cells from new Sca1 data
table(dgeSca1Somatic@ident)
#  3   4   5   7   8   9 
# 27  87  62  35 634 608
newSca1somatic=dgedata[,cells.use]
dim(newSca1somatic) # [1] 29786  1453
nCellperGene <- rowSums(newSca1somatic>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
newSca1somatic=newSca1somatic[genes.use,]
dim(newSca1somatic) # [1] 25187  1453
### new Sca1 somatic cell types
newSca1somaticID=rep(NA,length(dgeSca1Somatic@ident))
names(newSca1somaticID)=gsub("NewSca1","INT6",names(dgeSca1Somatic@ident))
### re-name NewSca1 data as INT6
tmp=gsub("NewSca1","INT6",colnames(newSca1somatic))
colnames(newSca1somatic)=tmp
### save new Sca1 somatic cells
write.table(as.matrix(newSca1somatic),file=paste0(home,"data_DGE/92783/INT6Somatic_1453cells_25187genes.dge.txt"), quote=FALSE, sep='\t',row.names=T, col.names=T)
### Read new Sca1 somatic cells
newSca1somatic=read.table(file=paste0(home,"data_DGE/92783/INT6Somatic_1453cells_25187genes.dge.txt"), row.names=1,header=T)
newSca1somatic=data.frame(GENE=rownames(newSca1somatic),newSca1somatic)
dim(newSca1somatic) # [1] 25187  1453+1

###### 3. Merge somatic cells of 24 ST datasets with INT6 (somatic cells from new Sca1 dataset)
somatic=merge(oldsomatic,newSca1somatic,by="GENE",all=TRUE)
dim(oldsomatic)     # 27504, 3628+1
dim(newSca1somatic) # 25187, 1453+1
dim(somatic)        # 29526 genes, 5081 cells + 1 column of gene name
somatic[is.na(somatic)] <- 0  # change NA to 0
row.names(somatic)=somatic[,1]
somatic[1:5,1:2]
somatic=somatic[,-1]
dim(somatic) # 29526 genes, 5081 cells
write.table(somatic, file=paste0(home,"data_DGE/MouseAdultST25/mergedMouseAdultST25Somatic_5081cells_29526genes.dge.txt"), quote=FALSE, sep='\t',row.names=T, col.names=T)
somatic=read.table(paste0(home,"data_DGE/MouseAdultST25/mergedMouseAdultST25Somatic_5081cells_29526genes.dge.txt"), header=T,row.names=1)
dgedata=somatic

###### 4. Filter genes for somatic group of 25 ST datasets
nUMIperGene <- rowSums(dgedata)
nCellperGene <- rowSums(dgedata>0)
dim(dgedata)[2]*0.0012 # 6.
dim(dgedata)[2]*0.0006 # 3.

### keep genes expressed in >3 cells
thres.nUMIperGene = 3
thres.nCellperGene = 3
length(nUMIperGene)
length(which(nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene))    # [1] 22734
genekeep=rownames(dgedata)[nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene]
dgedata2=dgedata[genekeep,]
dim(dgedata2) # 22734, 5081

###### 5. Normalization for cells, log-transformation and standardization for genes
dge <- new("seurat", raw.data = dgedata2)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Dec2017-MouseAdultST25Somatic",names.field = 1,names.delim = "_")
dge
# 22734 genes across 5081 samples.
which(table(gsub("_.*","",colnames(dgedata2))) != table(gsub("_.*","",colnames(dge@data)))) # named integer(0)

# change order of dataset in orig.ident
dataset=unique(gsub("_.*","",colnames(dge@data)))[c(1:16,25,17:24)]
length(dataset) # 25=6(ST)+2(1nDepleted)+3(SPG)+6(INT)+8(SER)
dge@data.info$orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
dge@ident=factor(dge@ident,levels=dataset)

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes) # 28
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

### Select Highly variable genes
pdf(paste(dgefile,"dge_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
length(dge@var.genes) # 2464


###### 6. PCA for somatic subset of merged 25 datasets 
Sys.time() # [1] "2017-12-04 10:56:59 EST"
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
Sys.time() # [1] "2017-12-04 11:33:37 EST"
dgeSomatic=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))

### Scree Plot for single PCA
numPCs=23 
i=1
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:200],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,20,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue # 1.909231
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,6),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
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


###### 7. tSNE and Louvain-Jaccard clustering using top PCs
numPCs=23 
j=1
dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,3,0.1), print.output = 0, save.SNN = T)
dgeSomaticall=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))

### check number of clusters for 25 ST somatic subset
print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)),length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.6)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2)),length(unique(dge@data.info$res.2.1)),length(unique(dge@data.info$res.2.2)),length(unique(dge@data.info$res.2.3)),length(unique(dge@data.info$res.2.4)),length(unique(dge@data.info$res.2.5)),length(unique(dge@data.info$res.2.6)),length(unique(dge@data.info$res.2.7)),length(unique(dge@data.info$res.2.8)),length(unique(dge@data.info$res.2.9)),length(unique(dge@data.info$res.3)) ))

# res.0.5 and res.0.6 both give 12 clusters, with exactly same cells for each cluster
table(dge@data.info$res.0.5,dge@data.info$res.0.6)
# use res.0.5 

### identify differentially-expressed markers
dge=SetAllIdent(dge,id="res.0.5")
markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)

###### Use res.0.5 -> 12 clusters -> merge to 7 major cell types
dge@data.info$CellType=rep(NA,length(dge@ident))
levels(dge@data.info$CellType)=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Unknown","Leydig","Sertoli")
dge@data.info$CellType[which(dge@data.info$res.0.5==10)] <- "InnateLymphoid"
dge@data.info$CellType[which(dge@data.info$res.0.5==9)] <- "Macrophage"
dge@data.info$CellType[which(dge@data.info$res.0.5==8)] <- "Endothelial"
dge@data.info$CellType[which(dge@data.info$res.0.5==11)] <- "Myoid"
dge@data.info$CellType[which(dge@data.info$res.0.5 %in% c(0,4,5))] <- "Unknown"
dge@data.info$CellType[which(dge@data.info$res.0.5==6)] <- "Leydig"
dge@data.info$CellType[which(dge@data.info$res.0.5 %in% c(1,2,3,7))] <- "Sertoli"
table(dge@data.info$CellType)
#   Endothelial InnateLymphoid         Leydig     Macrophage          Myoid 
#           179             64            314            139             49 
#       Sertoli        Unknown 
#          2131           2205
dge <- SetAllIdent(dge, id = "CellType")
levels(dge@ident)=levels(dge@data.info$CellType)
dgeSomaticall=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown")

### color scheme for 7 somatic cell types
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7)] 

### Visualize 7 somatic cell types in PCA for 25 ST somatic subset - Figure 1B right panel and Figure 5A
pdf(paste0(dgefile,"PCA_7CellTypes.pdf"),width=16,height=10)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)    # Figure 1B right panel
plot6=PCAPlot(dge,1,6,do.return = TRUE,pt.size = 1,do.label=F)
plot7=PCAPlot(dge,1,7,do.return = TRUE,pt.size = 1,do.label=F)
plot8=PCAPlot(dge,1,8,do.return = TRUE,pt.size = 1,do.label=F)
plot9=PCAPlot(dge,1,9,do.return = TRUE,pt.size = 1,do.label=F)
plot10=PCAPlot(dge,1,10,do.return = TRUE,pt.size = 1,do.label=F)  # Figure 5A
plott=TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5,plot6,plot7,plot8,plot9,plott),cols = 3)
dev.off()
## save as Figure 1B right panel and Figure 5A

### Visualize somatic cells of individual batch in PCA
orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
sets=levels(orig.ident)
plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
plotset[[i]]=PCAPlot(dge,1,5,do.return = TRUE,do.label=T,pt.size = 1,cells.use=rownames(dge@data.info)[which(gsub("_.*","",rownames(dge@data.info)) == set)])
}
pdf(paste(dgefile,"dge_PCA_SetPC5_samedotsize.pdf",sep=""),height=12,width=16)
par(mar=c(1,1,3,1))
MultiPlotList(plotset[1:12],cols = 4)
MultiPlotList(plotset[13:24],cols = 4)
MultiPlotList(plotset[c(25,1:11)],cols = 4)
dev.off()

### INT6 Vs. somatic cells of 24 ST datasets in PCA
sets=list("INT6",levels(orig.ident)[-17])
plotset=NULL
for(i in 1:2){
set=sets[[i]]
plotset[[i]]=PCAPlot(dge,1,5,do.return = TRUE,do.label=T,pt.size = 1,cells.use=rownames(dge@data.info)[which(gsub("_.*","",rownames(dge@data.info)) %in% set)])
}
pdf(paste(dgefile,"dge_PCA_INT6PC5_samedotsize.pdf",sep=""),height=4,width=8)
par(mar=c(1,1,3,1))
MultiPlotList(plotset[1:2],cols = 2)
dev.off()

### INT6 Vs. somatic cells of 24 ST datasets in tSNE
sets=list("INT6",levels(orig.ident)[-17])
plotset=NULL
for(i in 1:2){
set=sets[[i]]
plotset[[i]]=TSNEPlot(dge,do.return = TRUE,do.label=T,pt.size = 1,cells.use=rownames(dge@data.info)[which(gsub("_.*","",rownames(dge@data.info)) %in% set)])
}
pdf(paste(dgefile,"dge_tSNE_INT6_samedotsize.pdf",sep=""),height=4,width=8)
par(mar=c(1,1,3,1))
MultiPlotList(plotset,cols = 2)
dev.off()

### identify Differentially-Expressed Markers
library(FNN)
library(igraph)
markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)
table(markersall$cluster)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST25Somatic_MarkersAll_7celltypes_pct0.2_diffpct0.2_thresh2fold_1.8.2018.txt"),col.names=T,row.names=T,quote=F,sep="\t")


######## Visualize representative markers for somatic cell types in tSNE view - Figure 5B
### representative somatic markers for each somatic cell type
macro=c("Adgre1","Dab2","Apoe","Mrc1")
endo=c("Tie1","Vwf","Tek")
myoid=c("Acta2","Myh11","Acta1","Pcna")
unknown=c("Pdgfra","Ly6a","Thra","Thrb","Tcf21","Arx")
ley=c("Hsd3b1","Cyp17a1","Cyp11a1","Star","Insl3")
ser=c("Sox9","Clu","Vim")
innate=c("Id2","Il7r","Rora","Thy1","Ccl5","Cd52","Cd2")
allgene=c("Macrophage","Endothelia","Myoid","Unknown","Leydig","Sertoli","InnateLymphoid")
allgenelist=list(macro,endo,myoid,unknown,ley,ser,innate)
### color scheme
redblue100.alpha<-rgb(read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)
###
object=object2=dge
for(i in 1:length(allgene)){
setname=allgene[i]
genes=allgenelist[[i]]
for(j in 1:length(genes)){
feature=features.plot=genes[j]
dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
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
data.gene0 <- na.omit(data.frame(data.use[feature, ]))
data.plot$gene <- t(data.gene0)
st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]
### plot heatmap for expression of representative markers for each of the 7 somatic cell types
jpeg(paste0(dgename,setname,"_",feature,"_40-100redblue0.8.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[40:100][(z - min(z))/diff(range(z))*60 + 1]
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
jpeg(paste0(dgename,setname,"_",feature,"_40-100top1percent.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
zcolor <- redblue100.alpha[40:100][c((z - min(z))/diff(range(z))*60 + 1, rep(60,top1n))]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()
}
}
### save as Figure 5B


###### Heatmap for Rank correlation of all somatic cells (N=5081) 
### calculate rank correlation for all filtered somatic cells (N=5081) using HVG
nrho=cor(as.matrix(dge@data)[dge@var.genes,],method="spearman")
### order cells by 7 somatic cell types and randomly shuffling cells within each cell type
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
### rank correlation for somatic cells ordered by 7 somatic cell types
            data.use2=testcor[cells.use,cells.use]
            data.use2=minmax(data.use2,min=disp.min,max=disp.max)
### labeling for axis 
            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+15]=levels(cells.ident)
            row.lab2=gsub(".*_","",lab2)
            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+15]=levels(orig.ident)
            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) # draw a line between datasets
### color bar for axis
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
redblue100<-rgb(read.table('/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2) 
### heatmap of rank correlation for somatic cells in 7 somatic cell types
jpeg(file=paste(dgename,"RankedCorrelation_HVG_res.0.5.jpeg",sep=""),height=3500,width=3500,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=2,cexRow=2,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(8,6))                    # symm=F,symkey=F,symbreaks=F,
dev.off()

###### Heatmap of Jaccard distance for somatic cells in 7 somatic cell types - Figure S1E
### Jaccard distance for somatic cells in 7 somatic cell types
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
### heatmap of Jaccard Distance for somatic cells in 7 somatic cell types - Figure S1E
jpeg(file=paste(dgename,"SNN.jpeg",sep=""),height=3500,width=3500,res=300)
par(mar=c(10,4,1,2),mgp=c(2.5, 1, 0))
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black",sepwidth=c(0.01,0.01),RowSideColors=rlab2,ColSideColors=clab2,labCol=col.lab2,labRow=row.lab2,cexCol=2,cexRow=2,ColSideColorsSize = 1.5,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F,scale="none",margins=c(9,7))                    # symm=F,symkey=F,symbreaks=F,
dev.off()
### save as Figure S1E
