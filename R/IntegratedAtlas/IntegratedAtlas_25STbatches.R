### R script for analyzing merged 25 ST batches in Dec 2017 by Qianyi
### Related to Figure 1B-C, and Table S2

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

######### 25 Adult Mouse ST batches (24 batches + INT6)
### path
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"

### 1. Load filtered cells and all detected genes for old 24 ST batches
load(file =paste0(home,"data_DGE/MouseAdultST24genesall.Robj"))
dgeall # 36979 genes across 33180 samples

old24=dgeall@raw.data
old24=as.matrix(old24)
print(dim(old24)) # [1] 36979 33180
old24=data.frame(GENE=rownames(old24),old24)
old24[1:5,1:2]
                       GENE ST1_ATCTTGCACATC
0610005C13Rik 0610005C13Rik                0
0610007P14Rik 0610007P14Rik                0
0610009B22Rik 0610009B22Rik                6
0610009E02Rik 0610009E02Rik                0
0610009L18Rik 0610009L18Rik               11

### 2. load filtered cells and all detected genes for INT6 (somatic cells of new Sca1) dataset
newSca1somatic=read.table(file=paste0(home,"data_DGE/92783/INT6Somatic_1453cells_25187genes.dge.txt"), row.names=1,header=T)
newSca1somatic=data.frame(GENE=rownames(newSca1somatic),newSca1somatic)
dim(newSca1somatic) # [1] 25187  1453+1
newSca1somatic[1:5,1:2]
                       GENE INT6_GATTAGGCCATA
0610005C13Rik 0610005C13Rik                 1
0610007P14Rik 0610007P14Rik                 0
0610009B22Rik 0610009B22Rik                 0
0610009E02Rik 0610009E02Rik                 0
0610009L18Rik 0610009L18Rik                 1

### 3. Merge 24 ST datasets with INT6 dataset 
mouse=merge(old24,newSca1somatic,by="GENE",all=TRUE)
dim(mouse) # 37241 genes, 34633 cells + 1 column of gene name
mouse[is.na(mouse)] <- 0  # change NA to 0
row.names(mouse)=mouse[,1]
mouse[1:5,1:2]
mouse=mouse[,-1]
dim(mouse) # 37241 genes, 34633 cells
mouse[1:5,1:2]
              ST1_ATCTTGCACATC ST1_CCGCCTTAGCTN
0610005C13Rik                0                0
0610007P14Rik                0                0
0610009B22Rik                6                3
0610009E02Rik                0                0
0610009L18Rik               11                6

write.table(mouse, file=paste0(home,"data_DGE/MouseAdultST25/mergedMouseAdultST25_34633cells_37241genes.dge.txt"), quote=FALSE, sep='\t',row.names=T, col.names=T)

### Save as sparse matrix
library(Matrix)
sparse1=as(as.matrix(mouse),"sparseMatrix")   
writeMM(sparse1,paste0(home,"data_DGE/MouseAdultST25/matrix.mtx"))
write.table(colnames(sparse1), file=paste0(home,"data_DGE/MouseAdultST25/barcodes.tsv"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(cbind(rownames(sparse1),rownames(sparse1)), file=paste0(home,"data_DGE/MouseAdultST25/genes.tsv"), quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE)

### clear workspace
rm(list=ls())      # remove all objects from current workspace
gc()               # clear memory of removed objects


###### load merged gene expression matrix for 25 ST batches
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
mouse <- Read10X(paste0(home,"data_DGE/MouseAdultST25"))
dim(mouse) # [1] 37241 34633
mouse[1:5,1:2]
              ST1_ATCTTGCACATC ST1_CCGCCTTAGCTN
0610005C13Rik                0                0
0610007P14Rik                0                0
0610009B22Rik                6                3
0610009E02Rik                0                0
0610009L18Rik               11                6

dgedata=mouse
dgefile=dgename="figDec2017_MouseAdultST25/All_"
dataset=unique(gsub("_.*","",colnames(dgedata)))[c(1:16,25,17:24)]
length(dataset) # 25=6(ST)+2(1nDepleted)+3(SPG)+6(INT)+8(SER)

### 4. Filter genes for merged 25 ST datasets
nUMIperGene <- rowSums(dgedata)
nCellperGene <- rowSums(dgedata>0)

dim(dgedata)[2]*0.0012 # 41.
dim(dgedata)[2]*0.0006 # 20.

thres.nUMIperGene = 20
thres.nCellperGene = 15
length(nUMIperGene)
length(which(nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene))     #19356
# [1] 24923

gene1=rownames(dgedata)[nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene]
length(gene1)   # 24923

### add back known genes of interest with <=20 UMIs/15cells per gene
# 31 known genes of interest with low UMI/cells per gene
gene2=c("Ngf","T","Tbx20","Tert","Olig2","Prl","Hmx2","Itgb6","Mageb1","Neurod1","Neurod2","Ptchd4","Prl6a1","Fgf17","Hpx","Npy4r","Fgf6","Foxf2","Prokr2","Sry","Ccl3","Prf1","Pax8","Dcx","Ly6g5c","Ereg","Pgr","Sycp2l","Prnd","Cysltr1","Pdf")
length(gene2) # 31
length(which(gene2 %in% gene1)) # 7
gene2=gene2[which(!(gene2 %in% gene1))]
length(gene2) # 24
# 24 known genes of interest with <=20 UMIs/15cells per gene

### kept genes with >20 UMIs and expressed in >Cells, plus 24 genes of interest with UMI<=20
genekeep=c(gene1,gene2)
anyDuplicated(genekeep) # 0
length(genekeep) # 24947

### save the filtered cell and filtered genes
dgedata2=dgedata[genekeep,]
dim(dgedata2) # [1] 24947 34633


###### 5. Normalization for cells, log-transformation and standardization for genes
### Filtered cells, all detected genes
dge <- new("seurat", raw.data = dgedata)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Dec2017-MouseAdultST25all",names.field = 1,names.delim = "_")
dge
dgeall=dge
# 37241 genes across 34633 samples.
save(dgeall, file =paste0(home,"data_DGE/MouseAdultST25genesall.Robj"))

### Filtered cells, filtered genes
dge <- new("seurat", raw.data = dgedata2)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Dec2017-MouseAdultST25",names.field = 1,names.delim = "_")
dge
# 24947 genes across 34633 samples.
dge25=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))

which(table(gsub("_.*","",colnames(dgedata2))) != table(gsub("_.*","",colnames(dge@data)))) # named integer(0)

# change order of dataset in orig.ident
dataset=unique(gsub("_.*","",colnames(dge@data)))[c(1:16,25,17:24)]
dge@data.info$orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
dge@ident=factor(dge@ident,levels=dataset)

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes) # 32 dgeall, 27 dge25
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")


###### add major cell types from 24 ST batches to metadata of 25 ST batches
### note: merged all somatic cells as "Somatic" group
### 1 Somatic group (somatic group from 24 ST + INT6) and 4 Major germ cell types
load(file = paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge20=dge

celltype2=as.character(dge20@data.info$celltype2)
celltype1=rep("Somatic",length(which(dge25@data.info$orig.ident=="INT6")))
names(celltype2)=rownames(dge20@data.info)
names(celltype1)=rownames(dge25@data.info)[which(dge25@data.info$orig.ident=="INT6")]
length(celltype2) # 33180: somatic group and 4 major germ cell types from 24 ST
length(celltype1) # 1453: somatic group from INT6
oldcelltype=c(celltype2,celltype1)
length(oldcelltype)  # 34633
table(gsub("_.*","",names(oldcelltype)))
#INT1 INT2 INT3 INT4 INT5 INT6 SER1 SER2 SER3 SER4 SER5 SER6 SER7 SER8 SPG1 SPG2
#819  237 1210  653 1480 1453   81  161  451  716 2015 1399 3047 1979  957  937
#SPG3  ST1  ST2  ST3  ST4  ST5  ST6  ST7  ST8
#1212 2324 1642 2720 2363 1601 2341 1549 1286

celltype2=oldcelltype[rownames(dge25@data.info)]
celltype2=factor(oldcelltype,levels=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))
table(celltype2)
#       Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
#          5081           4239           8792           9923           6598 
sum(table(celltype2)) # 34633

### add major cell groups to metadata
dge=dge25
dge=AddMetaData(dge,celltype2,"celltype2")
dge@data.info$celltype2=factor(dge@data.info$celltype2,levels=levels(celltype2))
dge=SetAllIdent(dge,id="celltype2")
dge@ident=celltype2
table(dge@ident)
#       Somatic  Spermatogonia   Spermatocyte RoundSpermatid     Elongating 
#          5081           4239           8792           9923           6598 
dge25=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))



###### 6. PCA for merged 25 datasets 
dge=dge25
Sys.time() # [1] "2017-12-04 12:31:08 EST"
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
Sys.time() # [1] "2017-12-05 21:12:24 EST"
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))
dge25=dge

### PCA plot for merged 25 ST datasets (plot 1 somatic arm as grey)  - Figure 1B left panel
load(file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))
dge25=dge
dge <- SetAllIdent(dge, id = "celltype2")
dge@ident=factor(dge@ident,levels=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating"))
# flip PCs
dge@pca.obj[[1]]$x[,1]=-dge@pca.obj[[1]]$x[,1]
dge@pca.obj[[1]]$x[,2]=-dge@pca.obj[[1]]$x[,2]
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.rot[,2]=-dge@pca.rot[,2]
dge25=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))
# color scheme
library(RColorBrewer)
myBrewerPalette <- c("gray60",brewer.pal(12,"Paired")[1:4]) # used this on 12/20/2017
# PCA plot for merged 25 ST datasets with 1 somatic arm and 4 major germ cell types - Figure 1B left panel
pdf(paste0(dgefile,"PCA_5CellTypeswithMergedSomaticArm.pdf"),width=15,height=12)
plot2=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=F)
plot3=PCAPlot(dge,1,3,do.return = TRUE,pt.size = 1,do.label=F)
plot4=PCAPlot(dge,1,4,do.return = TRUE,pt.size = 1,do.label=F)
plot5=PCAPlot(dge,1,5,do.return = TRUE,pt.size = 1,do.label=F)
MultiPlotList(list(plot2,plot3,plot4,plot5),cols = 2)
dev.off()
## save as Figure 1B left panel

### Scree Plot for single PCA
numPCs=55 # 24947: genes with >20UMIs and >15cells plus 24 genes of interest
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

###### Louvain-Jaccard clustering of merged 25 datasets 
### store top 60 PC scores
# object@pca.rot only stored top 40 PCs
dim(dge@pca.rot)		# 34633    40
dim(dge@pca.obj[[1]]$x) # 34633 24947
dge@pca.rot[1:5,1:3]
dge@pca.obj[[1]]$x[1:5,1:3]
dge@pca.rot=as.data.frame(dge@pca.obj[[1]]$x[,1:60])

### Louvain-Jaccard clustering using top 55 PCs
Sys.time() # [1] "2017-12-18 09:07:27 EST"
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.1,2,0.1), print.output = 0, save.SNN = T)
Sys.time() # [1] "2017-12-18 10:37:14 EST"
Sys.time() # [1] "2017-12-18 11:53:12 EST"
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(2.1,4,0.1), print.output = 0, save.SNN = T)
Sys.time() # [1] "2017-12-18 13:09:06 EST"
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15clusters.Robj"))
dge25=dge

### tSNE using top 55 PCs
Sys.time() # [1] "2017-12-18 11:53:12 EST"
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)    # max_iter=2000
Sys.time() # [1] "2017-12-18 13:09:06 EST"

save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))
dge25=dge

library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7,1:4)] # 12/20/2017
pdf(paste0(dgefile,"tSNE_11CellTypes.pdf"),width=7,height=5)
TSNEPlot(dge,1,2,do.return = TRUE,do.label=F)
dev.off()
write.table(dge@tsne.rot,paste0(home,"data_DGE/MouseAdultST25_tSNE.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### number of clusters for 25 ST 
print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)),length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.6)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2)) ))

### cross tabulation with old cell types from merged 24 ST datasets and from INT6 somatic cells
table(dge@data.info$celltype2,dge@data.info$res.0.1)
# cluster 5,6,8 and 9 cells are the same cells as in somatic group, cluster 0-4 and 7 are the same as cells as in germ group. 

### decided to retain 1 somatic group and 4 major germ cell types in celltype2
### do focused clustering for somatic subset of 25 ST datasets (N=5081 cells) 

###### update Somatic Cell Types (New Somatic Cell Types from focused subset clustering of merged 25ST somatic cells) to metadata of merged 25 ST batches
### note: 11 major cell types = 4 major germ cell groups from 24 ST batches + 7 somatic cell types from re-clustering of all somatic cells of 25 ST batches
load(file = paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge20=dge

load(file =paste0(home,"data_DGE/MouseAdultST25Somatic.Robj"))
dgeSomaticall=dge

celltype2=as.character(dge20@data.info$celltype2)
newsomaticcelltype=as.character(dgeSomaticall@data.info$CellType)
names(celltype2)=rownames(dge20@data.info)
names(newsomaticcelltype)=rownames(dgeSomaticall@data.info)
germcelltype=celltype2[which(!(names(celltype2) %in% names(newsomaticcelltype)))]
length(germcelltype)       # 29552: 4 major germ cell groups from 24 ST batches
length(newsomaticcelltype) # 5081:  7 somatic cell types from re-clustering of all somatic cells of 25 ST batches
celltype=c(germcelltype,newsomaticcelltype)
length(celltype)           #  34633
table(gsub("_.*","",names(celltype)))
#INT1 INT2 INT3 INT4 INT5 INT6 SER1 SER2 SER3 SER4 SER5 SER6 SER7 SER8 SPG1 SPG2
#819  237 1210  653 1480 1453   81  161  451  716 2015 1399 3047 1979  957  937
#SPG3  ST1  ST2  ST3  ST4  ST5  ST6  ST7  ST8
#1212 2324 1642 2720 2363 1601 2341 1549 1286
sum(table(celltype)) # 34633
celltype=celltype[rownames(dge25@data.info)]
write.table(celltype,file =paste0(home,"data_DGE/11celltypes_MouseAdultST25.txt"),quote=F,row.names=T,col.names=F,sep="\t")

### load 11 major cell types
celltypes=t(read.table(paste0(home,"data_DGE/11celltypes_MouseAdultST25.txt"),row.names=1))
celltype=celltypes[1,]
names(celltype)=colnames(celltypes)
length(which(is.na(celltype)))	# 0
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
celltype=factor(celltype,levels=levels)
sum(table(celltype)) 			# 34633
table(celltype)
#InnateLymphoid     Macrophage    Endothelial          Myoid         Leydig 
#            64            139            179             49            314 
#       Sertoli        Unknown  Spermatogonia   Spermatocyte RoundSpermatid 
#          2131           2205           4239           8792           9923 
#    Elongating 
#          6598 

### add 11 major cell type to metadata
dgeall=AddMetaData(dgeall,celltype,"celltype")
dgeall@data.info$celltype=factor(dgeall@data.info$celltype,levels=levels)
dgeall=SetAllIdent(dge,id="celltype")
dgeall@ident=factor(dgeall@ident,levels=levels(dgeall@data.info$celltype))
save(dgeall, file =paste0(home,"data_DGE/MouseAdultST25genesall.Robj"))

dge=dge25
dge=AddMetaData(dge,celltype,"celltype")
dge@data.info$celltype=factor(dge@data.info$celltype,levels=levels)
dge=SetAllIdent(dge,id="celltype")
dge@ident=factor(dge@ident,levels=levels)
dge25=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST25genesUMI20cell15.Robj"))


######### Identify Differentially-expressed Markers for 11 major cell types - Table S2 and Figure 1C
load(file =paste0(home,"data_DGE/MouseAdultST25genesall.Robj"))
dgefile=dgename="figDec2017_MouseAdultST25/All_"
dge=dgeall
### Identify differentially-expressed markers using binomial likelihood test - Table S2
library(FNN)
library(igraph)
dge=dgeall
dge=SetAllIdent(dge,id="celltype")
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
dge@ident=factor(dge@ident,levels=levels)
Sys.time()  # [1] "2017-12-19 19:26:53 EST 
markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST25_MarkersAll_",genefilter[1],"_merged11celltypes_pct0.2_diffpct0.2_thresh2fold_clustersPC1-40_12.17.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")
Sys.time()  # [1] "2017-12-19 21:39:30 EST"

### Add average expression level of each markers for each cluster for 1) all cells in the cluster; and 2) the cells positively expressing the specific marker in the cluster
table(markersallin$cluster[!grepl("V",markersallin$cluster)])
markersall=markersallin[!grepl("V",markersallin$cluster),]
markersall$cluster=as.numeric(markersall$cluster)
meanmarker=matrix(0,nrow=dim(markersall)[1],ncol=3)
meanmarker[,1]=markersall$gene
for(markerj in 1:dim(markersall)[1]){
marker=markersall[markerj,6]
cluster=as.character(markersall[markerj,5])
tmp=dge@data[marker,which(dge@ident==cluster)]
meanmarker[markerj,2]=expMean(tmp)
meanmarker[markerj,3]=expMean(tmp[tmp!=0])
print(markerj)
}
markersallmean=cbind(markersall,meanmarker)
which(markersallmean[,6] != markersallmean[,7])   # integer(0)
write.table(markersallmean[,c(6,1:5,8,9)],paste0(home,"data_DGE/MouseAdultST25_MarkersAllMean_UMI20cell15_11celltypes_pct0.2_diffpct0.2_thresh2fold_12.20.2017.txt"),col.names=T,row.names=F,quote=F,sep="\t")
### save as Table S2

###### Heatmap of all markers for 11 major cell types - Figure 1C
### Calcualte gene expression for each cluster centroid
expMean <- function(x) {
  return(log(x = mean(x = exp(x = x) - 1) + 1))
}
centroid=matrix(,dim(dge@data)[1],length(unique(dge@ident)))
rownames(centroid)=rownames(dge@data)
colnames(centroid)=unique(dge@ident)
for(i in levels(dge@ident)){
    centroid[,i]=apply(dge@data[,which(dge@ident==i)],1,function(x) expMean(as.numeric(x))) 
    print(i)
}
centroid2=AverageExpression(dge)
which(round(centroid,2) != round(centroid2,2)) # integer(0)

### Standardize genes Across Cell Types
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,"11celltypes_centroid_allgenes_std.txt",row.names=T,col.names=T,quote=F,sep="\t")
centroid.std=read.table(paste0(home,"data_DGE/11celltypes_centroid_allgenes_std.txt"),header=T,row.names=1)

### separate each cell type by a vertical line
levels=levels(dge@ident)
colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
col.lab=rep("",length(levels))
col.lab=gsub(".*_","",levels)
### color scheme for 11 major cell types
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7,1:4)] # 12/20/2017
### add color bars for markers of each cell type 
ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")
### color scheme for heatmap
redblue100<-rgb(read.table(paste0(home,"data_DGE/redblue100.txt"),sep='\t',row.names=1,header=T))
col.use=redblue100
### plot heatmap of all markers for 11 major cell types - Figure 1C
data.use=centroid.std
row.lab=rownames(data.use)
jpeg(file=paste0(dgename,"centroid_std_allgenes_cluster.jpeg"),res=300,height=6000,width=1600)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(na.omit(data.use),dendrogram="row",Rowv=TRUE,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=0.8,cexRow=0.3,ColSideColorsSize = 2,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
library(pheatmap)
jpeg(file=paste0(dgefile,"clusterforgenes.jpeg"),res=300,height=10000,width=1600)
pheatmap(testcor,cluster_cols=FALSE,col=col.use2)
dev.off()
### save as Figure 1C left panel




