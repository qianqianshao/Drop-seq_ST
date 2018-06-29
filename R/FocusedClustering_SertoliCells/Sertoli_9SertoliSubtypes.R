### R script to generate 9 sertoli subtypes by Qianyi 
### Related to Figure 6B,D, Figure S6A, and Table S6: sertoli cells (N=20,646) with >1k detected genes 

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)

### load paths
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
dgefile=dgename="/home/qzm/figNov2017_MouseAdultST24_ReCluster45largecells/"
redblue100<-rgb(read.table('/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt',sep='\t',row.names=1,header=T))

### two sets of analysis: 1) all filtered sertoli cells, and 2) sertoli cells with >1k detected genes
sets=list(c(4:5),c(4:5))
setsname=c("Sertoli","SertoliGenes1k")
setslabel=c("Sertoli","SertoliGenes1k")
dgelist=list()

### load object of 24 ST datasets (INT6 did not contribute to sertoli cells)
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))

###### zoomed-in view of PCA using PCA of 24 ST datasets
j=1
cells.use=names(dge@ident[which(dge@ident %in% sets[[j]])])
pdf(paste(dgefile,"ZoomedInView_PCA_SertoliSubset.pdf",sep=""),height=5,width=6)
PCAPlot(dge,1,2,cells.use=cells.use,do.return = TRUE,pt.size = 1,do.label=T)
PCAPlot(dge,1,3,cells.use=cells.use,do.return = TRUE,pt.size = 1,do.label=T)
dev.off()

######### Focused clustering for Sertoli cells subset (N=2,099) on 6/28/2017
dgefile="figJun2017_MouseAdultST24_ReCluster45/"
###### PCA for sertoli cells subset
### extract cells in Sertoli clusters
dgetmp=SubsetData(dge,ident.use=sets[[j]])
### extract detected genes in sertoli cells
nCellperGene <- rowSums(as.matrix(dgetmp@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgetmp@data=dgetmp@data[genes.use,]
print(table(dgetmp@ident))
print(dim(dgetmp@data))
print(Sys.time())
dgetmp <- PCA(dgetmp, pc.genes = rownames(dgetmp@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dgetmp,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster_",setsname[j],".pdf",sep=""),height=5,width=6)
PCAPlot(dgetmp,1,2,do.return = TRUE,pt.size = 1,do.label=T)
PCAPlot(dgetmp,1,3,do.return = TRUE,pt.size = 1,do.label=T)
dev.off()
dgelist[[j]]=dgetmp

### Scree Plot for single PCA
numPCs=14
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:150],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,5),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.8,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()

###### Louvain-Jaccard clustering and tSNE using top PCs
dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,3,0.1), print.output = 0, save.SNN = T)

dgelist[[j]]=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))

### check number of clusters
print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)),length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.6)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2)),length(unique(dge@data.info$res.2.1)),length(unique(dge@data.info$res.2.2)),length(unique(dge@data.info$res.2.3)),length(unique(dge@data.info$res.2.4)),length(unique(dge@data.info$res.2.5)),length(unique(dge@data.info$res.2.6)),length(unique(dge@data.info$res.2.7)),length(unique(dge@data.info$res.2.8)),length(unique(dge@data.info$res.2.9)),length(unique(dge@data.info$res.3)) ))
res="res.1";i=1
ncluster=9

### Build Cluster tree 
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])
dge@data.info[,res[i]]=dge@ident

###### order the 9 sertoli cell clusters using optimal leaf reordering by Seriation
### for each cluster, calculate average normalized expression of each gene
centroid=matrix(,dim(dge@data)[1],length(unique(dge@ident)))
rownames(centroid)=rownames(dge@data)
colnames(centroid)=unique(dge@ident)
for(i in levels(dge@ident)){
    centroid[,i]=apply(dge@data[,which(dge@ident==i)],1,function(x) expMean(as.numeric(x))) 
    print(i)
}

### Reordering cluster centroid using Seriation
library(seriation)
n=ncluster=length(levels)
pdf(file=paste0(dgename,"Centroid_norm_Seriation_Dissimilarity_",res[resi],".pdf"))
tmp=centroid[,levels][,1:(sum(ncluster))]
colnames(tmp)=gsub(".*_","",colnames(tmp))
da <- dist(t(as.matrix(tmp)), method = "euclidean")
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

ordered="res.1order"

dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@data.info[,ordered]=factor(dge@data.info[,ordered])
dge=SetAllIdent(dge,ordered)
dgelist[[i]]=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))

###### visualize in PCA & tSNE
pdf(paste0(dgefile,"cluster_res.1.pdf"),width=5.5,height=5)
print(table(dge@data.info$clusters31_GeneFilter1,dge@ident))
PCAPlot(dge,pt.size = 1,do.label=TRUE)
TSNEPlot(dge,pt.size = 1,do.label=TRUE,label.size=8)
dgelist[[i]]=dge

### Contribution of each batch to each cluster
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefile,"Cluster",setsname[j],res[i],"_ncellspercluster_batch.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### Distribution of Number of UMIs, %MT in Each Cluster
  nUMI=MT=NULL
## nUMI per Cell
for(id in 1:ncluster){
  nUMI=rbind(nUMI,summary(dge@data.info$nUMIperCell2[which(dge@ident==id)]))
}
## %MT per cell
for(id in 1:ncluster){
  MT=rbind(MT,summary(dge@data.info$percent.mito[which(dge@ident==id)]))
}
print(nUMI)
print(MT)

pdf(file=paste(dgefile,"NGeneUMImt_Sertolicells.pdf",sep=""),height=4,width=11)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
dev.off()

 
######### Focused clustering for Sertoli cells with >1k genes (N=1,067) on 11/28/2017
dgefile=dgename="/home/qzm/figNov2017_MouseAdultST24_ReCluster45largecells/"
### remove small cells in Sertoli subset: remove all cells <= 1000 genes
j=2
  dge=dgelist[[1]]
  cells.use=rownames(dge@data.info)[which(dge@data.info$nGene>1000)]
  dgetmp=SubsetData(dge,cells.use=cells.use)
  dgetmp
  dgetmp=SetAllIdent(dgetmp,id="clusters31_GeneFilter1")
  print(table(dgetmp@ident))
  ncellsbatchcluster = table(dgetmp@data.info$orig.ident,dgetmp@ident)
  print(ncellsbatchcluster)
  print(prop.table(ncellsbatchcluster,2))
  print(prop.table(ncellsbatchcluster,1))
  dgetmp=SetAllIdent(dgetmp,id="res.1order")
  print(table(dgetmp@ident))
  ncellsbatchcluster = table(dgetmp@data.info$orig.ident,dgetmp@ident)
  print(ncellsbatchcluster)
  print(prop.table(ncellsbatchcluster,2))
  print(prop.table(ncellsbatchcluster,1))
### extract genes detected in sertoli cells with >1k genes
nCellperGene <- rowSums(as.matrix(dgetmp@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgetmp@data=dgetmp@data[genes.use,]
print(table(dgetmp@ident))
print(dim(dgetmp@data))
print(Sys.time())
dgetmp <- PCA(dgetmp, pc.genes = rownames(dgetmp@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dgetmp,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))
dgelist[[j]]=dgetmp

### plot PCA
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster_",setsname[j],".pdf",sep=""),height=5,width=6)
PCAPlot(dgetmp,1,2,do.return = TRUE,pt.size = 1,do.label=T)
PCAPlot(dgetmp,1,3,do.return = TRUE,pt.size = 1,do.label=T)
dev.off()

### Scree Plot for single PCA
numPCs=11
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:150],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(density(dge@pca.obj[[1]]$sdev),xlim=c(0,5),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+0.8,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()

###### Louvain-Jaccard clustering and tSNE using top PCs
dge=RunTSNE(dge,dims.use = 1:numPCs[j],do.fast=T)    # max_iter=2000
dge <- FindClusters(dge, pc.use = 1:numPCs[j], resolution = seq(0.1,3,0.1), print.output = 0, save.SNN = T)

dgelist[[j]]=dge
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))

### check number of clusters
print(c( length(unique(dge@data.info$res.0.1)),length(unique(dge@data.info$res.0.2)),length(unique(dge@data.info$res.0.3)),length(unique(dge@data.info$res.0.4)),length(unique(dge@data.info$res.0.5)),length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)),length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.6)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2)),length(unique(dge@data.info$res.2.1)),length(unique(dge@data.info$res.2.2)),length(unique(dge@data.info$res.2.3)),length(unique(dge@data.info$res.2.4)),length(unique(dge@data.info$res.2.5)),length(unique(dge@data.info$res.2.6)),length(unique(dge@data.info$res.2.7)),length(unique(dge@data.info$res.2.8)),length(unique(dge@data.info$res.2.9)),length(unique(dge@data.info$res.3)) ))

### Build Cluster tree 
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])
dge@data.info[,res[i]]=dge@ident

### use res.1.2 for 9 sertoli subtypes
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster45cellsGenes1k.Robj"))
res="res.1.2";resi=1
ncluster=9
dge=SetAllIdent(dge,id="res.1.2")
table(dge@ident)

### cross-tabulation with 9 clusters from clustering of all sertoli cells
table(dge@ident,dge@data,info$res.1order)

### modify the cluster IDs for 9 sertoli subtypes
id=as.numeric(dge@ident)
names(id)=names(dge@ident)
id[which(id==2)]<-"2A"
id[which(id==3)]<-"2B"
id[which(id==4)]<-"3A"
id[which(id==5)]<-"3B"
id[which(id==6)]<-"4A"
id[which(id==7)]<-"4C"
id[which(id==8)]<-"4D"
id[which(id==9)]<-"4B"
table(id)
table(dge@ident)
id=factor(id,levels=c("4D","4C","3B","1","2B","4B","4A","3A","2A"),order=T)
dge=AddMetaData(dge,id,"Cluster")
dge=SetAllIdent(dge,id="Cluster")
dge@ident=id
dgelist[[2]]=dge
save(dge,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster45cellsGenes1k.Robj"))


###### visualize in PCA & tSNE - Figure 6B
pdf(paste0(dgefile,"cluster_9subtypes.pdf"),width=5.5,height=5)
print(table(dge@data.info$clusters31_GeneFilter1,dge@ident))
PCAPlot(dge,pt.size = 1,do.label=TRUE)
TSNEPlot(dge,pt.size = 1,do.label=TRUE,label.size=8)  # Figure 6B
dev.off()
# save as Figure 6B

### Contribution of each batch to each cluster
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefile,"Cluster",setsname[j],res[i],"_ncellspercluster_batch.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### Distribution of Number of UMIs, %MT in Each Cluster - Figure S6A
  nUMI=MT=NULL
## nUMI per Cell
for(id in 1:ncluster){
  nUMI=rbind(nUMI,summary(dge@data.info$nUMIperCell2[which(dge@ident==id)]))
}
## %MT per cell
for(id in 1:ncluster){
  MT=rbind(MT,summary(dge@data.info$percent.mito[which(dge@ident==id)]))
}
print(nUMI)
print(MT)

library(RColorBrewer)
myBrewerPalette=brewer.pal(9,"Set1")[c(1:6,8,9,7)]
pdf(file=paste(dgefile,"NGeneUMImt.pdf",sep=""),height=4,width=11)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette,y.log=T)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
dev.off()
# save as Figure S6A

### % Cell Cycle genes in each cluster
G1S=c("Acd","Acyp1","Adamts1","Ankrd10","Apex2","Arglu1","Atad2","Bard1","Brd7","Capn7","Casp2","Casp8Ap2","Ccne1","Ccne2","Cdc25A","Cdc6","Cdca7","Cdca7L","Cep57","Chaf1A","Chaf1B","Clspn","Crebzf","Ctsd","Dis3","Dnajc3","Donson","Dscc1","Dtl","E2F1","Eif2A","Esd","Fam122A","Flad1","Gins2","Gins3","Gmnn","Gon7","Hells","Hoxb4","Hras","Hsf2","Insr","Ints8","Ivns1Abp","Lnpep","Luc7L3","Mcm2","Mcm4","Mcm5","Mcm6","Mdm1","Med31","Mri1","Msh2","Mturn","Nasp","Neat1","Nktr","Npat","Nup43","Orc1","Osbpl6","Otulin","Pank2","Pcdh7","Pcna","Plcxd1","Pms1","Pnn","Pold3","Rab23","Recql4","Rmi2","Rnf113a1","Rnpc3","Rsrp1","Sec62","Skp2","Slbp","Slc25A36","Snhg10","Srsf7","Ssr3","Taf15","Tipin","Topbp1","Tra2A","Ttc14","Ubr7","Uhrf1","Ung","Usp53","Vps72","Wdr76","Zfp367","Zmynd19","Zranb2")
S=c("Abcc5","Abhd10","Asf1B","Atad2","Bbs2","Bivm","Blm","Bmi1","Brca1","Brip1","Cald1","Calm2","Casp2","Ccdc14","Ccdc150","Ccdc84","Cdc45","Cdc7","Cdca5","Cdkn2Aip","Cenpm","Cenpq","Cenpu","Cers6","Chml","Coq9","Cpne8","Crebzf","Crls1","Ddias","Depdc7","Dhfr","Dna2","Dnajb4","Donson","Dscc1","Dync1Li2","E2F8","Eif4Ebp2","Esco2","Exo1","Ezh2","Fam178A","Fanca","Fanci","Fen1","Gclm","H1F0","Hells","Hist1H2Ac","Hist1H4C","Ints7","Kat2A","Kat2B","Kdelc1","Lmo4","Lyrm7","Man1A2","Map3K2","Mastl","Mbd4","Mcm8","Mycbp2","Nab1","Neat1","Nfe2L2","Nrd1","Nsun3","Nt5Dc1","Nup160","Ogt","Orc3","Osgin2","Phip","Phtf1","Phtf2","Pkmyt1","Pola1","Prim1","Ptar1","Rad18","Rad51","Rad51Ap1","Rbbp8","Reep1","Rfc2","Rhobtb3","Rmi1","Rpa2","Rrm1","Rrm2","Rsrc2","Sap30Bp","Shtn1","Slc38A2","Sp1","Srsf5","Svip","Top2A","Ttll7","Tyms","Ube2T","Ubl3","Usp1","Zbed5","Zwint")
G2M=c("4930427A07Rik","Anln","Ap3D1","Arhgap19","Arl4A","Armc1","Asxl1","Atl2","Aurkb","Bclaf1","Bora","Brd8","Bub3","C330027C09Rik","Casp3","Cbx5","Ccdc107","Ccna2","Ccnf","Cdc16","Cdc25C","Cdca2","Cdca3","Cdca8","Cdk1","Cdkn1B","Cdkn2C","Cdr2","Cenpl","Cep350","Cfd","Cflar","Chek2","Ckap2","Ckap2L","Cyth2","Dcaf7","Dhx8","Dnajb1","Entpd5","Espl1","Fadd","Fam83D","Fan1","Fancd2","G2E3","Gabpb1","Gas1","Gas2L3","H2Afx","Haus8","Hint3","Hipk2","Hjurp","Hmgb2","Hn1","Hp1Bp3","Hrsp12","Ifnar1","Iqgap3","Katna1","Kctd9","Kdm4A","Kif11","Kif20B","Kif22","Kif23","Kif5B","Kifc1","Klf6","Kpna2","Lbr","Lix1L","Lmnb1","Mad2L1","Malat1","Melk","Mgat2","Mid1","Mis18Bp1","Mnd1","Ncapd3","Ncaph","Ncoa5","Ndc80","Neil3","Nfic","Nipbl","Nmb","Nr3C1","Nucks1","Numa1","Nusap1","Pif1","Pknox1","Polq","Ppp1R2","Psmd11","Psrc1","Rangap1","Rccd1","Rdh11","Rnf141","Sap30","Ska3","Smc4","Stat1","Stil","Stk17B","Suclg2","Tfap2A","Timp1","Tmpo","Tnpo2","Top2A","Traip","Trim59","Trmt2A","Ttf2","Tuba1A","Tubb2A","Tubb4a","Tubb4B","Tubd1","Uaca","Ube2C","Vps25","Vta1","Wsb1","Znhit2")
M=c("Ahi1","Akirin2","Ankrd40","Anln","Anp32B","Anp32E","Arhgap19","Arl6Ip1","Asxl1","Atf7Ip","Aurka","Birc2","Birc5","Bub1","Cadm1","Ccdc88A","Ccdc90B","Ccna2","Ccnb2","Cdc20","Cdc25B","Cdc27","Cdc42Ep1","Cdca3","Cenpa","Cenpe","Cenpf","Cep55","Cflar","Cit","Ckap2","Ckap5","Cks1B","Cks2","Cnot10","Cntrob","Ctcf","Ctnna1","Ctnnd1","Depdc1A","Depdc1B","Diaph3","Dlgap5","Dnaja1","Dnajb1","Dr1","Dzip3","E2F5","Ect2","Fam64A","Foxm1","Fyn","G2E3","Gadd45A","Gas2L3","Got1","Grk6","Gtse1","Hcfc1","Hmg20B","Hmgb3","Hmmr","Hn1","Hp1Bp3","Hps4","Hs2St1","Hspa13","Hspa8","Inadl","Kif14","Kif20B","Kif2C","Kif5B","Klf9","Lbr","Lmna","Mcm4","Mdc1","Mis18Bp1","Mki67","Mllt4","Mzt1","Ncapd2","Ncoa5","Nek2","Nuf2","Nup35","Nup98","Nusap1","Odf2","Oraov1","Pbk","Pcf11","Plk1","Poc1A","Pom121","Ppp1R10","Prpsap1","Prr11","Psmg3","Ptp4A1","Ptpn9","Pwp1","Qrich1","Rad51C","Rangap1","Rbm8A","Rcan1","Rere","Rnf126","Rnf141","Rnps1","Rrp1","Sephs1","Setd8","Sfpq","Sgo2A","Shcbp1","Smarcb1","Smarcd1","Spag5","Sptbn1","Srf","Srsf3","Ss18","Suv420H1","Tacc3","Thrap3","Tle3","Tmem138","Tnpo1","Tomm34","Tpx2","Trip13","Tsg101","Tsn","Ttk","Tubb4B","Txndc9","Txnrd1","Ube2D3","Usp13","Usp16","Vangl1","Wibg","Wsb1","Ywhah","Zc3Hc1","Zfp207","Zfx","Zmym1")
MG1=c("Agfg1","Agpat3","Akap13","Amd1","Anp32E","Antxr1","Bag3","Btbd3","Cbx3","Cdc42","Cdk7","Cdkn3","Cep70","Cnih4","Ctr9","Cwc15","Dcp1A","Dctn6","Dexi","Dkc1","Dnajb6","Dsp","Dynll1","Eif4E","Elp3","Fam189B","Fam60A","Fopnl","Foxk2","Fxr1","G3Bp1","Gata2","Gnb1","Grpel1","Gspt1","Gtf3C4","Hif1A","Hmg20B","Hmgcr","Hsd17B11","Hspa8","Ilf2","Jmjd1C","Kdm5B","Kif5B","Kpnb1","Kras","Larp1","Larp7","Lrif1","Lyar","Morf4L2","Mrpl19","Mrps18B","Mrps2","Msl1","Mtpn","Ncoa3","Nfia","Nfic","Nucks1","Nufip2","Nup37","Odf2","Opn3","Pak1Ip1","Pbk","Pcf11","Plin3","Ppp2Ca","Ppp2R2A","Ppp6R3","Prc1","Psen1","Ptms","Pttg1","Rad21","Ran","Rheb","Rpl13A","Slc39A10","Snupn","Srsf3","Stag1","Syncrip","Taf9","Tcerg1","Tle3","Tmem138","Tob2","Top1","Troap","Tsc22D1","Tulp4","Ube2D3","Vangl1","Vcl","Wipf2","Wwc1","Yy1","Zbtb7A","Zcchc10","Zfp24","Zfp281","Zfp593")

print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))

G1S=G1S[which(G1S %in% rownames(dge@data))]
S=S[which(S %in% rownames(dge@data))]
G2M=G2M[which(G2M %in% rownames(dge@data))]
M=M[which(M %in% rownames(dge@data))]
MG1=MG1[which(MG1 %in% rownames(dge@data))]
print(c(length(G1S),length(S),length(G2M),length(M),length(MG1)))

labels=c("G1-S","S","G2-M","M","M-G1")
cellcyclegenes=list(G1S,S,G2M,M,MG1)

exp1=colSums(expm1(dge@data[cellcyclegene[[1]], ]))/colSums(expm1(dge@data))
exp2=colSums(expm1(dge@data[cellcyclegene[[2]], ]))/colSums(expm1(dge@data))
exp3=colSums(expm1(dge@data[cellcyclegene[[3]], ]))/colSums(expm1(dge@data))
exp4=colSums(expm1(dge@data[cellcyclegene[[4]], ]))/colSums(expm1(dge@data))
exp5=colSums(expm1(dge@data[cellcyclegene[[5]], ]))/colSums(expm1(dge@data))

## %Cell cycle genes per Cell
for(id in 1:ncluster){
	c1=rbind(c1,summary(exp1[which(dge@ident==id)]))
	c2=rbind(c2,summary(exp2[which(dge@ident==id)]))
	c3=rbind(c3,summary(exp3[which(dge@ident==id)]))
	c4=rbind(c4,summary(exp4[which(dge@ident==id)]))
	c5=rbind(c5,summary(exp5[which(dge@ident==id)]))
}
write.table(cbind(c1,c2,c3,c4,c5),paste0(dgefile,"Cluster_",setsname[j],"_cellcyclegenes.txt"),quote=F,row.names=F,col.names=F,sep="\t")

######### Finding differentially-expressed markers for 9 sertoli subtypes
###### 12.17.2017 differentially-expressed markers
j=2
dge=dgelist[[j]]
res="res.1.2";i=1
dge=SetAllIdent(dge,id=res[i])
markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_UMI20cell15_31clusters_ReCluster",setsname[j],"_",res[i],"_pct0.2_diffpct0.2_thresh2fold_clustersPC1-11_12.18.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")

### 12.17.2017 Update1: select markers based on fold change and p-value only
# change to min.pct=0 and min.diff.pct=0
j=2
dge=dgelist[[j]]
res="res.1.2";i=1
dge=SetAllIdent(dge,id=res[i])
markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0,min.diff.pct=0,thresh.use=log(2),test.use="bimod",do.print = TRUE)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_UMI20cell15_31clusters_ReCluster",setsname[j],"_",res[i],"_pct0_diffpct0_thresh2fold_clustersPC1-11_12.18.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")

###### 3.5.2018 find markers between subtypes
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster45cellsGenes1k.Robj"))
table(dge@ident)
# 4D  4C  3B   1  2B  4B  4A  3A  2A 
# 71 152 132  94  81  61 208 146 122

### 1. markers with different patterns
#2A Vs 2B
#3A Vs 3B
#4A Vs 4B
#4B Vs 4C-4D

markers2=FindMarkers(dge,"2A","2B",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers3=FindMarkers(dge,"3A","3B",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers4A=FindMarkers(dge,"4A","4B",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers4CD=FindMarkers(dge,"4B",c("4C","4D"),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markers2$cluster="2AVs2B";markers2$gene=rownames(markers2)
markers3$cluster="3AVs3B";markers3$gene=rownames(markers3)
markers4A$cluster="4AVs4B";markers4A$gene=rownames(markers4A)
markers4CD$cluster="4BVs4C&D";markers4CD$gene=rownames(markers4CD)

markersneighbor=rbind(markers2,markers3,markers4A,markers4CD)
head(markersneighbor)
table(markersneighbor$cluster)
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)

#  2AVs2B   3AVs3B   4AVs4B 4BVs4C&D 
#     135       48      132      157 
#      61       29       71       80 
#      74       19       61       77

write.table(markersneighbor,paste0(home,"data_DGE/MarkersAll_Diff1_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=T,quote=F,sep="\t")

markersneighbor$cluster=as.character(markersneighbor$cluster)
meanmarker=matrix(0,nrow=dim(markersneighbor)[1],ncol=3)
meanmarker[,1]=markersneighbor$gene
for(markerj in 1:dim(markersneighbor)[1]){
marker=markersneighbor[markerj,6]
cluster=NULL
clusters=strsplit(markersneighbor[markerj,5],"Vs")[[1]]
cluster[[1]]=strsplit(clusters[1],"&")[[1]]
cluster[[2]]=strsplit(clusters[2],"&")[[1]]
meanmarker[markerj,2]=expMean(dge@data[marker,which(dge@ident %in% cluster[[1]])])
meanmarker[markerj,3]=expMean(dge@data[marker,which(dge@ident %in% cluster[[2]])])
print(markerj)
}
markersneighbormean=cbind(markersneighbor,meanmarker)
which(markersneighbormean[,6] != markersneighbormean[,7])# make sure nothing
write.table(markersneighbormean[,c(6,1:5,8,9)],paste0(home,"data_DGE/MarkersAll_Diff1Mean_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

exps=data.frame(exp1=as.numeric(meanmarker[,2]),exp2=as.numeric(meanmarker[,3]))
exp=ifelse(exps$exp1>exps$exp2,exps$exp1,exps$exp2)
summary(exp)
highexp=markersneighbormean[which(exp>1.2),]
write.table(highexp[,c(6,1:5,8,9)],paste0(home,"data_DGE/MarkersAll_Diff1MeanHighExp_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
highexp1=highexp

table(highexp$cluster)
table(highexp[highexp$avg_diff>0,]$cluster)
table(highexp[highexp$avg_diff<0,]$cluster)
#  2AVs2B   3AVs3B   4AVs4B 4BVs4C&D 
#      87       34       85       94 
#      42       22       41       52 
#      45       12       44       42 


### overlapped markers among 2A-3A-4A or 2B-3B-4B
# overlapped markers for 4B between 4AVs4B and 4BVs4C&D
B41=highexp[which(highexp$cluster=="4AVs4B" & highexp$avg_diff<0),]
B42=highexp[which(highexp$cluster=="4BVs4C&D" & highexp$avg_diff>0),]
length(B41$gene) # 44
length(B42$gene) # 52
length(which(B41$gene %in% B42$gene)) # 17

# overlapped markers for 2A-3A-4A between 2AVs2B, 3AVs3B, 4AVs4B
A2=highexp[which(highexp$cluster=="2AVs2B" & highexp$avg_diff>0),]
A3=highexp[which(highexp$cluster=="3AVs3B" & highexp$avg_diff>0),]
A4=highexp[which(highexp$cluster=="4AVs4B" & highexp$avg_diff>0),]
length(A2$gene) 
length(A3$gene) 
length(A4$gene) 
Reduce(intersect, list(A2$gene,A3$gene,A4$gene)) # 2 
# [1] "Lgals1"  "Tspan17"

# overlapped markers for 2B-3B-4B between 2AVs2B, 3AVs3B, 4AVs4B
B2=highexp[which(highexp$cluster=="2AVs2B" & highexp$avg_diff<0),]
B3=highexp[which(highexp$cluster=="3AVs3B" & highexp$avg_diff<0),]
B4=highexp[which(highexp$cluster=="4AVs4B" & highexp$avg_diff<0),]
length(B2$gene) 
length(B3$gene) 
length(B4$gene) 
Reduce(intersect, list(B2$gene,B3$gene,B4$gene)) # 2 
# [1] "Lgals1"  "Tspan17"


### 2. markers with same patterns
#StageI-III:   4C-4D
#StageIX-XII:  2A-3A-4A
#StageIV-VI:   1-2B-3B-4B
#StageIV-VI:   1-2B

#2A-3A-4A Vs 1-2B-3B-4B
#4B Vs 2A-3A-4A
#4C-4D Vs 1-2B-4B-3B
#4D-4C-3B Vs 1-2B-4B
#4D-4C-3B Vs 1-2B
#1-2B Vs 4B
#1-2B Vs 4A

markersA=FindMarkers(dge,c("2A","3A","4A"),c("1","2B","3B","4B"),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markersA2=FindMarkers(dge,c("2A","3A","4A"),"4B",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers4CDa=FindMarkers(dge,c("4D","4C"),c("1","2B","3B","4B"),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers4CDb=FindMarkers(dge,c("4D","4C","3B"),c("1","2B","4B"),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers1b=FindMarkers(dge,c("1","2B"),"4B",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers1a=FindMarkers(dge,c("1","2B"),"4A",min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markersA$cluster="2A&3A&4AVs1&2B&3B&4B";markersA$gene=rownames(markersA)
markersA2$cluster="2A&3A&4AVs4B";markersA2$gene=rownames(markersA2)
markers4CDa$cluster="4D&4CVs1&2B&3B&4B";markers4CDa$gene=rownames(markers4CDa)
markers4CDb$cluster="4D&4C&3BVs1&2B&4B";markers4CDb$gene=rownames(markers4CDb)
markers1b$cluster="1&2BVs4B";markers1b$gene=rownames(markers1b)
markers1a$cluster="1&2BVs4A";markers1a$gene=rownames(markers1a)

markersneighbor=rbind(markersA,markersA2,markers4CDa,markers4CDb,markers1a,markers1b)
head(markersneighbor)
table(markersneighbor$cluster)
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)

#            1&2BVs4A             1&2BVs4B 2A&3A&4AVs1&2B&3B&4B 
#                 947                 1004                  186 
#                 483                  504                   30 
#                 464                  500                  156 

#        2A&3A&4AVs4B    4D&4C&3BVs1&2B&4B    4D&4CVs1&2B&3B&4B 
#                 118                  454                  183 
#                  39                  198                   37 
#                  79                  256                  146 

write.table(markersneighbor,paste0(home,"data_DGE/MarkersAll_Diff2_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=T,quote=F,sep="\t")

markersneighbor$cluster=as.character(markersneighbor$cluster)
meanmarker=matrix(0,nrow=dim(markersneighbor)[1],ncol=3)
meanmarker[,1]=markersneighbor$gene
for(markerj in 1:dim(markersneighbor)[1]){
marker=markersneighbor[markerj,6]
cluster=NULL
clusters=strsplit(markersneighbor[markerj,5],"Vs")[[1]]
cluster[[1]]=strsplit(clusters[1],"&")[[1]]
cluster[[2]]=strsplit(clusters[2],"&")[[1]]
meanmarker[markerj,2]=expMean(dge@data[marker,which(dge@ident %in% cluster[[1]])])
meanmarker[markerj,3]=expMean(dge@data[marker,which(dge@ident %in% cluster[[2]])])
print(markerj)
}
markersneighbormean=cbind(markersneighbor,meanmarker)
which(markersneighbormean[,6] != markersneighbormean[,7])# make sure nothing
write.table(markersneighbormean[,c(6,1:5,8,9)],paste0(home,"data_DGE/MarkersAll_Diff2Mean_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

exps=data.frame(exp1=as.numeric(meanmarker[,2]),exp2=as.numeric(meanmarker[,3]))
exp=ifelse(exps$exp1>exps$exp2,exps$exp1,exps$exp2)
summary(exp)
highexp=markersneighbormean[which(exp>1.2),]
write.table(highexp[,c(6,1:5,8,9)],paste0(home,"data_DGE/MarkersAll_Diff2MeanHighExp_MouseAdultST24_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
highexp2=highexp

table(highexp$cluster)
table(highexp[highexp$avg_diff>0,]$cluster)
table(highexp[highexp$avg_diff<0,]$cluster)
#            1&2BVs4A             1&2BVs4B 2A&3A&4AVs1&2B&3B&4B 
#                 653                  690                  136 
#                 340                  352                   22 
#                 313                  338                  114 

#        2A&3A&4AVs4B    4D&4C&3BVs1&2B&4B    4D&4CVs1&2B&3B&4B 
#                  76                  361                  143 
#                  18                  154                   37 
 #                 58                  207                  106 

###### 3.15.2018 calculate p-values for all genes for 9 sertoli subtypes
j=2
dge=dgelist[[j]]

### one Vs all others
markersvall=FindAllMarkers(dge,only.pos=FALSE,min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells = -Inf, return.thresh=Inf,test.use="bimod",do.print = TRUE)
# include negative markers by only.pos=FALSE, no fold-change cutoff
write.table(markersvall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_volcano_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
markersvall1=markersvall

table(markersvall$cluster)
#   4D    4C    3B     1    2B    4B    4A    3A    2A
#21219 21219 21219 21219 21219 21219 21219 21219 21219

markersv2A=markersvall[which(markersvall$cluster=="2A"),]
markersv2B=markersvall[which(markersvall$cluster=="2B"),]
markersv3A=markersvall[which(markersvall$cluster=="3A"),]
markersv3B=markersvall[which(markersvall$cluster=="3B"),]
markersv4A=markersvall[which(markersvall$cluster=="4A"),]
markersv4B=markersvall[which(markersvall$cluster=="4B"),]

dim(markersv2A)
dim(markersv2B)
dim(markersv3A)
dim(markersv3B)
dim(markersv4A)
dim(markersv4B)

markersv4CD=FindMarkers(dge,ident.1=c("4C","4D"),only.pos=FALSE,min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells = -Inf, test.use="bimod")
markersv4CD$cluster="4C4D";markersv4CD$gene=rownames(markersv4CD)
dim(markersv4CD)

### one Vs one
markersv2=FindMarkers(dge,"2A","2B",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv3=FindMarkers(dge,"3A","3B",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv4=FindMarkers(dge,"4A","4B",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv4ACD=FindMarkers(dge,"4A",c("4C","4D"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv4BCD=FindMarkers(dge,"4B",c("4C","4D"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")

markersv2$cluster="2AVs2B";markersv2$gene=rownames(markersv2)
markersv3$cluster="3AVs3B";markersv3$gene=rownames(markersv3)
markersv4$cluster="4AVs4B";markersv4$gene=rownames(markersv4)
markersv4ACD$cluster="4AVs4C4D";markersv4ACD$gene=rownames(markersv4ACD)
markersv4BCD$cluster="4BVs4C4D";markersv4BCD$gene=rownames(markersv4BCD)

write.table(markersv2,paste0(home,"data_DGE/MouseAdultST24_Markers2_volcano_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
write.table(markersv3,paste0(home,"data_DGE/MouseAdultST24_Markers3_volcano_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
write.table(markersv4,paste0(home,"data_DGE/MouseAdultST24_Markers4_volcano_UMI20cell15_31clusters_ReCluster45cellsGenes1k_res.1.2_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

dim(markersv2) # [1] 21219     6
dim(markersv3) # [1] 21219     6
dim(markersv4) # [1] 21219     6
dim(markersv4ACD) # [1] 21219     6
dim(markersv4BCD) # [1] 21219     6

### Combine matrix
markersvin=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(markersv2,markersv2A,markersv2B,markersv3,markersv3A,markersv3B,markersv4,markersv4A,markersv4B))
write.table(markersvin,paste0(home,"data_DGE/MouseAdultST24_MarkersIn_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

markersvinmean=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(meanmarker,markersv2,markersv2A,markersv2B,markersv3,markersv3A,markersv3B,markersv4,markersv4ACD,markersv4BCD,markersv4A,markersv4B,markersv4CD))
write.table(markersvin,paste0(home,"data_DGE/MouseAdultST24_MarkersIn_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.19.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

### calculate average expression for each gene
clusters=list("2A","2B","3A","3B","4A","4B","4C","4D",c("4C","4D"))
clusterid=c("2A","2B","3A","3B","4A","4B","4C","4D","4C4D")

meanmarker=matrix(0,nrow=dim(dge@data)[1],ncol=length(clusterid))
colnames(meanmarker)=clusterid
for(markerj in 1:dim(dge@data)[1]){
marker=rownames(dge@data)[markerj]
for(id in 1:length(clusterid)){
  cluster=clusters[[id]]
  tmp=dge@data[marker,which(dge@ident %in% cluster)]
  meanmarker[markerj,id]=expMean(tmp)
}
print(markerj)
}
meanmarker=data.frame(meanmarker)
meanmarker$gene=rownames(dge@data)
colnames(meanmarker)=c(clusterid,"gene")
dim(meanmarker)
write.table(meanmarker,paste0(home,"data_DGE/MouseAdultST24_MarkersMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.15.2018.txt"),col.names=T,row.names=T,quote=F,sep="\t")

### Combine matrix
markersvinmean=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(meanmarker,markersv2,markersv2A,markersv2B,markersv3,markersv3A,markersv3B,markersv4,markersv4A,markersv4B))
write.table(markersvinmean,paste0(home,"data_DGE/MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.15.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

markersvinmean=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(meanmarker,markersv2,markersv2A,markersv2B,markersv3,markersv3A,markersv3B,markersv4,markersv4ACD,markersv4BCD,markersv4A,markersv4B,markersv4CD))
dim(markersvinmean) # [1] 21219    70
markersvinmean[1,]
markersvinmean[1,c(13,14,15,20,25,28,29,30,35,40,43,44,45,48,49,50,53,54,55,60,65,70)]
#  pct.1.x pct.2.x cluster.x cluster.y cluster.x.1 pct.1.y pct.2.y cluster.y.1
#1       0       0    2AVs2B        2A          2B   0.007       0      3AVs3B
#  cluster.x.2 cluster.y.2 pct.1.x.1 pct.2.x.1 cluster.x.3 pct.1.y.1 pct.2.y.1
#1          3A          3B     0.014     0.016      4AVs4B     0.014     0.009
#  cluster.y.3    p_val.x avg_diff.x cluster.x.4 cluster.y.4 cluster.x.5
#1    4AVs4C4D 0.04046971 0.08087893    4BVs4C4D          4A          4B
#  cluster.y.5
#1        4C4D

markersvinmean=markersvinmean[,-c(13,14,15,20,25,28,29,30,35,40,43,44,45,48,49,50,53,54,55,60,65,70)]
dim(markersvinmean) # [1] 21219    48
markersvinmean[1,]
colnames(markersvinmean)=c("gene","Exp.2A","Exp.2B","Exp.3A","Exp.3B","Exp.4A","Exp.4B","Exp.4C","Exp.4D","Exp.4C4D",
"p.2Avs2B","diff.2Avs2B","p.2A","diff.2A","pct.2A","pct.no2A","p.2B","diff.2B","pct.2B","pct.no2B",
"p.3Avs3B","diff.3Avs3B","p.3A","diff.3A","pct.3A","pct.no3A","p.3B","diff.3B","pct.3B","pct.no3B",
"p.4Avs4B","diff.4Avs4B","p.4Avs4C4D","diff.4Avs4C4D","p.4Bvs4C4D","diff.4Bvs4C4D",
"p.4A","diff.4A","pct.4A","pct.no4A","p.4B","diff.4B","pct.4B","pct.no4B","p.4C4D","diff.4C4D","pct.4C4D","pct.no4C4D"
)
markersvinmean[1,]
write.table(markersvinmean,paste0(home,"data_DGE/MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.19.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")


### plot
setwd("C:/Users/qzm/Desktop/AdultMouseST/plot/figNov2017_MouseAdultST24_ReCluster45largecells")
a=read.table("MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.19.2018.txt",header=T)
library(scatterplot3d)

### Exp
s3d <- scatterplot3d(a$Exp.2A,a$Exp.3A,a$Exp.4A,xlab="Exp2A",ylab="Exp3A",zlab="Exp4A",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.2A>2.5 & a$Exp.3A>2.5 & a$Exp.4A>2.5),]
dim(aa) # 20
text(s3d$xyz.convert(aa$Exp.2A,aa$Exp.3A,aa$Exp.4A),col="red", labels=aa$gene,pos=4,cex=0.8) 

### Exp Vs %
plot(a$Exp.2A,a$pct.2A,cex=.5,pch=16,col=rgb(0,0,0,0.6),xlab="log Average Exp",ylab="% Cells")
points(a$Exp.2B,a$pct.2B,cex=.5,pch=16,col=rgb(1,0,0,0.5))
points(a$Exp.3A,a$pct.3A,cex=.5,pch=16,col=rgb(0,1,0,0.4))
points(a$Exp.3B,a$pct.3B,cex=.5,pch=16,col=rgb(0,0,1,0.3))
legend("bottomright",legend=c("2A","2B","3A","3B"),pch=16,col=c("black","red","green","blue"))

### Exp, fold change and p-value
s3d <- scatterplot3d(a$diff.2A,a$Exp.2A,-log10(a$p.2A),xlab="log Fold change in 2A",ylab="Average Exp in 2A",zlab="-log10(P) for 2A Vs Others",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.2A>3.5 & a$diff.2A > 0.7),]
dim(aa) # 20
text(s3d$xyz.convert(aa$diff.2A,aa$Exp.2A,-log10(aa$p.2A)),col="red", labels=aa$gene,pos=4,cex=0.8) 

### Exp 4C, 4D Vs 4C4D
library(scatterplot3d)
s3d <- scatterplot3d(a$Exp.4C,a$Exp.4D,a$Exp.4C4D,xlab="Exp4C",ylab="Exp4D",zlab="Exp4C4D",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.4C>2.5 & a$Exp.4D>2.5 & a$Exp.4C4D>2.5),]
dim(aa) # 49

s3d <- scatterplot3d(a$X2B,a$X3B,a$X4B,xlab="Exp2B",ylab="Exp3B",zlab="Exp4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$X2B>2.5 & a$X3B>2.5 & a$X4B>2.5),]
dim(aa) # 18
text(s3d$xyz.convert(aa$X2B,aa$X3B,aa$X4B), labels=aa$gene,pos=4,cex=0.8) 

s3d <- scatterplot3d(a$pct.1.x,a$pct.1.y.1,a$pct.1.x.3,xlab="%2A",ylab="%3A",zlab="%4A",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$pct.1.x>0.8 & a$pct.1.y.1>0.8 & a$pct.1.x.3>0.8),]
dim(aa) # 24
text(s3d$xyz.convert(aa$pct.1.x,aa$pct.1.y.1,aa$pct.1.x.3), labels=aa$gene,pos=4,cex=0.5) 

s3d <- scatterplot3d(a$pct.2.x,a$pct.2.y.1,a$pct.2.x.3,xlab="%2B",ylab="%3B",zlab="%4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$pct.2.x>0.8 & a$pct.2.y.1>0.8 & a$pct.2.x.3>0.8),]
dim(aa) # 23
text(s3d$xyz.convert(aa$pct.2.x,aa$pct.2.y.1,aa$pct.2.x.3), labels=aa$gene,pos=4,cex=0.5) 

### Fold change 2Avs2B, 3Avs3B, and 4Avs4B
s3d <- scatterplot3d(a$diff.2Avs2B,a$diff.3Avs3B,a$diff.4Avs4B,main="log fold change",xlab="2AVs2B",ylab="3AVs3B",zlab="4AVs4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$diff.2Avs2B>.5 & a$diff.3Avs3B>.5 & a$diff.4Avs4B>.5),]
dim(aa) # 18
text(s3d$xyz.convert(aa$diff.2Avs2B,aa$diff.3Avs3B,aa$diff.4Avs4B), col="red",labels=aa$gene,pos=4,cex=0.8) 
aa$gene
# [1] Ddx52    Dhrs7    Fez2     Fkbp4    Fnbp1    Gm2a    
# [7] Lgals1   Mtap7d3  Pfkl     Plp1     Ppm1m    Rora    
#[13] Sept7    Sfxn1    Slc25a13 Slc25a4  Top1     Tspan17

### p-value in 2Avs2B, 3Avs3B, and 4Avs4B
s3d <- scatterplot3d(-log10(a$p.2Avs2B),-log10(a$p.3Avs3B),-log10(a$p.4Avs4B),main="-log10(p-val)",xlab="2AVs2B",ylab="2AVs3B",zlab="4AVs4B",pch=16,color=rgb(0,0,0,0.5))
bb=a[which(-log10(a$p.2Avs2B)>10 & -log10(a$p.3Avs3B)>3 & -log10(a$p.2Avs2B)>10),]
dim(bb) # 9
text(s3d$xyz.convert(-log10(bb$p.2Avs2B),-log10(bb$p.3Avs3B),-log10(bb$p.4Avs4B)), col="red",labels=bb$gene,pos=4,cex=0.8) 
bb$gene
bb$gene
#[1] Fkbp4     Gm2a      Lgals1    Ppm1m     Ptgds    
#[6] Ptgr1     Rhox5     Serpina3a Slc25a4  

### fold change in 2Avs2B, 2AvsOthers, and 2BvsOthers
aa=a[which(a$p.2Avs2B<0.01 & a$diff.2Avs2B>.8 & a$p.2A<0.01 & a$diff.2A>.8 & a$Exp.2A > 1.5),]
dim(aa) # [1] 10 48
aa$gene
# [1] 1700019B21Rik 8030411F24Rik Bex4          Defb36       
# [5] Gm3880        Lgals1        Ndufv3        Rpl14        
# [9] Slc25a4       Sys1  

bb=a[which(a$p.2Avs2B<0.01 & a$diff.2Avs2B< -0.8 & a$diff.2B>.8  & a$p.2B<0.01 & a$Exp.2B > 1.5),]
dim(bb) # 8
bb$gene # 
#[1] BC051142 Cast     Cypt4    Eppin    Kif2b    Ptgds   
#[7] Rhox5    Stmn1

plot(a$diff.2Avs2B,a$diff.2A,pch=16,cex=a$Exp.2A/2,col=rgb(0,0,0,0.5),xlab="log Fold Change in 2A Vs 2B",ylab="Fold change in Vs Others")
points(a$diff.2Avs2B,-(a$diff.2B),pch=16,cex=a$Exp.2B/2,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("2A vs Others","Others vs 2B"),pch=16,col=c("black","red"))
text(aa$diff.2Avs2B,aa$diff.2A,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.2Avs2B,-(bb$diff.2B),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 3Avs3B, 3AvsOthers, and 3BvsOthers
aa=a[which(a$p.3Avs3B<0.01 & a$diff.3Avs3B>.7 & a$p.3A<0.01 & a$diff.3A>.7 & a$Exp.3A > 1.5),]
dim(aa) # [1] 3
aa$gene
#[1] Caskin1 Esyt3   Gpatch8 Itpr2   Mtap7d3 Sh3d19  Syne1

bb=a[which(a$p.3Avs3B<0.01 & a$diff.3Avs3B< -0.8 & a$diff.3B>.8  & a$p.3B<0.01 & a$Exp.3B > 1.5),]
dim(bb) # 5
bb$gene # 
#[1] Ddx26b Dnm3   Dst    Mical2 Ptprv 

plot(a$diff.3Avs3B,a$diff.3A,pch=16,cex=a$Exp.3A/3,col=rgb(0,0,0,0.5),xlab="log Fold Change in 3A Vs 3B",ylab="Fold change in Vs Others")
points(a$diff.3Avs3B,-(a$diff.3B),pch=16,cex=a$Exp.3B/3,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("3A vs Others","Others vs 3B"),pch=16,col=c("black","red"))
text(aa$diff.3Avs3B,aa$diff.3A,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.3Avs3B,-(bb$diff.3B),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 4Avs4B, 4AvsOthers, and 4BvsOthers
aa=a[which(a$p.4Avs4B<0.01 & a$diff.4Avs4B>.7 & a$p.4A<0.01 & a$diff.4A>.7 & a$Exp.4A > 1.5),]
dim(aa) # [1] 6
aa$gene
#[1] Acsl4   Gm2a    Pla2g7  Plp1    Rora    Slc25a4

bb=a[which(a$p.4Avs4B<0.01 & a$diff.4Avs4B< -0.8 & a$diff.4B>.8  & a$p.4B<0.01 & a$Exp.4B > 1.5),]
dim(bb) # 14
bb$gene # 
# [1] Ano1      Caskin1   Chp2      Ctsl      Dhcr24   
# [6] Drd4      Erp29     Mfge8     Nid1      Prnp     
#[11] Ptgds     Rhox5     Serpina3a Serpina5

plot(a$diff.4Avs4B,a$diff.4A,pch=16,cex=a$Exp.4A/4,col=rgb(0,0,0,0.5),xlab="log Fold Change in 4A Vs 4B",ylab="Fold change in Vs Others")
points(a$diff.4Avs4B,-(a$diff.4B),pch=16,cex=a$Exp.4B/4,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("4A vs Others","Others vs 4B"),pch=16,col=c("black","red"))
text(aa$diff.4Avs4B,aa$diff.4A,aa$gene,col="blue",pos=3,cex=0.8)
text(bb$diff.4Avs4B,-(bb$diff.4B),bb$gene,col="blue",pos=3,cex=0.8)

### fold change in 4Avs4C4D, 4AvsOthers, and 4C4DvsOthers
aa=a[which(a$p.4Avs4C4D<0.01 & a$diff.4Avs4C4D>.7 & a$p.4A<0.01 & a$diff.4A>.7 & a$Exp.4A > 1.5),]
dim(aa) # [1] 5
aa$gene
#[1] Pla2g7  Rora    Slc25a4 Timp3   Uchl1

bb=a[which(a$p.4Avs4C4D<0.01 & a$diff.4Avs4C4D< -0.8 & a$diff.4C4D>.8  & a$p.4C4D<0.01 & a$Exp.4C4D > 1.5),]
dim(bb) # 5
bb$gene # 
#[1] Aqp8      Fyb       Hmgn5     Itm2c     Serpina3a

plot(a$diff.4Avs4C4D,a$diff.4A,pch=16,cex=a$Exp.4A/4,col=rgb(0,0,0,0.5),xlab="log Fold Change in 4A Vs 4C4D",ylab="Fold change in Vs Others")
points(a$diff.4Avs4C4D,-(a$diff.4C4D),pch=16,cex=a$Exp.4C4D/4,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("4A vs Others","Others vs 4C4D"),pch=16,col=c("black","red"))
text(aa$diff.4Avs4C4D,aa$diff.4A,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.4Avs4C4D,-(bb$diff.4C4D),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 4Bvs4C4D, 4BvsOthers, and 4C4DvsOthers
aa=a[which(a$p.4Bvs4C4D<0.01 & a$diff.4Bvs4C4D>.8 & a$p.4B<0.01 & a$diff.4B>.8 & a$Exp.4B > 1.5),]
dim(aa) # [1] 14
aa$gene
# [1] Ano1    Arl6ip1 Caskin1 Cdkn1a  Chp2    Drd4    Erp29  
# [8] Mfge8   Nid1    P2rx2   Prnp    Ptgds   Tnni3   Wfdc10

bb=a[which(a$p.4Bvs4C4D<0.01 & a$diff.4Bvs4C4D< -0.7 & a$diff.4C4D>.7  & a$p.4C4D<0.01 & a$Exp.4C4D > 1.5),]
dim(bb) # 5
bb$gene # 
#[1] Fyb     Gas6    Itm2c   Slc26a2 Sparc

plot(a$diff.4Bvs4C4D,a$diff.4B,pch=16,cex=a$Exp.4B/4,col=rgb(0,0,0,0.5),xlab="log Fold Change in 4B Vs 4C4D",ylab="Fold change in Vs Others")
points(a$diff.4Bvs4C4D,-(a$diff.4C4D),pch=16,cex=a$Exp.4C4D/4,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("4B vs Others","Others vs 4C4D"),pch=16,col=c("black","red"))
text(aa$diff.4Bvs4C4D,aa$diff.4B,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.4Bvs4C4D,-(bb$diff.4C4D),bb$gene,col="blue",pos=1,cex=0.8)


###### 5.7.2018 build complete list for marker selection for 9 sertoli subtypes
#both 2A&2B Vs all others
#both 3A&3B Vs all others
#both 4A&4B Vs all others
#both 4C&4D Vs all others
#        2A  2B  3A  3B
#2Avs2B  ++  --  +   -
#2Avs3A  ++  +   --  -
#-----------------------
#3AVs3B  +   -   ++  --

j=2
dge=dgelist[[j]]
### one Vs all others
markersv2AB=FindMarkers(dge,ident.1=c("2A","2B"),only.pos=FALSE,min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells = -Inf, test.use="bimod")
markersv2AB$cluster="2A2B";markersv2AB$gene=rownames(markersv2AB)
dim(markersv2AB)
markersv3AB=FindMarkers(dge,ident.1=c("3A","3B"),only.pos=FALSE,min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells = -Inf, test.use="bimod")
markersv3AB$cluster="3A3B";markersv3AB$gene=rownames(markersv3AB)
dim(markersv3AB)
markersv4AB=FindMarkers(dge,ident.1=c("4A","4B"),only.pos=FALSE,min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells = -Inf, test.use="bimod")
markersv4AB$cluster="4A4B";markersv4AB$gene=rownames(markersv4AB)
dim(markersv4AB)

### one Vs one
markersvA=FindMarkers(dge,"2A","3A",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersvB=FindMarkers(dge,"2B","3B",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersvA2=FindMarkers(dge,"2A","4A",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersvB2=FindMarkers(dge,"2B","4B",min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv23=FindMarkers(dge,c("2A","2B"),c("3A","3B"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv24=FindMarkers(dge,c("2A","2B"),c("4A","4B"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv34=FindMarkers(dge,c("3A","3B"),c("4A","4B"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")
markersv4=FindMarkers(dge,c("4A","4B"),c("4C","4D"),min.pct=-Inf,min.diff.pct=-Inf,thresh.use=-Inf,min.cells =-Inf,test.use="bimod")

markersvA$cluster="2AVs3A";markersvA$gene=rownames(markersvA)
markersvB$cluster="2BVs3B";markersvB$gene=rownames(markersvB)
markersvA2$cluster="2AVs4A";markersvA2$gene=rownames(markersvA2)
markersvB2$cluster="2BVs4B";markersvB2$gene=rownames(markersvB2)
markersv23$cluster="2A2BVs3A3B";markersv23$gene=rownames(markersv23)
markersv24$cluster="2A2BVs4A4B";markersv24$gene=rownames(markersv24)
markersv34$cluster="3A3BVs4A4B";markersv34$gene=rownames(markersv34)
markersv4$cluster="4A4BVs4C4D";markersv4$gene=rownames(markersv4)

dim(markersvA) # [1] 21219     6
dim(markersvB) # [1] 21219     6
dim(markersvA2) # [1] 21219     6
dim(markersvB2) # [1] 21219     6
dim(markersv23) # [1] 21219     6
dim(markersv24) # [1] 21219     6
dim(markersv34) # [1] 21219     6
dim(markersv4) # [1] 21219     6

### Combine matrix
markersvin=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(markersv2AB,markersv3AB,markersv4AB,markersvA,markersvB,markersvA2,markersvB2,markersv23,markersv24,markersv34,markersv4))
write.table(markersvin,paste0(home,"data_DGE/MouseAdultST24_MarkersInAB_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

markersvin2=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(markersvA2,markersvB2,markersv23,markersv24,markersv34,markersv4))
write.table(markersvin2,paste0(home,"data_DGE/MouseAdultST24_MarkersInAB2_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
dim(markersvin2) # [1] 21219    31
markersvin2[1,]
markersvin2[1,c(1:3,7:8,12:13,17:18,22:23,27:28)]


markersvin2=markersvin2[,c(1:3,7:8,12:13,17:18,22:23,27:28)]
dim(markersvin2) # [1] 21219    13
markersvin2[1,]
colnames(markersvin2)=c("gene",
"p.2Avs4A","diff.2Avs4A","p.2Bvs4B","diff.2Bvs4B",
"p.2A2Bvs3A3B","diff.2A2Bvs3A3B","p.2A2Bvs4A4B","diff.2A2Bvs4A4B",
"p.3A3Bvs4A4B","diff.3A3Bvs4A4B","p.4A4Bvs4C4D","diff.4A4Bvs4C4D"
)
markersvin2[1,]
write.table(markersvin2,paste0(home,"data_DGE/MouseAdultST24_MarkersInAB2_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")


### calculate average expression for each gene
clusters=list(c("2A","2B"),c("3A","3B"),c("4A","4B"))
clusterid=c("2A2B","3A3B","4A4B")

meanmarker=matrix(0,nrow=dim(dge@data)[1],ncol=length(clusterid))
colnames(meanmarker)=clusterid
for(markerj in 1:dim(dge@data)[1]){
marker=rownames(dge@data)[markerj]
for(id in 1:length(clusterid)){
  cluster=clusters[[id]]
  tmp=dge@data[marker,which(dge@ident %in% cluster)]
  meanmarker[markerj,id]=expMean(tmp)
}
print(markerj)
}
meanmarker=data.frame(meanmarker)
meanmarker$gene=rownames(dge@data)
colnames(meanmarker)=c(clusterid,"gene")
dim(meanmarker)
write.table(meanmarker,paste0(home,"data_DGE/MouseAdultST24_MeanAB_volcano_ReCluster45cellsGenes1k_res.1.2_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")


### Combine matrix
markersvinmean2=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(meanmarker,markersv2AB,markersv3AB,markersv4AB,markersvA,markersvB))
write.table(markersvinmean2,paste0(home,"data_DGE/MouseAdultST24_MarkersInMeanAB_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

dim(markersvinmean2) # [1] 21219    29
markersvinmean2[1,]
markersvinmean2[1,-c(9,14,19,22:24,27:29)]
#           gene 2A2B        3A3B       4A4B   p_val.x  avg_diff.x pct.1.x
#1 0610005C13Rik    0 0.008950384 0.04132305 0.3976381 -0.01828608       0
#  pct.2.x cluster.x p_val.y   avg_diff.y pct.1.y pct.2.y cluster.y p_val.x.1
#1   0.008      2A2B       1 -0.007946779   0.004   0.008      3A3B 0.1543284
#  avg_diff.x.1 pct.1.x.1 pct.2.x.1 cluster.x.1 p_val.y.1 avg_diff.y.1
#1    0.0355805     0.015     0.004        4A4B 0.7493013  -0.01697413
#  cluster.y.1 p_val avg_diff cluster
#1      2AVs3A     1        0  2BVs3B

markersvinmean2=markersvinmean2[,-c(9,14,19,22:24,27:29)]
dim(markersvinmean2) # [1] 21219    20
markersvinmean2[1,]
colnames(markersvinmean2)=c("gene","Exp.2A2B","Exp.3A3B","Exp.4A4B",
"p.2A2B","diff.2A2B","pct.2A2B","pct.no2A2B",
"p.3A3B","diff.3A3B","pct.3A3B","pct.no3A3B",
"p.4A4B","diff.4A4B","pct.4A4B","pct.no4A4B",
"p.2Avs3A","diff.2Avs3A","p.2Bvs3B","diff.2Bvs3B"
)
markersvinmean2[1,]
write.table(markersvinmean2,paste0(home,"data_DGE/MouseAdultST24_MarkersInMeanAB_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

markersvinmean=read.table(paste0(home,"data_DGE/MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_3.19.2018.txt"),header=T,row.names=1)
a=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(markersvinmean,markersvinmean2))
dim(a) # 67
a[1,]
a=a[,c(1,2,3,49,4,5,50,6,7,51,8:20,52:55,21:30,56:59,64:67,31:44,60:63,45:48)]
write.table(a,paste0(home,"data_DGE/MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
a=Reduce(function(x,y)merge(x,y,by="gene",all=TRUE),list(a,markersvin2))
write.table(a,paste0(home,"data_DGE/MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")
dim(a) # 79
colnames(a)
# [1] "gene"          "Exp.2A"        "Exp.2B"        "Exp.2A2B"
# [5] "Exp.3A"        "Exp.3B"        "Exp.3A3B"      "Exp.4A"
# [9] "Exp.4B"        "Exp.4A4B"      "Exp.4C"        "Exp.4D"
#[13] "Exp.4C4D"      "p.2Avs2B"      "diff.2Avs2B"   "p.2A"
#[17] "diff.2A"       "pct.2A"        "pct.no2A"      "p.2B"
#[21] "diff.2B"       "pct.2B"        "pct.no2B"      "p.2A2B"
#[25] "diff.2A2B"     "pct.2A2B"      "pct.no2A2B"    "p.3Avs3B"
#[29] "diff.3Avs3B"   "p.3A"          "diff.3A"       "pct.3A"
#[33] "pct.no3A"      "p.3B"          "diff.3B"       "pct.3B"
#[37] "pct.no3B"      "p.3A3B"        "diff.3A3B"     "pct.3A3B"
#[41] "pct.no3A3B"    "p.2Avs3A"      "diff.2Avs3A"   "p.2Bvs3B"
#[45] "diff.2Bvs3B"   "p.4Avs4B"      "diff.4Avs4B"   "p.4Avs4C4D"
#[49] "diff.4Avs4C4D" "p.4Bvs4C4D"    "diff.4Bvs4C4D" "p.4A"
#[53] "diff.4A"       "pct.4A"        "pct.no4A"      "p.4B"
#[57] "diff.4B"       "pct.4B"        "pct.no4B"      "p.4A4B"
#[61] "diff.4A4B"     "pct.4A4B"      "pct.no4A4B"    "p.4C4D"
#[65] "diff.4C4D"     "pct.4C4D"      "pct.no4C4D"

###### plot
setwd("C:/Users/qzm/Desktop/AdultMouseST/plot/figNov2017_MouseAdultST24_ReCluster45largecells")
a=read.table("MouseAdultST24_MarkersInMean_volcano_ReCluster45cellsGenes1k_res.1.2_combined_5.7.2018.txt",header=T)
library(scatterplot3d)

### Exp
s3d <- scatterplot3d(a$Exp.2A2B,a$Exp.3A3B,a$Exp.4A4B,xlab="Exp2A2B",ylab="Exp3A3B",zlab="Exp4A4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.2A2B>2.5 & a$Exp.3A3B>2.5 & a$Exp.4A4B>2.5),]
dim(aa) # 20
text(s3d$xyz.convert(aa$Exp.2A2B,aa$Exp.3A3B,aa$Exp.4A4B),col="red", labels=aa$gene,pos=4,cex=0.8) 

### Exp Vs %
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(a$Exp.2A2B,a$pct.2A2B,cex=.5,pch=16,col=rgb(0,0,0,0.6),xlab="log Average Exp",ylab="% Cells")
points(a$Exp.3A3B,a$pct.3A3B,cex=.5,pch=16,col=rgb(1,0,0,0.4))
points(a$Exp.4A4B,a$pct.4A4B,cex=.5,pch=16,col=rgb(0,1,0,0.5))
points(a$Exp.4C4D,a$pct.4C4D,cex=.5,pch=16,col=rgb(0,0,1,0.2))
legend("bottomright",legend=c("2A2B","3A3B","4A4B","4C4D"),pch=16,col=c("black","red","green","blue"))

### Exp, fold change and p-value
s3d <- scatterplot3d(a$diff.2A,a$Exp.2A,-log10(a$p.2A),xlab="log Fold change in 2A",ylab="Average Exp in 2A",zlab="-log10(P) for 2A Vs Others",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.2A>3.5 & a$diff.2A > 0.7),]
dim(aa) # 20
text(s3d$xyz.convert(aa$diff.2A,aa$Exp.2A,-log10(aa$p.2A)),col="red", labels=aa$gene,pos=4,cex=0.8) 

### Exp 4C, 4D Vs 4C4D
library(scatterplot3d)
s3d <- scatterplot3d(a$Exp.4C,a$Exp.4D,a$Exp.4C4D,xlab="Exp4C",ylab="Exp4D",zlab="Exp4C4D",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$Exp.4C>2.5 & a$Exp.4D>2.5 & a$Exp.4C4D>2.5),]
dim(aa) # 49


s3d <- scatterplot3d(a$X2B,a$X3B,a$X4B,xlab="Exp2B",ylab="Exp3B",zlab="Exp4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$X2B>2.5 & a$X3B>2.5 & a$X4B>2.5),]
dim(aa) # 18
text(s3d$xyz.convert(aa$X2B,aa$X3B,aa$X4B), labels=aa$gene,pos=4,cex=0.8) 

s3d <- scatterplot3d(a$pct.1.x,a$pct.1.y.1,a$pct.1.x.3,xlab="%2A",ylab="%3A",zlab="%4A",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$pct.1.x>0.8 & a$pct.1.y.1>0.8 & a$pct.1.x.3>0.8),]
dim(aa) # 24
text(s3d$xyz.convert(aa$pct.1.x,aa$pct.1.y.1,aa$pct.1.x.3), labels=aa$gene,pos=4,cex=0.5) 

s3d <- scatterplot3d(a$pct.2.x,a$pct.2.y.1,a$pct.2.x.3,xlab="%2B",ylab="%3B",zlab="%4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$pct.2.x>0.8 & a$pct.2.y.1>0.8 & a$pct.2.x.3>0.8),]
dim(aa) # 23
text(s3d$xyz.convert(aa$pct.2.x,aa$pct.2.y.1,aa$pct.2.x.3), labels=aa$gene,pos=4,cex=0.5) 

### Fold change 2Avs2B, 3Avs3B, and 4Avs4B
s3d <- scatterplot3d(a$diff.2Avs2B,a$diff.3Avs3B,a$diff.4Avs4B,main="log fold change",xlab="2AVs2B",ylab="3AVs3B",zlab="4AVs4B",pch=16,color=rgb(0,0,0,0.5))
aa=a[which(a$diff.2Avs2B>.5 & a$diff.3Avs3B>.5 & a$diff.4Avs4B>.5),]
dim(aa) # 18
text(s3d$xyz.convert(aa$diff.2Avs2B,aa$diff.3Avs3B,aa$diff.4Avs4B), col="red",labels=aa$gene,pos=4,cex=0.8) 
aa$gene
# [1] Ddx52    Dhrs7    Fez2     Fkbp4    Fnbp1    Gm2a    
# [7] Lgals1   Mtap7d3  Pfkl     Plp1     Ppm1m    Rora    
#[13] Sept7    Sfxn1    Slc25a13 Slc25a4  Top1     Tspan17

### p-value in 2A2B, 3A3B, and 2Avs3A
s3d <- scatterplot3d(-log10(a$p.2A2B),-log10(a$p.3A3B),-log10(a$p.2Avs3A),main="-log10(p-val)",xlab="2A2B",ylab="3A3B",zlab="2AVs3A",pch=16,color=rgb(0,0,0,0.5))
quantile(-log10(a$p.2A2B),.95)
quantile(-log10(a$p.3A3B),.95)
quantile(-log10(a$p.2Avs3A),.95)
bb=a[which(-log10(a$p.2A2B)>50 & -log10(a$p.3A3B)>50 & -log10(a$p.2Avs3A)>25),]
dim(bb) # 9
text(s3d$xyz.convert(-log10(bb$p.2A2B),-log10(bb$p.3A3B),-log10(bb$p.2Avs3A)), col="red",labels=bb$gene,pos=4,cex=0.8) 
bb$gene
bb$gene
#[1] Fkbp4     Gm2a      Lgals1    Ppm1m     Ptgds    
#[6] Ptgr1     Rhox5     Serpina3a Slc25a4  

### fold change in 2Avs3A, 2AvsOthers, and 3AvsOthers
aa=a[which(a$p.2Avs3A<0.01 & a$diff.2Avs3A>1.5 & a$p.2A<0.01 & a$diff.2A>.8 & a$Exp.2A > 1.5),]
dim(aa) # [1] 20 67
aa$gene
# [1] 1700011M02Rik 1700019B21Rik Atox1         Bex4         
# [5] Chchd10       Chchd2        Cycs          Defb19       
# [9] Defb36        Etd           Fth1          Gm3880       
#[13] Gm9112        Minos1        Rpl14         Rplp1        
#[17] Rps12         Rps29         Rps3          Tmsb4x  

bb=a[which(a$p.2Avs3A<0.01 & a$diff.2Avs3A< -1.5 & a$diff.3A>.8  & a$p.3A<0.01 & a$Exp.3A > 1.5),]
dim(bb) # 28
bb$gene # 
# [1] 4632427E13Rik Caskin1       Chka          Clcn2        
# [5] Clk1          Dusp11        Fus           Gpatch8      
# [9] Itpr2         Kcnq1ot1      Kmt2c         Macf1        
#[13] Malat1        Mlxipl        Nfat5         Nktr         
#[17] Pabpn1        Phip          Pnisr         Rbm26        
#[21] Scarb1        Son           Srsf10        Syne1        
#[25] Trank1        Ubn2          Upf3b         Zcchc7

par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(a$diff.2Avs3A,a$diff.2A,pch=16,cex=a$Exp.2A/2,col=rgb(0,0,0,0.5),xlab="log Fold Change in 2A Vs 3A",ylab="Fold change in Vs Others")
points(a$diff.2Avs3A,-(a$diff.3A),pch=16,cex=a$Exp.2B/2,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("2A vs Others","Others vs 3A"),pch=16,col=c("black","red"))
text(aa$diff.2Avs3A,aa$diff.2A,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.2Avs3A,-(bb$diff.3A),bb$gene,col="blue",pos=1,cex=0.8)


### fold change in 2Bvs3B, 2BvsOthers, and 3BvsOthers
aa=a[which(a$p.2Bvs3B<0.01 & a$diff.2Bvs3B>.7 & a$p.2B<0.01 & a$diff.2B>.7 & a$Exp.2B > 1.5),]
dim(aa) # [1] 32
aa$gene
# [1] Aard      Akr1a1    Cast      Chp2      Cypt4     Defb45   
# [7] Eppin     Etd       Gatm      Gsta2     H2afj     Hmgb3    
#[13] Map1lc3b  Marcksl1  Mrfap1    Mt1       Ndufa1    Ptgds    
#[19] Rhox5     Rhox8     Rnf7      Rpl41     Rplp0     Rps21    
#[25] Serpina3a Smpx      Stmn1     Tmsb4x    Tomm7     Tsx      
#[31] Wfdc10    Yipf1

bb=a[which(a$p.2Bvs3B<0.01 & a$diff.2Bvs3B< -0.7 & a$diff.3B>.7  & a$p.3B<0.01 & a$Exp.3B > 1.5),]
dim(bb) # 30
bb$gene # 
# [1] Ccnt2   Daam2   Ddx17   Ddx26b  Dnm3    Dpysl4  Dst    
# [8] Etl4    Fyb     Gm2044  Hdc     Iffo1   Malat1  Mical2 
#[15] Mpi     Nfat5   Nxf3    Pabpn1  Phf20l1 Pnisr   Pnn    
#[22] Ppp1r9a Ptprv   Rsrp1   Snhg20  Son     Sptbn1  Srsf1  
#[29] Tpp2    Trim7 

plot(a$diff.2Bvs3B,a$diff.2B,pch=16,cex=a$Exp.2B/3,col=rgb(0,0,0,0.5),xlab="log Fold Change in 2B Vs 3B",ylab="Fold change in Vs Others")
points(a$diff.2Bvs3B,-(a$diff.3B),pch=16,cex=a$Exp.3B/3,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("2B vs Others","Others vs 3B"),pch=16,col=c("black","red"))
text(aa$diff.2Bvs3B,aa$diff.2B,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.2Bvs3B,-(bb$diff.3B),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 2A2Bvs3A3B, 2A2BvsOthers, and 3A3BvsOthers
aa=a[which(a$p.2A2Bvs3A3B<0.01 & a$diff.2A2Bvs3A3B>.8 & a$p.2A2B<0.01 & a$diff.2A2B>.8 & a$Exp.2A2B > 1.5),]
dim(aa) # [1] 53
aa$gene
# [1] 1500011B03Rik 1700011M02Rik Atox1         Atp5h        
# [5] Atp6v0e       Chchd10       Chchd2        Cited1       
# [9] Cox6a1        Cox6b1        Cox6c         Cox8a        
#[13] Cycs          Dbi           Defb19        Defb45       
#[17] Etd           Fth1          Gm3880        Gm9112       
#[21] Gstp1         Lgals1        Marcksl1      Mt1          
#[25] Myeov2        Myl9          Ndufa1        Ndufa3       
#[29] Ndufa4        Ppia          Rhox5         Romo1        
#[33] Rpl13         Rpl14         Rpl23a        Rpl27a       
#[37] Rpl41         Rplp0         Rplp1         Rplp2        
#[41] Rps12         Rps17         Rps19         Rps27a       
#[45] Rps29         Rps3          Sin3b         Tma7         
#[49] Tmsb4x        Tpt1          Tsx           Ubb          
#[53] Uqcr11 

bb=a[which(a$p.2A2Bvs3A3B<0.01 & a$diff.2A2Bvs3A3B< -0.8 & a$diff.3A3B>.8  & a$p.3A3B<0.01 & a$Exp.3A3B > 1.5),]
dim(bb) # 31
bb$gene # 
# [1] Chka     Clcn2    Clk1     Dusp11   Fus      Garnl3  
# [7] Gm2044   Gpatch8  Kcnq1ot1 Kmt2c    Macf1    Malat1  
#[13] Nfat5    Nisch    Nktr     Pabpn1   Pnisr    Pnn     
#[19] Ppp1r9a  Prpf38b  Rbm26    Rsrp1    Snhg20   Son     
#[25] Srsf10   Trim7    Ttc14    Ubn2     Upf3b    Zcchc7  
#[31] Zkscan3 
 

plot(a$diff.2A2Bvs3A3B,a$diff.2A2B,pch=16,cex=a$Exp.2A2B/3,col=rgb(0,0,0,0.5),xlab="log Fold Change in 2A2B Vs 3A3B",ylab="Fold change in Vs Others")
points(a$diff.2A2Bvs3A3B,-(a$diff.3A3B),pch=16,cex=a$Exp.3A3B/3,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("2A2B vs Others","Others vs 3A3B"),pch=16,col=c("black","red"))
text(aa$diff.2A2Bvs3A3B,aa$diff.2A2B,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.2A2Bvs3A3B,-(bb$diff.3A3B),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 4Avs4C4D, 4AvsOthers, and 4C4DvsOthers
aa=a[which(a$p.4Avs4C4D<0.01 & a$diff.4Avs4C4D>.7 & a$p.4A<0.01 & a$diff.4A>.7 & a$Exp.4A > 1.5),]
dim(aa) # [1] 5
aa$gene
#[1] Pla2g7  Rora    Slc25a4 Timp3   Uchl1

bb=a[which(a$p.4Avs4C4D<0.01 & a$diff.4Avs4C4D< -0.8 & a$diff.4C4D>.8  & a$p.4C4D<0.01 & a$Exp.4C4D > 1.5),]
dim(bb) # 5
bb$gene # 
#[1] Aqp8      Fyb       Hmgn5     Itm2c     Serpina3a

plot(a$diff.4Avs4C4D,a$diff.4A,pch=16,cex=a$Exp.4A/4,col=rgb(0,0,0,0.5),xlab="log Fold Change in 4A Vs 4C4D",ylab="Fold change in Vs Others")
points(a$diff.4Avs4C4D,-(a$diff.4C4D),pch=16,cex=a$Exp.4C4D/4,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("4A vs Others","Others vs 4C4D"),pch=16,col=c("black","red"))
text(aa$diff.4Avs4C4D,aa$diff.4A,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.4Avs4C4D,-(bb$diff.4C4D),bb$gene,col="blue",pos=1,cex=0.8)

### fold change in 4Bvs4C4D, 4BvsOthers, and 4C4DvsOthers
aa=a[which(a$p.4Bvs4C4D<0.01 & a$diff.4Bvs4C4D>.8 & a$p.4B<0.01 & a$diff.4B>.8 & a$Exp.4B > 1.5),]
dim(aa) # [1] 14
aa$gene
# [1] Ano1    Arl6ip1 Caskin1 Cdkn1a  Chp2    Drd4    Erp29  
# [8] Mfge8   Nid1    P2rx2   Prnp    Ptgds   Tnni3   Wfdc10

bb=a[which(a$p.4Bvs4C4D<0.01 & a$diff.4Bvs4C4D< -0.7 & a$diff.4C4D>.7  & a$p.4C4D<0.01 & a$Exp.4C4D > 1.5),]
dim(bb) # 5
bb$gene # 
#[1] Fyb     Gas6    Itm2c   Slc26a2 Sparc

plot(a$diff.4Bvs4C4D,a$diff.4B,pch=16,cex=a$Exp.4B/4,col=rgb(0,0,0,0.5),xlab="log Fold Change in 4B Vs 4C4D",ylab="Fold change in Vs Others")
points(a$diff.4Bvs4C4D,-(a$diff.4C4D),pch=16,cex=a$Exp.4C4D/4,col=rgb(1,0,0,0.5))
legend("bottomright",legend=c("4B vs Others","Others vs 4C4D"),pch=16,col=c("black","red"))
text(aa$diff.4Bvs4C4D,aa$diff.4B,aa$gene,col="blue",pos=1,cex=0.8)
text(bb$diff.4Bvs4C4D,-(bb$diff.4C4D),bb$gene,col="blue",pos=1,cex=0.8)


