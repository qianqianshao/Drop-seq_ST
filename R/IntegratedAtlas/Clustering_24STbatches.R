### R script for analyzing merged 24 ST batches in June 2017 by Qianyi, Related to Figure 1B and 2A

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)

### 24 Adult Mouse ST batches (without INT6) = 6 ST + 2 1n-Depleted + 3 SPG-enriched + 5 INT-enriched + 8 SER-enriched
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
topcells=read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/topcellST24", stringsAsFactors=F)
infile=topcells[,1]
redblue100<-rgb(read.table('/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt'),sep='\t',row.names=1,header=T))

###### 1. Merge all files:
#####  1.1 merge files for mouse cells and genes
### read in files in mouse cells and genes
mouse1=read.table(paste0(home,"data_DGE/",infile[1],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse2=read.table(paste0(home,"data_DGE/",infile[2],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse3=read.table(paste0(home,"data_DGE/",infile[3],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse4=read.table(paste0(home,"data_DGE/",infile[4],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse5=read.table(paste0(home,"data_DGE/",infile[5],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse6=read.table(paste0(home,"data_DGE/",infile[6],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse7=read.table(paste0(home,"data_DGE/",infile[7],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse8=read.table(paste0(home,"data_DGE/",infile[8],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse9=read.table(paste0(home,"data_DGE/",infile[9],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse10=read.table(paste0(home,"data_DGE/",infile[10],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse11=read.table(paste0(home,"data_DGE/",infile[11],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse12=read.table(paste0(home,"data_DGE/",infile[12],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse13=read.table(paste0(home,"data_DGE/",infile[13],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse14=read.table(paste0(home,"data_DGE/",infile[14],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse15=read.table(paste0(home,"data_DGE/",infile[15],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse16=read.table(paste0(home,"data_DGE/",infile[16],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse17=read.table(paste0(home,"data_DGE/",infile[17],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse18=read.table(paste0(home,"data_DGE/",infile[18],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse19=read.table(paste0(home,"data_DGE/",infile[19],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse20=read.table(paste0(home,"data_DGE/",infile[20],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse21=read.table(paste0(home,"data_DGE/",infile[21],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse22=read.table(paste0(home,"data_DGE/",infile[22],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse23=read.table(paste0(home,"data_DGE/",infile[23],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)
mouse24=read.table(paste0(home,"data_DGE/",infile[24],"/mouse_gene_exon_tagged_cleaned5000cells.dge.txt.gz"),header=T, stringsAsFactors=F)

### add batch name to Cell barcode names
names(mouse1)[-1]=paste(topcells[1,3],names(mouse1)[-1],sep="_")
names(mouse2)[-1]=paste(topcells[2,3],names(mouse2)[-1],sep="_")
names(mouse3)[-1]=paste(topcells[3,3],names(mouse3)[-1],sep="_")
names(mouse4)[-1]=paste(topcells[4,3],names(mouse4)[-1],sep="_")
names(mouse5)[-1]=paste(topcells[5,3],names(mouse5)[-1],sep="_")
names(mouse6)[-1]=paste(topcells[6,3],names(mouse6)[-1],sep="_")
names(mouse7)[-1]=paste(topcells[7,3],names(mouse7)[-1],sep="_")
names(mouse8)[-1]=paste(topcells[8,3],names(mouse8)[-1],sep="_")
names(mouse9)[-1]=paste(topcells[9,3],names(mouse9)[-1],sep="_")
names(mouse10)[-1]=paste(topcells[10,3],names(mouse10)[-1],sep="_")
names(mouse11)[-1]=paste(topcells[11,3],names(mouse11)[-1],sep="_")
names(mouse12)[-1]=paste(topcells[12,3],names(mouse12)[-1],sep="_")
names(mouse13)[-1]=paste(topcells[13,3],names(mouse13)[-1],sep="_")
names(mouse14)[-1]=paste(topcells[14,3],names(mouse14)[-1],sep="_")
names(mouse15)[-1]=paste(topcells[15,3],names(mouse15)[-1],sep="_")
names(mouse16)[-1]=paste(topcells[16,3],names(mouse16)[-1],sep="_")
names(mouse17)[-1]=paste(topcells[17,3],names(mouse17)[-1],sep="_")
names(mouse18)[-1]=paste(topcells[18,3],names(mouse18)[-1],sep="_")
names(mouse19)[-1]=paste(topcells[19,3],names(mouse19)[-1],sep="_")
names(mouse20)[-1]=paste(topcells[20,3],names(mouse20)[-1],sep="_")
names(mouse21)[-1]=paste(topcells[21,3],names(mouse20)[-1],sep="_")
names(mouse22)[-1]=paste(topcells[22,3],names(mouse22)[-1],sep="_")
names(mouse23)[-1]=paste(topcells[23,3],names(mouse23)[-1],sep="_")
names(mouse24)[-1]=paste(topcells[24,3],names(mouse24)[-1],sep="_")

# note: data.table is much fater and memory-saving than data.frame, highly recommended!!!
library(data.table)
datatablemouse1=data.table(mouse1,key="GENE")
datatablemouse2=data.table(mouse2,key="GENE")
datatablemouse3=data.table(mouse3,key="GENE")
datatablemouse4=data.table(mouse4,key="GENE")
datatablemouse5=data.table(mouse5,key="GENE")
datatablemouse6=data.table(mouse6,key="GENE")
datatablemouse7=data.table(mouse7,key="GENE")
datatablemouse8=data.table(mouse8,key="GENE")
datatablemouse9=data.table(mouse9,key="GENE")
datatablemouse10=data.table(mouse10,key="GENE")
datatablemouse11=data.table(mouse11,key="GENE")
datatablemouse12=data.table(mouse12,key="GENE")
datatablemouse13=data.table(mouse13,key="GENE")
datatablemouse14=data.table(mouse14,key="GENE")
datatablemouse15=data.table(mouse15,key="GENE")
datatablemouse16=data.table(mouse16,key="GENE")
datatablemouse17=data.table(mouse17,key="GENE")
datatablemouse18=data.table(mouse18,key="GENE")
datatablemouse19=data.table(mouse19,key="GENE")
datatablemouse20=data.table(mouse20,key="GENE")
datatablemouse21=data.table(mouse21,key="GENE")
datatablemouse22=data.table(mouse22,key="GENE")
datatablemouse23=data.table(mouse23,key="GENE")
datatablemouse24=data.table(mouse24,key="GENE")

### merge 24 ST batches
mouse=Reduce(function(x,y) merge(x,y,all=TRUE), list(datatablemouse1,datatablemouse2,datatablemouse3,datatablemouse4,datatablemouse5,datatablemouse6,datatablemouse7,datatablemouse8,datatablemouse9,datatablemouse10,datatablemouse11,datatablemouse12,datatablemouse13,datatablemouse14,datatablemouse15,datatablemouse16,datatablemouse17,datatablemouse18,datatablemouse19,datatablemouse20,datatablemouse21,datatablemouse22,datatablemouse23,datatablemouse24))

### clear data to save memory
rm(mouse1,mouse2,mouse3,mouse4,mouse5,mouse6,mouse7,mouse8,mouse9)
rm(list=ls(pattern="mouse1"))
rm(list=ls(pattern="mouse2"))
rm(list=ls(pattern="datatable"))
gc() # gc() release memory after remove the object

### convert back to data.frame, then change NA to 0
Sys.time()
mousetmp=as.data.frame(mouse)
mousetmp[is.na(mousetmp)] <- 0
Sys.time()                   # 14 min passed
### set gene names as row.names
row.names(mousetmp)=mousetmp[,1]
mousetmp[1:5,1:2]
mouse=mousetmp[,-1]
dim(mouse)             #[1] 36979 57854
### save the merged file
write.table(mouse, file=paste0(home,"data_DGE/MouseAdultST24/mergedMouseAdultST24_",dim(mouse)[2],"cells_",dim(mouse)[1],"genes.dge.txt"), quote=FALSE, sep='\t',row.names=T, col.names=T)

### Convert my data to sparse matrix
library(Matrix)
sparse1=as(as.matrix(mouse),"sparseMatrix")   
object.size(mouse)       #  17,124,594,432 bytes
object.size(sparse1)     #     914,601,664 bytes

### Write sparse matrix
writeMM(sparse1,paste0(home,"data_DGE/MouseAdultST24/matrix.mtx"))
write.table(colnames(sparse1), file=paste0(home,"data_DGE/MouseAdultST24/barcodes.tsv"), quote=FALSE, sep='\t', row.names=FALSE, col.names=FALSE)
write.table(cbind(rownames(sparse1),rownames(sparse1)), file=paste0(home,"data_DGE/MouseAdultST24/genes.tsv"), quote=FALSE, sep='\t',row.names=FALSE, col.names=FALSE)


###### load merged gene expression matrix for 24 ST batches
mouse <- Read10X(paste0(home,"data_DGE/MouseAdultST24"))
mouse[1:5,1:5]
object.size(mouse)     #   1,216,909,040 bytes

dgedata=mouse
name=name2="MouseAdultST24"
dgefile=dgename="figJun2017_MouseAdultST24/"
dataset=unique(gsub("_.*","",colnames(dgedata)))

###### Filter for cells (>500 detected genes and <10% MT)
nUMIperCell <- colSums(dgedata)
nGeneperCell <- colSums(dgedata>0)
dgedata.tmp=dgedata[,nGeneperCell>500]
mito.genes <- grep("^mt-", rownames(dgedata.tmp), value = T) # mouse
percent.mito <- colSums(dgedata.tmp[mito.genes, ])/colSums(dgedata.tmp)
length(which(percent.mito>0.1))       # 545
dgedata.tmp=dgedata[,nGeneperCell>500][,percent.mito<0.1]

### check nUMIs, nGenes, %MT, %x and %y after cell filter and before any gene filtering
dge <- new("seurat", raw.data = dgedata.tmp)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jun2017-MouseAdultST24",names.field = 1,names.delim = "_")
dge                # 36979 genes across 33180 samples.

# change order of dataset in orig.ident
dge@data.info$orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
dge@ident=factor(dge@ident,levels=dataset)

# add nUMIperCell to the dge datainfo
dge <- AddMetaData(dge, nUMIperCell, "nUMIperCell")

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes)
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

# add %ChrX and %ChrY genes to the dge datainfo
x=read.table(paste0(home,"data_DGE/mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(paste0(home,"data_DGE/mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570

chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
length(chrx.genes)  # 2045
length(chry.genes)  # 642
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")

# save dge for gene plot
dgeall=dge
save(dgeall, file =paste0(home,"data_DGE/MouseAdultST24genesall.Robj"))

# Average UMI, gene, %mito, %x, %y in each batch
c(mean(dge@data.info$nUMI),
mean(dge@data.info$nGene),
mean(dge@data.info$percent.mito),
mean(percent.x),
mean(percent.y) )

for(i in 1:length(dataset)){
  print(c(dataset[i],
mean(dge@data.info$nUMI[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])]),
mean(dge@data.info$nGene[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])])))
}

for(i in 1:length(dataset)){
  print(c(
mean(dge@data.info$percent.mito[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])]),
mean(percent.x[which(gsub("_.*","",names(percent.x))==dataset[i])]),
mean(percent.y[which(gsub("_.*","",names(percent.y))==dataset[i])])))
}

pdf(file=paste(dgefile,"NGeneUMImt.pdf",sep=""),height=4,width=11)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
dev.off()

###### To determine the criteria for gene filtering
nUMIperGene <- rowSums(dgedata.tmp)
nCellperGene <- rowSums(dgedata.tmp>0)

dim(dgedata.tmp)[2]*0.0006  # 19.9
dim(dgedata.tmp)[2]*0.00045 # 14.9

### double-check expression of lowly-expressed genes
perGene=cbind(nUMIperGene,nCellperGene)
perGene=perGene[order(perGene[,1],perGene[,2]),]
write.table(perGene,paste0(home,"data_DGE/MouseAdultST24_per",dim(perGene)[1],"Genesummary.txt"),quote=F,row.names=T,col.names=T,sep="\t")

pdf(file=paste(dgefile,"nUMICellperGene.pdf",sep=""),height=4,width=4)
plot(density(nUMIperGene), main="#UMIs per Gene")
plot(density(nCellperGene), main="#Non0 Cells per Gene")
plot(density(log2(nUMIperGene)), main="log2(#UMIs per Gene)")
plot(density(log2(nCellperGene)), main="log2(#Non0 Cells per Gene)")
plot(density(log10(nUMIperGene)), main="log10(#UMIs per Gene)")
plot(density(log10(nCellperGene)), main="log10(#Non0 Cells per Gene)")
dev.off()

##### rank all genes in each cell according to the number of normalized UMIs; then prioritize genes based on its abundance ranking; calculate mean and sd of gene ranking within each dataset and between datasets
### normalize and log-transform, require cells with >500 genes
dge <- new("seurat", raw.data = dgedata.tmp)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jun2017-MouseAdultST24",names.field = 1,names.delim = "_",do.scale=F,do.center=F)
dge                # 36979 genes across 33180 samples.
dge@data[1:5,1:5]

### Rank each gene for each cell based on the number of normalized UMIs
# give the most abundant gene rank = 1, the second most abundant rank = 2, and so on
# note for: rank(,ties.method)
# "average", "first", "last", "random", "max", "min"
length(rownames(dge@data))
# [1] 36979
rankexp=data.frame(matrix(0,dim(dge@data)[1],dim(dge@data)[2]))
names(rankexp)=colnames(dge@data)
rownames(rankexp)=rownames(dge@data)
for(i in 1:dim(dge@data)[2]){
cell=colnames(dge@data)[i]
rankexp[,cell]=rank(-dge@data[,cell],ties.method="average")
}
### Calculate mean and SD of gene rank for each dataset and for all datasets
rankexpmean=rankexpsd=data.frame(matrix(0,dim(dge@data)[1],length(unique(gsub("_.*","",colnames(dge@data))))+1))
names(rankexpmean)=names(rankexpsd)=c("Overall",unique(gsub("_.*","",colnames(dge@data))))
rownames(rankexpmean)=rownames(rankexpsd)= rownames(rankexp)
rankexpmean[,"Overall"]=apply(rankexp,1,mean)
rankexpsd[,"Overall"]=apply(rankexp,1,sd)
for(dataset in unique(gsub("_.*","",colnames(dge@data)))){
tmp=rankexp[,grep(dataset,names(rankexp))]
rankexpmean[,dataset]=apply(tmp,1,mean)
rankexpsd[,dataset]=apply(tmp,1,sd)
}
### Calculate SD of average rank and averaged SD among all datasets
rankexpsd=rankexpsd[,-("SD")]
rankexpsd$Mean=apply(rankexpsd[,-1],1,mean)
rankexpmean$SD=apply(rankexpmean[,-1],1,sd)
rankexpmean$OverallSD=rankexpsd$Overall
rankexpsd$OverallMean=rankexpmean$Overall
##### Prioritize top 50 highly-expressed (most abundant) genes (genes ranked top by overall mean abundance)
rankexpmeanorder=rankexpmean[order(rankexpmean$Overall),]
rankexpmeanorder[1:15,]
rankexpmeanorder[16:30,]
rankexpmeanorder[31:45,]
rankexpmeanorder[46:50,]
# check the SD for the top 50 genes
rankexpsd[rownames(rankexpmeanorder[1:30,] ),c(1,12)]
rankexpsd[rownames(rankexpmeanorder[31:50,] ),c(1,12)]
##### Prioritize top 50 most stably-expressed genes among those highly-expressed genes (genes with smallest overall SD for abundance)
# note: genes with 0 UMIs or very lowly-expressed would have the smallest overall SD, so I have to exclude genes that are very lowly-expressed first
### Plot Overall SD vs. overall Mean for each gene
length(which(rankexpsd$Overall>=1000))			# 20901
length(which(rankexpsd$OverallMean<=10000))		# 601
jpeg(file=paste(dgename,"generank.jpeg",sep=""),res=300,height=1600,width=2000)
plot(rankexpsd$OverallMean,rankexpsd$Overall,col=rgb(0,0,0,0.2),xlab="Average Rank of Each Gene",ylab="SD")
abline(v=10000,lty=2)
abline(h=1000,lty=2)
dev.off()


###### Filter for genes (Gene filter #1: genes with >20 UMIs/gene and >15 cells/gene, plus 31 additional genes of interest)
### keep genes with >20 UMIs/gene and >15 cells/gene
thres.nUMIperGene = 20
thres.nCellperGene = 15
length(which(nUMIperGene>thres.nUMIperGene)) # 24475
length(which(nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene))
gene1=rownames(dgedata.tmp)[nUMIperGene>thres.nUMIperGene & nCellperGene>thres.nCellperGene]
length(gene1)   # 24451

### add back 31 known genes of interest with <=20 UMIs/15cells per gene
gene2=c("Ngf","T","Tbx20","Tert","Olig2","Prl","Hmx2","Itgb6","Mageb1","Neurod1","Neurod2","Ptchd4","Prl6a1","Fgf17","Hpx","Npy4r","Fgf6","Foxf2","Prokr2","Sry","Ccl3","Prf1","Pax8","Dcx","Ly6g5c","Ereg","Pgr","Sycp2l","Prnd","Cysltr1","Pdf")
length(gene2) # 31

genekeep=c(gene1,gene2)
anyDuplicated(genekeep) # 0
length(genekeep) # 24482
# kept UMI>20,Cell>15,plus 31 genes of interest with UMI<=20

### save the filtered cell and filtered genes
dgedata2=dgedata.tmp[genekeep,]
dim(dgedata2)


###### Setup Seurat object with filtered cells and filtered genes for 24 ST batches (without INT6)
### setup object, normalize for cells, log-transform, standardize for genes
dge <- new("seurat", raw.data = dgedata2)
dge <- Setup(dge, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jun2017-MouseAdultST24",names.field = 1,names.delim = "_")
dge
# 24482 genes across 33180 cells.
table(gsub("_.*","",colnames(dgedata2)))
table(gsub("_.*","",colnames(dge@data)))

# change order of dataset in orig.ident
dataset=unique(gsub("_.*","",colnames(dge@data)))
dge@data.info$orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
dge@ident=factor(dge@ident,levels=dataset)

# add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes) # 27
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

# add %ChrX and %ChrY genes to the dge datainfo
x=read.table(paste0(home,"data_DGE/mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(paste0(home,"data_DGE/mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570

chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
length(chrx.genes)  # 752
length(chry.genes)  # 12
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")

# Calculate Average UMI, gene, %mito, %x, %y in each batch after gene filter
c(mean(dge@data.info$nUMI),
mean(dge@data.info$nGene),
mean(dge@data.info$percent.mito),
mean(percent.x),
mean(percent.y) )

for(i in 1:length(dataset)){
  print(c(dataset[i],
mean(dge@data.info$nUMI[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])]),
mean(dge@data.info$nGene[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])])))
}

for(i in 1:length(dataset)){
  print(c(
mean(dge@data.info$percent.mito[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])]),
mean(percent.x[which(gsub("_.*","",names(percent.x))==dataset[i])]),
mean(percent.y[which(gsub("_.*","",names(percent.y))==dataset[i])])))
}

pdf(file=paste(dgefile,"NGeneUMImt_filteredgenes20.pdf",sep=""),height=4,width=11)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
dev.off()

###### Selection of Highly Variable Genes (HVG)
pdf(paste(dgefile,"dge_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
length(dge@var.genes) # 2187

###### Perform PCA for Gene Filter #1
### PCA using all genes
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
Sys.time() # "2017-06-12 12:23:04 EDT"
dge <- PCA(dge, pc.genes = rownames(dge@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
Sys.time() # [1] "2017-06-13 19:05:51 EDT"
# flip PCs
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.rot[,2]=-dge@pca.rot[,2]
dge@pca.rot[,3]=-dge@pca.rot[,3]
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))
dge20=dge # Gene Filter #1, PCA using all genes

### PCA using highly-variable genes
dge <- PCA(dge, pc.genes=dge@var.genes, do.print = TRUE, pcs.print = 5, genes.print = 5)
# Project
dge=ProjectPCA(dge,do.print=FALSE)
# flip PCs
dge@pca.rot[,1]=-dge@pca.rot[,1]
dge@pca.rot[,2]=-dge@pca.rot[,2]
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15HVG.Robj"))
dge20HVG=dge # Gene Filter #1, PCA using HVG


###### Filter for genes (Gene Filter #2: genes with >1k UMIs) and do PCA
genekeep2=rownames(dgedata.tmp)[nUMIperGene>=1000]
dgedata22=dgedata.tmp[genekeep2,]
dge2 <- new("seurat", raw.data = dgedata22)
dge2 <- Setup(dge2, min.cells = 0, min.genes = 0, do.logNormalize = T, total.expr = 1e4, project = "Jun2017-MouseAdultST24",names.field = 1,names.delim = "_")
dge2
# change order of dataset in orig.ident
dataset=unique(gsub("_.*","",colnames(dge2@data)))
dge2@data.info$orig.ident=factor(dge2@data.info$orig.ident,levels=dataset)
dge2@ident=factor(dge2@ident,levels=dataset)
dge2 <- PCA(dge2, pc.genes = rownames(dge2@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
dge=dge2
dge@pca.rot[,1]=-dge@pca.rot[,1]
save(dge, file =paste0(home,"data_DGE/MouseAdultST24genesUMI1k.Robj"))
dge1k=dge # Gene Filter #2, PCA using all genes


###### Organize all gene filters
genefilter=c("UMI20cell15","UMI20cell15HVG","UMI1k")
dgelist=list(dge20,dge20HVG,dge1k)
dgefiles=c("figJun2017_MouseAdultST24/UMI20cell15_","figJun2017_MouseAdultST24/UMI20cell15HVG_","figJun2017_MouseAdultST24/UMI1k_")

### plot top PCs with different gene filters
dgefile="figJun2017_MouseAdultST24/"
plots=NULL
pdf(paste(dgefile,"dge_PCA_all2.pdf",sep=""),height=6,width=24)
for(pc in 1:20){
print(pc)
for(i in 1:3){
	dge=dgelist[[i]]
    plots[[i]]=PCAPlot(dge,1,pc,do.return = TRUE,pt.size = 1)
}
MultiPlotList(plots,cols = 3)
}
dev.off()

### analyze correlation of top PCs with different gene filters
compare=list(c(1,3),c(1,2),c(2,3))
pdf(paste(dgefile,"dge_PCA_cor.pdf",sep=""),height=4,width=12)
for(pc in 1:5){
par(mfrow=c(1,3),mar=c(4,4,2,1),mgp=c(2.5, 1, 0))
for(i in 1:3){
comp=compare[[i]]
plot(dgelist[[comp[1]]]@pca.rot[,pc],dgelist[[comp[2]]]@pca.rot[,pc],xlab=paste0("PC",pc," - ",genefilter[comp[1]]),col=rgb(0,0,0,0.1),pch=16,ylab=paste0("PC",pc," - ",genefilter[comp[2]]),cex.lab=1.5)
abline(lm(dgelist[[comp[2]]]@pca.rot[,pc]~dgelist[[comp[1]]]@pca.rot[,pc]),col="green3")
mtext(paste0("cor=",round(cor(dgelist[[comp[1]]]@pca.rot[,pc],dgelist[[comp[2]]]@pca.rot[,pc]),2),",p=",cor.test(dgelist[[comp[1]]]@pca.rot[,pc],dgelist[[comp[2]]]@pca.rot[,pc])$p.val),3,col="green3")
}
}
dev.off()

###### Determining the number of top PCs using Scree Plot and density plot of eigenvalues
### Scree plot
numPCs=c(40,26,40)
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""),height=4,width=12)
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in 1:3){
dge=dgelist[[i]]
plot(dge@pca.obj[[1]]$sdev[1:120],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,3,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=genefilter[i],cex=1.5)
}
dev.off()

### Density plot of Eigenvalues
pdf(paste(dgefile,"dge_PCA_Variablel_eigenvalue.pdf",sep=""),height=4,width=12)
par(mfrow=c(1,3),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
for(i in 1:3){
dge=dgelist[[i]]
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
print(eigenvalue)
plot(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue",main="",cex.lab=1.5)
polygon(density(dge@pca.obj[[1]]$sdev),col="black")
lines(density(dge@pca.obj[[1]]$sdev),col="red",lwd=2,xlab="Eigenvalue")
abline(v=eigenvalue,col="red",lwd=2,lty=2)
text(eigenvalue+1.2,0.8,col="red",paste(numPCs[i],"PCs"))
legend("topright",legend=genefilter[i],cex=1.5)
}
dev.off() 


###### tSNE using top PCs
for(i in 1:3){
print(i)
dge=dgelist[[i]]
dge=RunTSNE(dge,dims.use = 1:numPCs[i],do.fast=T)    # max_iter=2000
dgelist[[i]]=dge
}

### Application of Louvain-Jaccard Clustering
i=1
dge=dgelist[[i]]
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.6,4,0.1), print.output = 0, save.SNN = T)
dgelist[[i]]=dge
dge20=dge

i=2
dge=dgelist[[i]]
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.6,4,0.1), print.output = 0, save.SNN = T)
dgelist[[i]]=dge
dge20HVG=dge

i=3
dge=dgelist[[i]]
dge <- FindClusters(dge, pc.use = 1:numPCs[i], resolution = seq(0.6,4,0.1), print.output = 0, save.SNN = T)
dgelist[[i]]=dge
dge1k=dge

### tSNE plot
plott=list()
for(i in 1:3){
print(i)
dge=dgelist[[i]]
#dge <- SetAllIdent(dge, id = "orig.ident")
plott[[i]]=TSNEPlot(dge,do.return=TRUE,pt.size = 1)
}
pdf(paste0(dgefile,"tSNE.pdf"),width=24,height=6)
MultiPlotList(plott,cols = 3)
dev.off()

###### cluster attributes from Louvain-Jaccard Clustering
for(i in 1:3){
print(i)
dge=dgelist[[i]]
print(c( length(unique(dge@data.info$res.0.6)),length(unique(dge@data.info$res.0.7)),length(unique(dge@data.info$res.0.8)),length(unique(dge@data.info$res.0.9)),length(unique(dge@data.info$res.1)),length(unique(dge@data.info$res.1.1)),length(unique(dge@data.info$res.1.2)),length(unique(dge@data.info$res.1.3)),length(unique(dge@data.info$res.1.4)),length(unique(dge@data.info$res.1.5)),length(unique(dge@data.info$res.1.6)),length(unique(dge@data.info$res.1.7)),length(unique(dge@data.info$res.1.8)),length(unique(dge@data.info$res.1.9)),length(unique(dge@data.info$res.2)),length(unique(dge@data.info$res.2.1)),length(unique(dge@data.info$res.2.2)),length(unique(dge@data.info$res.2.3)),length(unique(dge@data.info$res.2.4)),length(unique(dge@data.info$res.2.5)),length(unique(dge@data.info$res.2.6)),length(unique(dge@data.info$res.2.7)),length(unique(dge@data.info$res.2.8)),length(unique(dge@data.info$res.2.9)),length(unique(dge@data.info$res.3)),length(unique(dge@data.info$res.3.1)),length(unique(dge@data.info$res.3.2)),length(unique(dge@data.info$res.3.3)),length(unique(dge@data.info$res.3.4)),length(unique(dge@data.info$res.3.5)),length(unique(dge@data.info$res.3.6)),length(unique(dge@data.info$res.3.7)),length(unique(dge@data.info$res.3.8)),length(unique(dge@data.info$res.3.9)),length(unique(dge@data.info$res.4)) ))
}

### plot PCA and tSNE for the same number of clusters for each gene filter
res=paste0("res.",c(0.8,1.7,1.1))  # 31 clusters
n=nclusters=31

dgefile="figJun2017_MouseAdultST24/"
plot2=plott=NULL
### re-order clusters by hierarchical clustering
pdf(paste0(dgefile,"cluster_",n,"_tree.pdf"),width=12,height=3)
par(mfrow=c(1,3),mar=c(0.5,0.5,0.5,0.5),mgp=c(2.5, 1, 0))
for(i in 1:3){
dge=dgelist[[i]]
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])
#PlotClusterTree(dge)
dge@data.info[,res[i]]=dge@ident
dgelist[[i]]=dge
plot2[[i]]=PCAPlot(dge,1,2,do.return = TRUE,pt.size = 1,do.label=T)
plott[[i]]=TSNEPlot(dge,do.return=TRUE,pt.size = 1,do.label=T,label.size=8)
dev.off()
pdf(paste0(dgefile,"cluster_",n,".pdf"),width=24,height=6)
MultiPlotList(plot2,cols = 3)
MultiPlotList(plott,cols = 3)
dev.off()

### plot more PC dimensions for 31 clusters
dgefile="figJun2017_MouseAdultST24/"
plots=NULL
pdf(paste0(dgefile,"cluster_",n,"_morePCs.pdf"),width=24,height=6)
for(pc in 3:10){
print(pc)
for(i in 1:3){
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
    plots[[i]]=PCAPlot(dge,1,pc,do.return = TRUE,pt.size = 1,do.label=T)
}
MultiPlotList(plots,cols = 3)
}
dev.off()

### Cross tabulation of cluster IDs
cc=100
for(j in 1:length(nclusters)){
n=nclusters[j]
res=resall[[j]]
table=NULL
for(i in 1:3){
dge=dgelist[[i]]
dge <- SetAllIdent(dge, id = res[i])
table[[i]]=dge@ident
}
comp<-table(table[[1]],table[[2]])
comp<-table(table[[1]],table[[3]])
new=NULL
rowid=NULL
colid=as.numeric(colnames(comp))
rowc=1
for(j in 1:dim(comp)[2]){
### For each row, sort numbers and select numbers >cc in the column, then fix these rows and move to next row and next column
maxc=comp[order(comp[,j],decreasing=T),]
rowc=max(1,length(which(maxc[,j]>cc)))
rowid=c(rowid,rownames(maxc[1:(rowc+1),])[1:rowc])
comp=maxc[-c(1:rowc),]
new=rbind(new,maxc[1:rowc,])
if(is.vector(comp)){
break
}
}
new=rbind(new,comp)
rowid=c(rowid,rownames(maxc)[-1])
rownames(new)=rowid
print(new)
}

### save dge
for(i in 1:3){
dge=dgelist[[i]]
save(dge,file =paste0(home,"data_DGE/MouseAdultST24genes",genefilter[i],".Robj"))
}

### Contribution of each batch to each cluster
for(i in 1:3){
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  print(ncellsbatchcluster)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  print(percentcellsbatch)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  print(percentcellscluster)
}
for(i in 1:3){
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefiles[i],"ncellspercluster_batch.txt"),quote=F,row.names=T,col.names=T,sep="\t")
}

### Save 31 clusters from gene filter #1 with HVG and from gene filter #2 in the dge@data.info of gene filter #1 with all genes
# dge=dgelist[[1]] is furthur used in marker gene analysis
i=1
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
  tmp=dge@ident
  i=1
  dge=dgelist[[i]]
  dge=AddMetaData(dge,tmp,"clusters31_GeneFilter1")
  dgelist[[1]]=dge
i=2
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
  tmp=dge@ident
  i=1
  dge=dgelist[[i]]
  dge=AddMetaData(dge,tmp,"clusters31_GeneFilter1HVG")
  dgelist[[1]]=dge
i=3
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])
  tmp=dge@ident
  i=1
  dge=dgelist[[i]]
  dge=AddMetaData(dge,tmp,"clusters31_GeneFilter2")
  dgelist[[1]]=dge

dge@data.info[1:5,-(1:41)]
table(dge@data.info$clusters31_GeneFilter1,dge@data.info$clusters31_GeneFilter1HVG)
save(dge,file =paste0(home,"data_DGE/MouseAdultST24genes",genefilter[1],".Robj"))
dge20=dge

### Comparison with 14 clusters of 6 ST
ident14=read.table(paste0(home,"data_DGE/ident_14clusters_ST6.txt"))
ident2=ident14[,1]
names(ident2)=rownames(ident14)
names(ident2)=gsub("T-7wk-A","ST1",names(ident2))
names(ident2)=gsub("T-7wk-B","ST2",names(ident2))
names(ident2)=gsub("T-9wk-1","ST3",names(ident2))
names(ident2)=gsub("T-9wk-2","ST4",names(ident2))
names(ident2)=gsub("T-9wk-A","ST5",names(ident2))
names(ident2)=gsub("T-9wk-B","ST6",names(ident2))
comp=table(dge@ident[names(ident2)],ident2)
sum(comp) # 12991
print(comp)

### Distribution of Number of UMIs and %MT in Each Cluster
i=1
res=paste0("res.",c(0.8,1.7,1.1))  # 31 clusters
  dge=dgelist[[i]]
  dge <- SetAllIdent(dge, id = res[i])

pdf(file=paste(dgefile,"NGeneUMImt_31clusters.pdf",sep=""),height=4,width=11)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
dev.off()

ncluster=31
## nUMI per Cell
for(j in 1:ncluster){
	print(summary(dge@data.info$nUMIperCell2[which(dge@ident==j)]))
}
## %MT per cell
for(j in 1:ncluster){
	print(summary(dge@data.info$percent.mito[which(dge@ident==j)]))
}


###### Order 31 clusters from Gene Filter #1 with all genes
### order cell by cluster ID 
dge=dgelist[[1]]
dge=SetAllIdent(dge,id="res.0.8")
ident=dge@ident
levels=levels(ident)
cells=sort(ident)
cells.use=NULL
for(i in 1:length(levels)){
   tmpname=names(cells[which(cells == levels[i])])
   cells.use=c(cells.use,tmpname)
}
cells.ident=as.factor(cells)
names(cells.ident)=cells.use
levels(cells.ident)=levels

### Calculate centroid of each cluster
tmpdge=data.frame(t(as.matrix(dgeall@data[,cells.use])))
## make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dgeall@data)[cells.use])  # nothing
mouseclustersall=data.frame(ident=cells.ident,t(as.matrix(dgeall@data[,cells.use])))
### for each cluster, calculate average normalized expression of each gene
genecountsall=matrix(,dim(dgeall@data)[1],length(unique(mouseclustersall$ident)))
rownames(genecountsall)=rownames(dgeall@data)
colnames(genecountsall)=unique(mouseclustersall$ident)
for(i in unique(mouseclustersall$ident)){
    genecountsall[,i]=apply(mouseclustersall[mouseclustersall$ident==i,-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Ordering the 31 cluster centroids using Seriation by OLO
library(seriation)
n=ncluster=31
tmp=genecountsall[,levels][,1:ncluster]
da <- dist(t(as.matrix(tmp)), method = "euclidean")
pdf(file=paste(dgename,"Centroid_norm_Seriation_Dissimilarity.pdf"))
# plot original matrix
dissplot(da, method = NA,options = list(main = paste("24ST Dissimilarity")))
# plot with seriation using method “ARSA”
dissplot(da, method="ARSA", options = list(main = paste("24ST Dissimilarity with seriation ARSA")))
# plot with seriation using method “OLO”
dissplot(da, method="OLO",options = list(main = paste("24ST Dissimilarity with seriation OLO")))
dev.off()

### obtain ordered 31 cluster IDs
levelss=get_order(seriate(da,method="OLO"))
levels=levelss
levelss #8 10 7 9 6 5 4 2 1 3 11 13 15 16 19 18 17 12 14 22 23 24 21 30 20 31 29 28 27 26 25

### Calculate ranked correlation for each normalized cluster centroid
cc=cor(as.matrix(genecountsall),method="spearman")
data.use=cc[levels,levels]
## add separation line between each cluster and add labeling of cluster IDs
colsep.use=cumsum(table(gsub("_.*","",levels)))
col.lab=rep("",length(levels))
col.lab[round(cumsum(table(gsub("_.*","",levels)))-table(gsub("_.*","",levels))/2)+0]=unique(gsub("_.*","",levels))
row.lab=gsub(".*_","",levels)
## add color bar
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=rep(myBrewerPalette[1:31],5)[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
colnames(clab)=rownames(rlab)=c("","Cluster")
## color scheme
col.use=redblue100
### plot correlation of ordered 31 cluster centroid
col.use=redblue100
pdf(file=paste(dgename,"Centroid_RankedCorrelation_indiv_31clusters_seriation.pdf",sep=""),height=5.5,width=5)
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,rowsep=colsep.use,sepcolor="black", labCol=col.lab,labRow=row.lab, symm=F,symkey=F,symbreaks=F, scale="none")
dev.off()

### reorder cells using ordered 31 cluster IDs and shuffle cells randomly within each cluster
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
            data.use2=testcor[cells.use,cells.use]
            data.use2=minmax(data.use2,min=disp.min,max=disp.max)
#3 add labels of batch names for row or column 
            lab2=rep("",length(cells.use))
            lab2[round(cumsum(table(cells.ident)[levels(cells.ident)])-table(cells.ident)[levels(cells.ident)]/2)+15]=levels(cells.ident)
            row.lab2=gsub(".*_","",lab2)
            orig.ident=factor(gsub("_.*","",cells.ident),levels=unique(gsub("_.*","",cells.ident)))
            col.lab2=rep("",length(cells.use))
            col.lab2[round(cumsum(table(orig.ident)[levels(orig.ident)])-table(orig.ident)[levels(orig.ident)]/2)+15]=levels(orig.ident)
## add separation line between each dataset
            colsep.use2=cumsum(table(cells.ident)[levels(cells.ident)])
            colsep.use2=cumsum(table(orig.ident)[levels(orig.ident)]) 
## add color bar
sidecol2=do.call(rbind,strsplit(as.character(cells.ident),"_"))
sidecol2=cbind(sidecol2,sidecol2)
for(rep in 1:length(unique(sidecol2[,1]))){
a=unique(sidecol2[,1])[rep]
sidecol2[which(sidecol2[,1]==a),2]<-rep
}
ncluster=31
clusters=1:ncluster
nbatch=1
for(i in 1:nbatch){
start=cumsum(c(0,ncluster[-nbatch]))[i]
for(j in 1:ncluster[i]){
sidecol2[which(cells.ident==unique(cells.ident)[start+j]),2]<-j
}
}
rlab2=rbind(c("white")[as.numeric(sidecol2[,1])],myBrewerPalette[as.numeric(sidecol2[,2])])
clab2=cbind(rlab2[2,],rlab2[1,])
colnames(clab2)=rownames(rlab2)=c("","Cluster")
## color scheme
midrange2=midrange
maxrange2=midrange2*2
col.use2=redblue100[c(rep(1:50,each=round(midrange2*100)),rep(50:100,each=round((maxrange2-midrange2)*100)),rep(100,50*round((1-maxrange2)*100)))]
length(col.use2)
### plot correlation of all cells in ordered 31 clusters
jpeg(file=paste(dgename,"RankedCorrelation_indiv_31clusters_3k_seriation.jpeg",sep=""),height=3000,width=3000,res=300)
heatmap.3(data.use2,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use2,colsep = colsep.use2,rowsep=colsep.use2,sepcolor="black", labCol=col.lab2,labRow=row.lab2, symm=F,symkey=F,symbreaks=F,scale="none")                    
dev.off()

### save the ordered 31 cluster IDs
which(unique(sidecol2[,1])!=get_order(do)) # integer(0)
cells.ident.ordered=factor(as.numeric(sidecol2[,2]),ordered=TRUE)
names(cells.ident.ordered)=names(cells.ident)

ordered="Cluster31Seriation_GeneFilter1"
dge=dgelist[[1]]
dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge@data.info[,ordered]=factor(dge@data.info[,ordered])
dge=SetAllIdent(dge,ordered)
# save the dge file
dgelist[[i]]=dge
dge20=dge
save(dge,file = paste0(home,"data_DGE/MouseAdultST24genes",genefilter[1],".Robj"))
which(dgelist[[i]]@ident != dgelist[[i]]@data.info$ClusterSeriation)

### Visualize the 31 ordered clusters in PCA and tSNE
pdf(paste(dgefile,"Cluster31Seriation_GeneFilter1.pdf"),width=9,height=7.5)
PCAPlot(dge,do.label=TRUE)
TSNEPlot(dge,do.label=TRUE,label.size=8)
dev.off()


###### Finding differentially expressed genes for all clusters
library(Seurat)
library(FNN)
library(igraph)
library(dplyr)

home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
genefilter=c("UMI20cell15","UMI20cell15HVG","UMI1k")
dgefiles=c("figJun2017_MouseAdultST24/UMI20cell15_","figJun2017_MouseAdultST24/UMI20cell15HVG_","figJun2017_MouseAdultST24/UMI1k_")
i=1
load(file =paste0(home,"data_DGE/MouseAdultST24genes",genefilter[1],".Robj"))
dge=dgelist[[1]]
dge=SetAllIdent(dge,id="res.0.8") # 31 clusters
# the order is: 8 10 7 9 6 5 4 2 1 3 11 13 15 16 19 18 17 12 14 22 23 24 21 30 20 31 29 28 27 26 25

markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod",do.print = TRUE)
table(markersall$cluster)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_",genefilter[i],"_pct0.2_diffpct0.2_thresh2fold_clustersPC1-40_6.19.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")

# 6.26.2017 check markers in neighboring clusters
markers29=FindMarkers(dge,29,c(20,28),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers2=FindMarkers(dge,2,c(1,3),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers11=FindMarkers(dge,11,c(4,5,13),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers12=FindMarkers(dge,12,c(17,14),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers14=FindMarkers(dge,14,c(17,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers20=FindMarkers(dge,20,29,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers28=FindMarkers(dge,28,29,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers202=FindMarkers(dge,20,c(21,29),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers282=FindMarkers(dge,28,c(29,25),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markers29$cluster="29Vs20&28";markers29$gene=rownames(markers29)
markers2$cluster="2Vs1&3";markers2$gene=rownames(markers2)
markers11$cluster="11Vs4&5&13";markers11$gene=rownames(markers11)
markers12$cluster="12Vs17&14";markers12$gene=rownames(markers12)
markers14$cluster="14Vs17&22";markers14$gene=rownames(markers14)
markers20$cluster="20Vs29";markers20$gene=rownames(markers20)
markers28$cluster="28Vs29";markers28$gene=rownames(markers28)
markers202$cluster="20Vs21&29";markers202$gene=rownames(markers202)
markers282$cluster="28Vs29&25";markers282$gene=rownames(markers282)

markersneighbor=rbind(markers2,markers11,markers12,markers14,markers29)
markersneighbor=rbind(markers20,markers202,markers28,markers282)

table(markersneighbor$cluster)
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)

# 6.30.2017
markers2=FindMarkers(dge,2,c(1,3),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers22=FindMarkers(dge,2,c(1,3,11),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers1=FindMarkers(dge,1,c(2,3),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers12=FindMarkers(dge,1,c(2,3,11),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers3=FindMarkers(dge,3,c(1,11),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers32=FindMarkers(dge,3,c(1,2,11,13),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers11=FindMarkers(dge,11,c(3,13),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers112=FindMarkers(dge,11,c(1,2,3,13,15),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers113=FindMarkers(dge,11,c(1,2,3,15),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers13=FindMarkers(dge,13,c(11,15),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers132=FindMarkers(dge,13,c(1,2,3,11,15),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers133=FindMarkers(dge,13,c(1,2,3,15),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markers2$cluster="2Vs1&3";markers2$gene=rownames(markers2)
markers22$cluster="2Vs1&3&11";markers22$gene=rownames(markers22)
markers1$cluster="1Vs2&3";markers1$gene=rownames(markers1)
markers12$cluster="1Vs2&3&11";markers12$gene=rownames(markers12)
markers3$cluster="3Vs1&11";markers3$gene=rownames(markers3)
markers32$cluster="3Vs1&2&11&13";markers32$gene=rownames(markers32)
markers11$cluster="11Vs3&13";markers11$gene=rownames(markers11)
markers112$cluster="11Vs1&2&3&13&15";markers112$gene=rownames(markers112)
markers113$cluster="11Vs1&2&3&15";markers113$gene=rownames(markers113)
markers13$cluster="13Vs11&15";markers13$gene=rownames(markers13)
markers132$cluster="13Vs1&2&3&11&15";markers132$gene=rownames(markers132)
markers133$cluster="13Vs1&2&3&15";markers133$gene=rownames(markers133)

markersneighbor=rbind(markers2,markers22,markers1,markers12,markers3,markers32,markers11,markers112,markers113,markers13,markers132,markers133)
print(cbind(table(markersneighbor$cluster)[unique(markersneighbor$cluster)],
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)[unique(markersneighbor$cluster)],
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)[unique(markersneighbor$cluster)] ))

markers12=FindMarkers(dge,12,c(17,14),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers122=FindMarkers(dge,12,c(17,14,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers123=FindMarkers(dge,12,c(13,17,14,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers14=FindMarkers(dge,14,c(12,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers141=FindMarkers(dge,14,c(17,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers142=FindMarkers(dge,14,c(17,12,22),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers143=FindMarkers(dge,14,c(17,12,22,23),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers22=FindMarkers(dge,22,c(14,23),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers222=FindMarkers(dge,22,c(12,14,23,24),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers223=FindMarkers(dge,22,c(14,21),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers23=FindMarkers(dge,23,c(22,24),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers232=FindMarkers(dge,23,c(14,22,24,21),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers233=FindMarkers(dge,23,c(14,21),min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markers12$cluster="12Vs17&14";markers12$gene=rownames(markers12)
markers122$cluster="12Vs17&14&22";markers122$gene=rownames(markers122)
markers123$cluster="12Vs13&17&14&22";markers123$gene=rownames(markers123)
markers14$cluster="14Vs12&22";markers14$gene=rownames(markers14)
markers141$cluster="14Vs17&22";markers141$gene=rownames(markers141)
markers142$cluster="14Vs17&12&22";markers142$gene=rownames(markers142)
markers143$cluster="14Vs17&12&22&23";markers143$gene=rownames(markers143)
markers22$cluster="22Vs14&23";markers22$gene=rownames(markers22)
markers222$cluster="22Vs12&14&23&24";markers222$gene=rownames(markers222)
markers223$cluster="22Vs14&21";markers223$gene=rownames(markers223)
markers23$cluster="23Vs22&24";markers23$gene=rownames(markers23)
markers232$cluster="23Vs14&22&24&21";markers232$gene=rownames(markers232)
markers233$cluster="23Vs14&21";markers233$gene=rownames(markers233)

markersneighbor=rbind(markers12,markers122,markers123,markers14,markers141,markers142,markers143,markers22,markers222,markers223,markers23,markers232,markers233)
print(cbind(table(markersneighbor$cluster)[unique(markersneighbor$cluster)],
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)[unique(markersneighbor$cluster)],
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)[unique(markersneighbor$cluster)] ))


######### merging 31 ordered clusters to represent major cell types
dge=dge20
dge=SetAllIdent(dge,id="res.0.8") # 31 clusters
# the order for the 31 clusters is: 8 10 7 9 6 5 4 2 1 3 11 13 15 16 19 18 17 12 14 22 23 24 21 30 20 31 29 28 27 26 25

somatic=c(8,10,7,9,6,5,4)
spg=c(2,1,3,11,13) # the cluster 2-1-3 was SPG, while cluster 11 & 13 were transitions between SPG and Scytes
spermatocyte=c(15,16,19,18,17,12,14)
round=c(22,23,24,21,30,20,31,29)
elongating=c(28,27,26,25)

dge@data.info$celltype2=rep(NA,length(dge@data.info$res.0.8))
levels(dge@data.info$celltype2)=c("Somatic","Spermatogonia","Spermatocyte","RoundSpermatid","Round/Elongating")
dge@data.info$celltype2[which(dge@data.info$res.0.8 %in% somatic)] <- "Somatic"
dge@data.info$celltype2[which(dge@data.info$res.0.8 %in% spg)] <- "Spermatogonia"
dge@data.info$celltype2[which(dge@data.info$res.0.8 %in% spermatocyte)] <- "Spermatocyte"
dge@data.info$celltype2[which(dge@data.info$res.0.8 %in% round)] <- "RoundSpermatid"
dge@data.info$celltype2[which(dge@data.info$res.0.8 %in% elongating)] <- "Round/Elongating"

table(dge@data.info$celltype2)[levels(dge@data.info$celltype2)]
### these include a somatic cell group and 4 major germ cell types, the latter of which are colored in Fig1B

save(dge, file = paste0(home,"data_DGE/MouseAdultST24genes",genefilter[1],".Robj"))
dge20=dge


