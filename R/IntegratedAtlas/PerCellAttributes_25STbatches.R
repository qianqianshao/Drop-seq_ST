### R script to calculate per-cell attributes for 35k cells of merged 25 ST batches in Dec 2017 by Qianyi
### Related to Figure 1D, Figure S1F-G, Table S1 and PerCellAttributes text file in GEO: GSE112393

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(gplots)

############ 25 Adult Mouse ST batches
### path
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
dgefile=dgename="figDec2017_MouseAdultST25/All_"
### Load filtered cells and all detected genes for 25 ST batches
load(file =paste0(home,"data_DGE/MouseAdultST25genesall.Robj"))
dgeall # 37241 genes across 34633 samples.
dge=dgeall

######### Per-cell attributes for nGene, nUMI, %MT, %ChrX, %ChrY and %Autosome - Figure 1D left panels
###### add percent.mito to the dge datainfo
mito.genes <- grep("^mt-", rownames(dge@data), value = T)
length(mito.genes) # 32 dgeall
percent.mito <- colSums(expm1(dge@data[mito.genes, ]))/colSums(expm1(dge@data))
dge <- AddMetaData(dge, percent.mito, "percent.mito")

###### add %ChrX and %ChrY genes to the metadata of 25 ST datasets
### load genes from each Chr
x=read.table(paste0(home,"data_DGE/mouseChrXgenes"),stringsAsFactors=FALSE)[,1]
y=read.table(paste0(home,"data_DGE/mouseChrYgenes"),stringsAsFactors=FALSE)[,1]
autosome=read.table(paste0(home,"data_DGE/mouseautosomegenes"),stringsAsFactors=FALSE)
chr1=read.table(paste0(home,"data_DGE/mouseChr1genes"),stringsAsFactors=FALSE)[,1]
chr3=read.table(paste0(home,"data_DGE/mouseChr3genes"),stringsAsFactors=FALSE)[,1]
chr8=read.table(paste0(home,"data_DGE/mouseChr8genes"),stringsAsFactors=FALSE)[,1]
chr10=read.table(paste0(home,"data_DGE/mouseChr10genes"),stringsAsFactors=FALSE)[,1]
length(x) # 2605
length(y) # 1570
dim(autosome) # 41390
table(autosome[,2])
length(chr1) # 3419
length(chr3) # 2874
length(chr8) # 1782
length(chr10) # 1619
### keep genes from each Chr that were detected in 25 ST datasets
chrx.genes <- x[which(x %in% rownames(dge@data))]
chry.genes <- y[which(y %in% rownames(dge@data))]
autosome.genes <- autosome[which(autosome[,1] %in% rownames(dge@data)),]
chr1.genes <- chr1[which(chr1 %in% rownames(dge@data))]
chr3.genes <- chr3[which(chr3 %in% rownames(dge@data))]
chr8.genes <- chr8[which(chr8 %in% rownames(dge@data))]
chr10.genes <- chr10[which(chr10 %in% rownames(dge@data))]
length(chrx.genes)  # 2062
length(chry.genes)  # 643
dim(autosome.genes) # 34502
table(autosome.genes[,2])
length(chr1.genes)  # 2995
length(chr3.genes)  # 2483
length(chr8.genes)  # 1488
length(chr10.genes) # 1343
### calculate the proportion of each Chr expression
percent.x <- colSums(expm1(dge@data[chrx.genes, ]))/colSums(expm1(dge@data))
percent.y <- colSums(expm1(dge@data[chry.genes, ]))/colSums(expm1(dge@data))
percent.autosome <- colSums(expm1(dge@data[autosome.genes[,1], ]))/colSums(expm1(dge@data))
percent.1 <- colSums(expm1(dge@data[chr1.genes, ]))/colSums(expm1(dge@data))
percent.3 <- colSums(expm1(dge@data[chr3.genes, ]))/colSums(expm1(dge@data))
percent.8 <- colSums(expm1(dge@data[chr8.genes, ]))/colSums(expm1(dge@data))
percent.10 <- colSums(expm1(dge@data[chr10.genes, ]))/colSums(expm1(dge@data))
### add to metadata
dge <- AddMetaData(dge, percent.x, "percent.x")
dge <- AddMetaData(dge, percent.y, "percent.y")
dge <- AddMetaData(dge, percent.autosome, "percent.autosome")
dge <- AddMetaData(dge, percent.1, "percent.Chr1")
dge <- AddMetaData(dge, percent.3, "percent.Chr3")
dge <- AddMetaData(dge, percent.8, "percent.Chr8")
dge <- AddMetaData(dge, percent.10, "percent.Chr10")
dgeall=dge
write.table(dge@data.info,file =paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),quote=F,row.names=T,col.names=T,sep="\t")

###### Per-Dataset Attributes
# change order of dataset in orig.ident
dataset=unique(gsub("_.*","",colnames(dge@data)))[c(1:16,25,17:24)]
dge@data.info$orig.ident=factor(dge@data.info$orig.ident,levels=dataset)
dge@ident=factor(dge@ident,levels=dataset)
# Calculate Average UMI, gene, %mito, %x, %y in each batch 
c(mean(dge@data.info$nUMI),
mean(dge@data.info$nGene),
mean(dge@data.info$percent.mito),
mean(percent.x),
mean(percent.y) )

dataset=as.character(unique(dge@data.info$orig.ident))
for(i in 1:length(dataset)){
  print(c(dataset[i],
mean(dge@data.info$nUMI[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])]),
mean(dge@data.info$nGene[which(gsub("_.*","",rownames(dge@data.info))==dataset[i])])))
}
### saved in Table S1 (two columns on the right)


###### Distribution of nGenes, nUMIs, %MT, %ChrX and %ChY in 11 major cell types - Figure 1D left panels
### 11 major cell types
dge=SetAllIdent(dge,id="celltype")
dge@ident=factor(dge@ident,levels=levels(dge@data.info$celltype))
### color scheme for 11 major cell types 
library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7,1:4)] # started using this on 12/20/2017
### Violin plot for per-cell attributes
pdf(file=paste(dgefile,"NGeneUMI_11clusters.pdf",sep=""),height=4,width=8)
VlnPlot(dge, "nGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nUMI", nCol = 1,cols.use=myBrewerPalette,y.log=T)
VlnPlot(dge, "nReads", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nReads", nCol = 1,cols.use=myBrewerPalette,y.log=T)
dev.off()
pdf(file=paste(dgefile,"percentMTChr_11clusters.pdf",sep=""),height=4,width=8)
VlnPlot(dge, "percent.mito", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.x", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.y", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.autosome", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.Chr1", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.Chr3", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.Chr8", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "percent.Chr10", nCol = 1,cols.use=myBrewerPalette)
dev.off()
### save as Figure 1D left panels


######### select highly variable genes for merged 25 ST datasets
pdf(paste(dgefile,"dge_VariableGenes0.2.pdf",sep=""),height=7.5,width=11)
dge=MeanVarPlot(dge,y.cutoff = 0.2,x.low.cutoff = 0.2,x.high.cutoff=10,y.high.cutoff=30,fxn.x = expMean,fxn.y = logVarDivMean,do.text=FALSE)    # x-axis: average expression; y-axis: dispersion, SD
points(data.x[pass.cutoff],data.norm.y[pass.cutoff],col="green3",pch=20,cex=0.6)
legend("topright",pch=20,cex=1.5,col="green3",legend=paste(length(dge@var.genes),"Variable Genes"))
dev.off()
length(dge@var.genes) # 2521 HVG y.cutoff=0.2 for dgeall, used this
dgeall=dge


######### Gini index for each cell across 11 major cell types - Figure 1D right panel and Figure S1F
###### Calculate Lorenz Curve and Gini index for each cell
# Gini coefficient of nUMI for each cell: represents uneven distribution of genes for each cell
library(ineq)
### calculate Gini using all ~37k genes
GiniAll=apply( dge@raw.data,2,function(x) ineq(x,type=c("Gini")) )
length(GiniAll)  # [1] 34633
summary(GiniAll)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9091  0.9573  0.9735  0.9706  0.9847  0.9953
write.table(GiniAll,file =paste0(home,"data_DGE/MouseAdultST25_Gini_genesall.txt"),row.names=T,col.names=F,sep="\t",quote=F)
### calculate Gini using detected genes for each cell 
GiniNon0=apply( dge@raw.data,2,function(x) ineq(x[which(x!=0)],type=c("Gini")) )
length(GiniNon0) # [1] 34633
summary(GiniNon0)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.09293 0.29420 0.41300 0.40290 0.51030 0.75250
write.table(GiniNon0,file =paste0(home,"data_DGE/MouseAdultST25_Gini_genesNon0.txt"),row.names=T,col.names=F,sep="\t",quote=F)
### calculate Gini using highly-variable genes
length(dge@var.genes) # 2521
GiniHVG=apply( dge@raw.data[dge@var.genes,],2,function(x) ineq(x,type=c("Gini")) )
length(GiniHVG)   # [1] 34633
summary(GiniHVG)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.9055  0.9546  0.9715  0.9683  0.9828  0.9944 
write.table(Gini1k,file =paste0(home,"data_DGE/MouseAdultST25_Gini_genesHVG.txt"),row.names=T,col.names=F,sep="\t",quote=F)

### Add Gini index to metadata
GiniAll=t(read.table(paste0(home,"data_DGE/MouseAdultST25_Gini_genesall.txt"),header=F,row.names=1))[1,]
GiniNon0=t(read.table(paste0(home,"data_DGE/MouseAdultST25_Gini_genesNon0.txt"),header=F,row.names=1))[1,]
GiniHVG=t(read.table(paste0(home,"data_DGE/MouseAdultST25_Gini_genesHVG.txt"),header=F,row.names=1))[1,]
dge=AddMetaData(dge,GiniAll,"GiniAll")
dge=AddMetaData(dge,GiniNon0,"GiniNon0")
dge=AddMetaData(dge,GiniHVG,"GiniHVG")
dgeall=dge
write.table(dge@data.info,file =paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),quote=F,row.names=T,col.names=T,sep="\t")

### plot example Lorenz curve to double-check if Gini calculation makes sense
GiniUMI=GiniAll
## for cells with high Gini and low Gini, and for cells with high exp and low exp
which(GiniUMI==max(GiniUMI)) # SER8_TGAACTCATATA
which(GiniUMI==min(GiniUMI)) # SER8_GCGCGTTCCAAG
test=as.matrix(dge@raw.data)[,which(GiniUMI==max(GiniUMI) | GiniUMI==min(GiniUMI) | dge@data.info$nUMI==max(dge@data.info$nUMI) | dge@data.info$nUMI==min(dge@data.info$nUMI) )]
write.table(test,"test4cells.txt",quote=F,sep="\t",row.names=T,col.names=T)

pdf(paste0(dgefile,"LorenzCurve.pdf"))
plot.Lc(dge@raw.data[,which(GiniUMI==max(GiniUMI))],main=paste("Cell with Highest Gini",max(GiniUMI)))
plot.Lc(dge@raw.data[,which(GiniUMI==min(GiniUMI))],main=paste("Cell with Lowest Gini",min(GiniUMI)))
plot.Lc(dge@raw.data[,which(dge@data.info$nUMI==max(dge@data.info$nUMI))],main=paste("Cell with Highest Exp",max(dge@data.info$nUMI)))
plot.Lc(dge@raw.data[,which(dge@data.info$nUMI==min(dge@data.info$nUMI))],main=paste("Cell with Lowest Exp",min(dge@data.info$nUMI)))
dev.off()

test=read.table("test4cells.txt",header=T,row.names=1)

low=paste0("LowestExpCell TotalUMI=566 Gini=",round(ineq(test[,1]),2)," %Non-0Genes=",sprintf("%.2g%%",sum(test[,1]!=0)/n*100))
high=paste0("HighestExpCell TotalUMI=183145 Gini=",round(ineq(test[,2]),2)," %Non-0Genes=",sprintf("%.2g%%",sum(test[,2]!=0)/n*100))
Lc.p <- Lc(test[,1])
Lc.u <- Lc(test[,2])
plot(Lc.p,col=4)
lines(Lc.u, col=2)
legend("topleft",legend=c("LowestExpCell","HighestExpCell"),col=c(4,2),text.col=c(4,2),lty=1,bty="n")

plot(Lc.p, general=TRUE,col=4)
lines(Lc.u, general=TRUE, col=2)

low=paste0("LowestGiniCell Gini=",round(ineq(test[,3]),3)," TotalUMI=25070"," %Non-0Genes=",sprintf("%.2g%%",sum(test[,3]!=0)/n*100))
high=paste0("HighestGiniCell Gini=",round(ineq(test[,4]),3)," TotalUMI=2628"," %Non-0Genes=",sprintf("%.2g%%",sum(test[,4]!=0)/n*100))
Lc.p <- Lc(test[,3])
Lc.u <- Lc(test[,4])
plot(Lc.p,col=3)
lines(Lc.u, col=6)
legend("topleft",legend=c("LowestGiniCell","HighestGiniCell"),col=c(3,6),text.col=c(3,6),lty=1,bty="n")

### Try to Compute standard Lorenz curve for Gini=0, 0.2, 0.5, 0.8, 1
dim(test)
n=dim(test)[1]
a=c(rep(0,n*0.9),(n*0.9):n)         # Gini=0.9
b=1:log(n)    # Gini=0.3
c=(1:n)^2     # Gini=0.5
ineq(a)
plot(Lc(a))
lines(Lc(b),col=2)
lines(Lc(c),col=3)
legend("topleft",legend=paste0("Gini=",round(c(ineq(a),ineq(b),ineq(c)),2)),col=1:3,lty=1)

a=sqrt(1:n)   # Gini=0.2
b=1:log10(n)  # Gini=0.25
c=(1:n)^3     # Gini=0.6
ineq(a)
plot(Lc(a))
lines(Lc(b),col=2)
lines(Lc(c),col=3)
legend("topleft",legend=paste0("Gini=",round(c(ineq(a),ineq(b),ineq(c)),2)),col=1:3,lty=1)

a=rep(1,n)   # Gini=0
b=(1:n)^6    # Gini=0.75
c=(1:n)^8    # Gini=0.8
ineq(a)
plot(Lc(a))
lines(Lc(b),col=2)
lines(Lc(c),col=3)
legend("topleft",legend=paste0("Gini=",round(c(ineq(a),ineq(b),ineq(c)),2)),col=1:3,lty=1)

a=c(rep(0,n-1),n) # Gini=1
b=(1:n)^1.35    # Gini=0.4
c=(1:n)^8    # Gini=0.8
ineq(a)
plot(Lc(a))
lines(Lc(b),col=2)
lines(Lc(c),col=3)
legend("topleft",legend=paste0("Gini=",round(c(ineq(a),ineq(b),ineq(c)),2)),col=1:3,lty=1)

### Plot Standard Lorenz Curve for Gini=0,0.2,0.4,0.5,0.6,0.8,1
a=rep(1,n)    # Gini=0
b=sqrt(1:n)   # Gini=0.2
c=(1:n)^1.35  # Gini=0.4
d=(1:n)^2     # Gini=0.5
e=(1:n)^3     # Gini=0.6
f=(1:n)^8     # Gini=0.8
g=c(rep(0,n*0.9),(n*0.9):n)         # Gini=0.9
h=c(rep(0,n-1),n) # Gini=1
plot(Lc(a))
lines(Lc(b),col=2)
lines(Lc(c),col=3)
lines(Lc(d),col=4)
lines(Lc(e),col=5)
lines(Lc(f),col=6)
lines(Lc(g),col="gold")
lines(Lc(h),col="gray60")
GiniCoef=paste0("Gini=",round(c(ineq(a),ineq(b),ineq(c),ineq(d),ineq(e),ineq(f),ineq(g),ineq(h)),2))
formula=c("rep(1,n)","sqrt(1:n)","(1:n)^1.35","(1:n)^2","(1:n)^3","(1:n)^8","c(rep(0,n*0.9),(n*0.9):n)","c(rep(0,n-1),n)")
legend("topleft",legend=paste(GiniCoef," ",formula),col=c(1:6,"gold","gray60"),text.col=c(1:6,"gold","gray60"),lty=1)

GiniGeneVsUMI=GiniUMIVsReads=rep(0,10)
for(j in 1:10){
  GiniGeneVsUMI[j]=ineq((datainfo$nGene/datainfo$nUMI)[which(datainfo$celltype==j)],type=c("Gini")) 
  GiniUMIVsReads[j]=ineq((datainfo$nUMI/datainfo$nReads)[which(datainfo$celltype==j)],type=c("Gini")) 
}
barplot(GiniGeneVsUMI,col=myBrewerPalette)
barplot(GiniUMIVsReads,col=myBrewerPalette)


###### Plot distribution of Gini index for each cell using Violin Plot - Figure 1D right-most panel
datainfo$nUMIVsnGene=datainfo$nUMI/datainfo$nGene
datainfo$nReadsVsnUMI=datainfo$nReads/datainfo$nUMI
pdf(file=paste(dgefile,"Gini_11clusters.pdf",sep=""),height=4,width=8)
VlnPlot(dge, "nUMIVsnGene", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "nReadsVsnUMI", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "GiniAll", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "GiniNon0", nCol = 1,cols.use=myBrewerPalette)
VlnPlot(dge, "GiniHVG", nCol = 1,cols.use=myBrewerPalette)
dev.off()
### save as Figure 1D right-most panel

###### Correlation between Gini index and cell size factor - Figure S1F
datainfo=dge@data.info
datainfo$Gini=GiniUMI
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
datainfo$celltype=factor(datainfo$celltype,levels=levels)

### Gini index Vs nUMI/nGene
## scatter plot of nUMI/nGene (x-axis) and Gini index (y-axis) for each cell type 
# plotted each cell type separately, adding all others as background
datainfo_bg=datainfo[,! names(datainfo) %in% "celltype"]
ggplot(datainfo,aes(x=nUMIVsnGene,y=Gini,colour=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2)+
geom_point() +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')

## correlation between Gini index and nUMI/nGene
cor.test(datainfo$Gini,datainfo$nUMIVsnGene)
# within-cell type correlation
for(i in 1:11){
print(cor(datainfo$Gini[which(datainfo$celltype==levels[i])],datainfo$nUMIVsnGene[which(datainfo$celltype==levels[i])]))
}
for(i in 1:11){
print(cor.test(datainfo$Gini[which(datainfo$celltype==levels[i])],datainfo$nUMIVsnGene[which(datainfo$celltype==levels[i])])$p.val)
}
# between cell type correlation
tmpGini=tmp1=tmp2=rep(0,11)
for(i in 1:11){
  tmpGini[i]=mean(datainfo[which(datainfo$celltype==levels[i]),]$Gini)
  tmp1[i]=mean(datainfo[which(datainfo$celltype==levels[i]),]$nUMIVsnGene)
  tmp2[i]=mean(datainfo[which(datainfo$celltype==levels[i]),]$nUMI)
}
cor.test(tmpGini,tmp1)
cor.test(tmpGini,tmp2)

### Gini index Vs nUMI - Figure S1F
## scatter plot of nUMI (x-axis) and Gini index (y-axis) for each cell type - Figure S1F
# plotted each cell type separately, adding all others as background
pdf(file=paste(dgefile,"GiniVsnUMI_11clusters.pdf",sep=""),height=4,width=6.5)
# Gini for all genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniAll,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~celltype,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for non-0 genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniNon0,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~celltype,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
# Gini for highly-variable genes of each cell
ggplot(datainfo,aes(x=nUMI,y=GiniHVG,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.4)+
geom_point(alpha=0.8,size=0.4) +
facet_wrap(~celltype,ncol=4) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
dev.off()
## Save as Figure S1F

## Correlation between Gini index and nUMI for each cell
cor.test(datainfo$Gini,datainfo$nUMI)
for(i in 1:11){
print(cor(datainfo$Gini[which(datainfo$celltype==levels[i])],datainfo$nUMI[which(datainfo$celltype==levels[i])]))
}
for(i in 1:11){
print(cor.test(datainfo$Gini[which(datainfo$celltype==levels[i])],datainfo$nUMI[which(datainfo$celltype==levels[i])])$p.val)
}

### scatter plot of Gini index Vs nGene
ggplot(datainfo,aes(x=nGene,y=Gini,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2)+
geom_point(size=1) +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()

### Correlation bewteen nUMI Vs nGene
## scatter plot of nUMI (x-axis) and nGene (y-axis)
ggplot(datainfo,aes(x=nUMI,y=nGene,color=celltype))+
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
scale_x_continuous(trans = 'log10')
pdf(file=paste(dgefile,"NGeneUMImt_11clusters2.pdf",sep=""),height=4,width=6)
# nGene ~ nUMI
ggplot(dge@data.info, aes(x=nUMI,y=nGene,color=celltype)) +
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(breaks=c(566,50000,100000,150000)) 
# nGene ~ log10 nUMI
ggplot(dge@data.info, aes(x=nUMI,y=nGene,color=celltype)) +
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(trans = 'log10',
                   breaks = trans_breaks('log10', function(x) 10^x),
                   labels = trans_format('log10', math_format(10^.x)))
dev.off()
## scatter plot of nGene Vs nUMI for each cell type
pdf(file=paste(dgefile,"NGeneUMImt_11clusters_facet.pdf",sep=""),height=6,width=15)
ggplot(dge@data.info, aes(x=nUMI,y=nGene,color=celltype)) +
geom_point(size=1) +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(breaks=c(566,50000,100000,150000)) 
dev.off()
## correlation between nUMI and nGene
cor.test(datainfo$nUMI,datainfo$nGene)

### scatter plot of nReads and nGene
pdf(file=paste(dgefile,"NGeneReads_11clusters2.pdf",sep=""),height=4,width=6)
# nGene ~ nReads
ggplot(dge@data.info, aes(x=nReads,y=nGene,color=celltype)) +
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(breaks=c(821,seq(250000,750000,250000))) 
# nGene ~ log10 nReads
ggplot(dge@data.info, aes(x=nReads,y=nGene,color=celltype)) +
geom_point(size=1) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(trans = 'log10',
                   breaks = trans_breaks('log10', function(x) 10^x),
                   labels = trans_format('log10', math_format(10^.x)))
dev.off()
# nGene Vs nReads for each cell type
pdf(file=paste(dgefile,"NGeneReads_11clusters_facet.pdf",sep=""),height=6,width=15)
ggplot(dge@data.info, aes(x=nReads,y=nGene,color=celltype)) +
geom_point(size=1) +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette) +
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
scale_y_continuous(breaks=c(501,1000,seq(2500,10000,2500))) +
scale_x_continuous(trans = 'log10',
                   breaks = trans_breaks('log10', function(x) 10^x),
                   labels = trans_format('log10', math_format(10^.x)))
dev.off()


### Per-cell type summary
## nUMI per Cell
for(j in 1:11){
	print(summary(dge@data.info$nUMI[which(dge@ident==j)]))
}
## %MT per cell
for(j in 1:11){
	print(summary(dge@data.info$percent.mito[which(dge@ident==j)]))
}
## nGene per cell
for(j in 1:11){
  print(summary(dge@data.info$nGene[which(dge@ident==j)]))
}
## %X per cell
for(j in 1:11){
  print(summary(dge@data.info$percent.x[which(dge@ident==j)]))
}



######### 5.8.2018 %ChrX and %ChrY distribution in round spermatids and other cell types - Figure S1G
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
dgefile="/home/qzm/figDec2017_MouseAdultST25/All_"
datainfo=read.table(paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),header=T,row.names=1,stringsAsFactors=F)
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
datainfo$celltype=factor(datainfo$celltype,levels=levels)

library(RColorBrewer)
myBrewerPalette <- c(brewer.pal(12,"Paired"))[c(12,11,10,6,9,5,7,1:4)] # 1
library(ggplot2)
library(ggExtra)

###### plot distribution of %ChrX and %ChrY across 11 major cell types
### plot each cell type separately, adding all others as background
datainfo_bg=datainfo[,! names(datainfo) %in% "celltype"]
### plot in linear scale
jpeg(file=paste0(dgefile,"ChrXY.jpeg"),res=300,height=1400,width=2600)
ggplot(datainfo,aes(x=percent.x,y=percent.y,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(size=0.8) +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
xlab("%ChrX")+
ylab("%ChrY")
dev.off()
# linear scale cannot show the details
### plot in log10 scale
jpeg(file=paste0(dgefile,"ChrXY_log10.jpeg"),res=300,height=1400,width=2600)
ggplot(datainfo,aes(x=percent.x,y=percent.y,color=celltype))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(size=0.8) +
facet_wrap(~celltype,ncol=5) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="Cell Type",override.aes = list(size=2))) +
theme_bw()+
xlab("%ChrX")+
ylab("%ChrY")+
scale_x_continuous(trans = 'log10')+
scale_y_continuous(trans = 'log10')
dev.off()
# Warning: Transformation introduced infinite values in continuous axis 

### replace 0 by 5E-6 to avoid infinite values after log-transformation
min(datainfo$percent.x[which(datainfo$percent.x!=0)]) # [1] 2.97274e-05
min(datainfo$percent.y[which(datainfo$percent.y!=0)]) # [1] 1.083459e-05
datainfo$percent.x[which(datainfo$percent.x==0)]<-5E-6
datainfo$percent.y[which(datainfo$percent.y==0)]<-5E-6
datainfo_bg=datainfo[,! names(datainfo) %in% "celltype"]

### plot in log10 scale, semi-transparent dots, marginal density plot - Figure S1G
# in order to visualize the density of cell from specific cell type, need to have two groups: this cell type vs all other cells
for(i in 1:length(levels)){

	level=levels[i]
	col=myBrewerPalette[i]
	tmp1=datainfo[which(datainfo$celltype==level),]
	tmp2=datainfo
	tmp1$group=level
	tmp2$group="All"
	tmp=rbind(tmp1,tmp2)
	tmp$group=factor(tmp$group,levels=c(level,"All"))

	p1 <- ggplot(tmp,aes(x=percent.x,y=percent.y,colour=group,fill=group))+
	geom_point(alpha=0.1,size=0.8)+
	scale_color_manual(values = c(col,"grey"))+
	scale_fill_manual(values = c(col,"grey"))+
	geom_point(data=tmp1,colour=col,fill=col,alpha=0.2,size=0.8) +
	theme_bw()+
	ggtitle(level)+
	xlab("%ChrX")+
	ylab("%ChrY")+
	scale_x_continuous(trans = 'log10')+
	scale_y_continuous(trans = 'log10')+
	theme(legend.position="none")

	jpeg(file=paste0(dgefile,"ChrXY_log10_",level,".jpeg"),res=300,height=1000,width=950)
	ggMarginal(p=p1,groupColour=TRUE,groupFill=TRUE,alpha=0.5)
	dev.off()

}
## save as Figure S1G left-most panel

### plot %ChrX and %ChrY for each germ cell type with 4 stratified UMI groups: <1k, 1k-3k, 3k-5k, >5k - Figure S1G
for(i in 8:11){
for(j in c(500,1000,3000,5000)){
	if(j==500){
		j2=1000
	} else if(j==1000){
		j2=3000;label="1k-3k"
	} else if(j==3000){
		j2=5000;label="3k-5k"
	} else {
		j2=max(datainfo$nUMI);label=">5k"
	}

	datainfotmp=datainfo[which(datainfo$nUMI >j & datainfo$nUMI <= j2),]
	level=levels[i]
	col=myBrewerPalette[i]
	tmp1=datainfotmp[which(datainfotmp$celltype==level),]
	if(j==500){
		label=paste0(min(tmp1$nUMI),"-1k")
	}
	tmp2=datainfo
	tmp1$group=level
	tmp2$group="All"
	tmp=rbind(tmp1,tmp2)
	tmp$group=factor(tmp$group,levels=c(level,"All"))

	p1 <- ggplot(tmp,aes(x=percent.x,y=percent.y,colour=group,fill=group))+
	geom_point(alpha=0.1,size=0.8)+
	scale_color_manual(values = c(col,"grey"))+
	scale_fill_manual(values = c(col,"grey"))+
	geom_point(data=tmp1,colour=col,fill=col,alpha=0.5,size=0.8) +
	theme_bw()+
	ggtitle(paste(level,label))+
	xlab("%ChrX")+
	ylab("%ChrY")+
	scale_x_continuous(trans = 'log10')+
	scale_y_continuous(trans = 'log10')+
	theme(legend.position="none")

	jpeg(file=paste0(dgefile,"ChrXY_log10_",level,"_nUMI",label,".jpeg"),res=300,height=1000,width=950)
	# ggMarginal(p=p1,colour=col,fill=col,alpha=0.5)
	ggMarginal(p=p1,groupColour=TRUE,groupFill=TRUE,alpha=0.5)
	dev.off()
}
}
## save as Figure S1G right panels

### calculate number and %cells with no detectable ChrX or ChrY genes for each germ cell type
datainfo=read.table(paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),header=T,row.names=1,stringsAsFactors=F)
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
datainfo$celltype=factor(datainfo$celltype,levels=levels)
for(i in 8:11){
for(j in c(500,1000,3000,5000)){
	if(j==500){
		j2=1000
	} else if(j==1000){
		j2=3000
	} else if(j==3000){
		j2=5000
	} else {j2=max(datainfo$nUMI)}
	datainfotmp=datainfo[which(datainfo$nUMI >j & datainfo$nUMI <= j2),]
	# print(c(j,j2,dim(datainfotmp)))
	level=levels[i]
	tmp1=datainfotmp[which(datainfotmp$celltype==level),]
	x0=length(which(tmp1$percent.x==0))
	x0f=length(which(tmp1$percent.x==0))/nrow(tmp1)
	y0=length(which(tmp1$percent.y==0))
	y0f=length(which(tmp1$percent.y==0))/nrow(tmp1)
	print(c(x0,x0f,y0,y0f))
}
}
## save the number and % in Figure S1G

### check odds ratio for round spermatids cells with no detectable ChrX or ChrY genes
i=10
	datainfotmp=datainfo
	level=levels[i]
	col=myBrewerPalette[i]
	tmp1=datainfotmp[which(datainfotmp$celltype==level),]
	tmp2=datainfotmp[which(datainfotmp$celltype!=level),]

	m0=length(which(tmp1$percent.x==0))
	m1=length(which(tmp1$percent.x!=0))
	a0=length(which(tmp2$percent.x==0))
	a1=length(which(tmp2$percent.x!=0))
	print(c(m0,m1,a0,a1))
#[1]     4  9919   788 23922

	m0=length(which(tmp1$percent.y==0))
	m1=length(which(tmp1$percent.y!=0))
	a0=length(which(tmp2$percent.y==0))
	a1=length(which(tmp2$percent.y!=0))
	print(c(m0,m1,a0,a1))
#[1]  1390  8533 14468 10242


##### for Scyte with >5k UMIs, divide cells into 3 groups based on %ChrX 
i=9
	j=5000;j2=max(datainfo$nUMI);label=">5k"

datainfo=read.table(paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),header=T,row.names=1,stringsAsFactors=F)
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
datainfo$celltype=factor(datainfo$celltype,levels=levels)

	datainfotmp=datainfo[which(datainfo$nUMI >j & datainfo$nUMI <= j2),]
	level=levels[i]
	col=myBrewerPalette[i]
	tmp1=datainfotmp[which(datainfotmp$celltype==level),]

## divide Scyte with >5k UMIs to 4 groups based on %ChrX (3 groups with detectable ChrX genes, 1 group with no detectable ChrX genes)	
	summary(tmp1$percent.x[which(tmp1$percent.x!=0)])
	#2.973e-05 2.550e-04 4.623e-04 1.536e-03 1.493e-03 1.510e-02
	q1=quantile(tmp1$percent.x[which(tmp1$percent.x!=0)],.35)
	q2=quantile(tmp1$percent.x[which(tmp1$percent.x!=0)],.70)
	q1=quantile(tmp1$percent.x[which(tmp1$percent.x!=0)],.7)
	q2=quantile(tmp1$percent.x[which(tmp1$percent.x!=0)],.88)
	c1=rownames(tmp1)[which(tmp1$percent.x==0)]
	c2=rownames(tmp1)[which(tmp1$percent.x>0 & tmp1$percent.x<=q1)]
	c3=rownames(tmp1)[which(tmp1$percent.x>q1 & tmp1$percent.x<=q2)]
	c4=rownames(tmp1)[which(tmp1$percent.x>q2)]
	scyte5kUMIs=c(c1,c2,c3,c4)
	xgroup=c(rep(1,length(c0)),rep(2,length(c2)),rep(3,length(c3)),rep(4,length(c4)))
	names(xgroup)=scyte5kUMIs
	xgroup=as.factor(xgroup)
	table(xgroup)
 #  1    2    3    4
 #159 1970 1969 1688
 #  1    2    3    4
 #159 3939 1012  676

## scatterplot and density plot for stratfied %ChrX for Scyte with >5k UMIs
	tmp1$xgroup=xgroup[rownames(tmp1)]
 	summary(tmp1$percent.x[which(tmp1$xgroup==1)])
  	summary(tmp1$percent.x[which(tmp1$xgroup==2)])
 	summary(tmp1$percent.x[which(tmp1$xgroup==3)])
 	summary(tmp1$percent.x[which(tmp1$xgroup==4)])

	tmp2=datainfo
	tmp2$xgroup="All"
	tmp=rbind(tmp1,tmp2)
	tmp$xgroup=factor(tmp$xgroup,levels=c(1:4,"All"))
	tmp$percent.x[which(tmp$percent.x==0)]<-5E-6
	tmp$percent.y[which(tmp$percent.y==0)]<-5E-6
	tmp1$percent.x[which(tmp1$percent.x==0)]<-5E-6
	tmp1$percent.y[which(tmp1$percent.y==0)]<-5E-6


	p1 <- ggplot(tmp,aes(x=percent.x,y=percent.y,colour=xgroup,fill=xgroup))+
	geom_point(alpha=0.1,size=0.8)+
	scale_color_manual(values = c("yellow","green","blue","red","grey"))+
	scale_fill_manual(values = c("yellow","green","blue","red","grey"))+
	geom_point(data=tmp1,alpha=0.5,size=0.8) +
	theme_bw()+
	ggtitle("Scyte >5kUMIs %ChrX")+
	xlab("%ChrX")+
	ylab("%ChrY")+
	scale_x_continuous(trans = 'log10')+
	scale_y_continuous(trans = 'log10')+
	theme(legend.position="none")

	jpeg(file=paste0(dgefile,"ChrXY_log10_",level,"_nUMI",label,"_ChrXgroups.jpeg"),res=300,height=1000,width=950)
	ggMarginal(p=p1,groupColour=TRUE,groupFill=TRUE,alpha=0.5)
	dev.off()


###### PCA plot for stratfied %ChrX for Scyte with >5k UMIs
PCs=read.table(paste0(home,"MouseAdultST25/Germ1kGenes/Germ_PC1-3.txt"),row.names=1)
others=rownames(PCs)[which(!rownames(PCs) %in% scyte5kUMIs)]

pdf(paste(dgefile,"ScyteOver5kUMIs_4ChrXgroups.pdf",sep=""),height=8,width=8)
par(mfrow=c(2,2),mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(PCs[others,1],PCs[others,3],col=rgb(0,0,0,0.1),pch=16,cex=0.5,xlab="PC1",ylab="PC2")
points(PCs[c1,1],PCs[c1,3],col=rgb(1,1,0,0.6),pch=16,cex=0.5)
legend("bottomleft",pch=16,col=rgb(1,1,0,0.6),legend="%ChrX == 0")
plot(PCs[others,1],PCs[others,3],col=rgb(0,0,0,0.1),pch=16,cex=0.5,xlab="PC1",ylab="PC2")
points(PCs[c2,1],PCs[c2,3],col=rgb(0,1,0,0.3),pch=16,cex=0.5)
legend("bottomleft",pch=16,col=rgb(0,1,0,0.6),legend=paste("%ChrX: 0 -",sprintf("%.2f%%",q1*100)))
plot(PCs[others,1],PCs[others,3],col=rgb(0,0,0,0.1),pch=16,cex=0.5,xlab="PC1",ylab="PC2")
points(PCs[c3,1],PCs[c3,3],col=rgb(0,0,1,0.3),pch=16,cex=0.5)
legend("bottomleft",pch=16,col=rgb(0,0,1,0.6),legend=paste("%ChrX:",sprintf("%.2f%%",q1*100),"-",sprintf("%.2f%%",q2*100)))
plot(PCs[others,1],PCs[others,3],col=rgb(0,0,0,0.1),pch=16,cex=0.5,xlab="PC1",ylab="PC2")
points(PCs[c4,1],PCs[c4,3],col=rgb(1,0,0,0.3),pch=16,cex=0.5)
legend("bottomleft",pch=16,col=rgb(1,0,0,0.6),legend=paste("%ChrX:",sprintf("%.2f%%",q2*100),"- 1.5%"))
dev.off()

pdf(paste(dgefile,"ScyteOver5kUMIs_4ChrXgroups1.pdf",sep=""),height=5,width=5)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(PCs[others,1],PCs[others,3],col=rgb(0,0,0,0.1),pch=16,cex=0.5,xlab="PC1",ylab="PC2")
points(PCs[scyte5kUMIs,1],PCs[scyte5kUMIs,3],col=c(rgb(1,1,0,0.6),rgb(0,1,0,0.3),rgb(0,0,1,0.3),rgb(1,0,0,0.3))[xgroup],pch=16,cex=0.5)
legend("topright",pch=16,col=c("yellow","green","blue","red"),legend=c("%ChrX == 0",paste("%ChrX: 0 -",sprintf("%.2f%%",q1*100)),paste("%ChrX:",sprintf("%.2f%%",q1*100),"-",sprintf("%.2f%%",q2*100)),paste("%ChrX:",sprintf("%.2f%%",q2*100),"- 1.5%")))
dev.off()


###### heatmap for %ChrX in PCA
datainfo=read.table(paste0(home,"data_DGE/MouseAdultST25_datainfo.txt"),header=T,row.names=1,stringsAsFactors=F)
levels=c("InnateLymphoid","Macrophage","Endothelial","Myoid","Leydig","Sertoli","Unknown","Spermatogonia","Spermatocyte","RoundSpermatid","Elongating")
datainfo$celltype=factor(datainfo$celltype,levels=levels)

PCs=read.table(paste0(home,"MouseAdultST25/Germ1kGenes/Germ_PC1-3.txt"),row.names=1)

dgefile=dgename="/home/qzm/figDec2017_MouseAdultST25/"
redblue100.alpha<-rgb(read.table("/scratch/junzli_flux/qzm/Dropseq_analysis/data_DGE/redblue100.txt",sep='\t',row.names=1,header=T),alpha=0.8)

setname="Germ1kGenes"
dim.1 = 1; dim.2 = 3; cells.use = NULL; pt.size = 1;
pch.use = 16; reduction.use = "PCA";
use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE

data.plot=PCs[,c(1,3)]
data.plot$x=data.plot[,1]
data.plot$y=data.plot[,2]
data.plot$pt.size <- pt.size

data.use <- data.frame(t(datainfo[rownames(data.plot),]$percent.x))
rownames(data.use)=feature=features.plot="percent.x"

data.gene0 <- na.omit(data.frame(data.use[features.plot, ]))
data.plot$gene <- t(data.gene0)

st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

color=redblue100.alpha[c(1:50,rep(51:80,each=3),rep(81:101,each=5))]
ll=length(color)-1
jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8_red10.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- color[(z - min(z))/diff(range(z))*ll + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()
image(cbind(seq(0,1,by=0.01),seq(0,1,by=0.01)),col=zcolor)


##### use red-to-blue for spermatocytes GC4-8 for %ChrX heatmap
ident=read.table(paste0(home,"MouseAdultST25/Germ1kGenes/Germ_ident.txt"),row.names=1)

st0<- data.plot

scyte=rownames(ident)[which(ident[,1]>=4 & ident[,1]<=8)]
length(scyte) # 7729
st6=st0[scyte,]
z=st6[,6]

zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8_scyteBlueToRed.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
plot(st0[,1],st0[,2], col=rgb(0,0,0,0.1),pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
points(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3)
dev.off()

# convert to log2 scale
z=st6[,6]
min(z[which(z!=0)]) # [1] 2.97274e-05
z[which(z==0)] <- 2.5E-5
z=log2(z)
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
jpeg(paste0(dgename,setname,"_",feature,"log2_redblue0.8_scyteBlueToRed.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
plot(st0[,1],st0[,2], col=rgb(0,0,0,0.1),pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
points(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3)
dev.off()

# convert to log10 scale
z=st6[,6]
min(z[which(z!=0)]) # [1] 2.97274e-05
z[which(z==0)] <- 2.5E-5
z=log10(z)
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
jpeg(paste0(dgename,setname,"_",feature,"log10_redblue0.8_scyteBlueToRed.jpeg"),height=1700,width=1600,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
plot(st0[,1],st0[,2], col=rgb(0,0,0,0.1),pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
points(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3)
dev.off()
image(cbind(seq(min(z),max(z),by=0.01),seq(min(z),max(z),by=0.01)),col=redblue100.alpha)
summary(z)
10^quantile(z,0) # note: I converted 0 to 2.5E-5; so this should be 0
10^quantile(z,0.2)*100
10^quantile(z,0.4)*100
10^quantile(z,0.6)*100
10^quantile(z,0.8)*100
10^quantile(z,1)*100


###### plot distribution of %ChrX and %ChrY across 12 Germ Cell Clusters (GC1-12)
home="/scratch/junzli_flux/qzm/Dropseq_analysis/"
ident=read.table(paste0(home,"MouseAdultST25/Germ1kGenes/Germ_ident.txt"),row.names=1)
cluster=factor(paste0("GC",ident[,1]),levels=paste0("GC",1:12),order=T)
names(cluster)=rownames(ident)
datainfo=datainfo[names(cluster),]
datainfo$cluster=cluster
library(RColorBrewer)
myBrewerPalette=brewer.pal(12,"Paired")

# plot in linear scale
jpeg(file=paste0(dgefile,"ChrXY_GC1-12.jpeg"),res=300,height=900,width=2600)
ggplot(datainfo,aes(x=percent.x,y=percent.y,color=cluster))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(size=0.8) +
facet_wrap(~cluster,ncol=6) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="",override.aes = list(size=2))) +
theme_bw()+
xlab("%ChrX")+
ylab("%ChrY")
dev.off()

# plot in log10 scale
jpeg(file=paste0(dgefile,"ChrXY_log10_GC1-12.jpeg"),res=300,height=900,width=2600)
ggplot(datainfo,aes(x=percent.x,y=percent.y,color=cluster))+
geom_point(data=datainfo_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(size=0.8) +
facet_wrap(~cluster,ncol=6) +
scale_color_manual(values=myBrewerPalette)+
guides(color=guide_legend(title="",override.aes = list(size=2))) +
theme_bw()+
xlab("%ChrX")+
ylab("%ChrY")+
scale_x_continuous(trans = 'log10')+
scale_y_continuous(trans = 'log10')
dev.off()

