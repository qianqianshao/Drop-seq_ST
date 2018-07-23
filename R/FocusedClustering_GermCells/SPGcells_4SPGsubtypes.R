### R script to generate 4 SPG subtypes by Qianyi on 6/28/2017
### Related to Figure 4 and Table S5: SPG cells (N=2,484) with >1k UMIs

### load libraries
library(Seurat)
library(dplyr)
library(Matrix)

### analysis for SPG cells with >1k UMIs
sets=list(c(1:3));j=1
setsname="123cellsUMI1k"
setslabel="SpermatogoniaUMI1k"
dgefile=dgename="figJul2017_MouseAdultST24_ReCluster123filter_6_9/"

### load object of 24 ST datasets (INT6 dataset did not contribute to SPG cells)
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15.Robj"))

##### zoomed-in view of PCA for SPG subset with >1k UMIs using PCA of 24 ST
cells.use=cells.use[which(dge@data.info$nUMIperCell2[which(dge@ident %in% sets[[j]])]>1000)]
pdf(paste(dgefile,"ZoomedInView_PCA_SPGsubset.pdf",sep=""),height=10,width=18)
PCAPlot(dge,1,2,cells.use=cells.use,do.return = TRUE,pt.size = 1,do.label=T)
PCAPlot(dge,1,3,cells.use=cells.use,do.return = TRUE,pt.size = 1,do.label=T)
dev.off()

######### Focused clustering for SPG cells subset (N=1,067) with >1k UMIs on 6/28/2017
###### Extract cells in SPG subset
dgetmp=SubsetData(dge,ident.use=sets[[j]])
### Remove SPG cells with <=1k UMIs in SPG subset
  cells.use=rownames(dge@data.info)[which(dge@data.info$nUMIperCell2>1000)]
  dgetmp=SubsetData(dge,cells.use=cells.use)
  dgetmp
  dgetmp=SetAllIdent(dgetmp,id="clusters31_GeneFilter1")
  print(table(dgetmp@ident))
  ncellsbatchcluster = table(dgetmp@data.info$orig.ident,dgetmp@ident)
  print(ncellsbatchcluster)
  print(prop.table(ncellsbatchcluster,2))
  print(prop.table(ncellsbatchcluster,1))
### Extract genes detected in SPG cells with >1k UMIs
nCellperGene <- rowSums(as.matrix(dgetmp@data)>0)
genes.use=names(nCellperGene[which(nCellperGene!=0)])
dgetmp@data=dgetmp@data[genes.use,]
print(table(dgetmp@ident))
print(dim(dgetmp@data))
### PCA for SPG subset
print(Sys.time())
dgetmp <- PCA(dgetmp, pc.genes = rownames(dgetmp@data), do.print = TRUE, pcs.print = 5, genes.print = 5)
print(Sys.time())
save(dgetmp,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster",setsname[j],".Robj"))

### plot PCA
pdf(paste(dgefile,"ZoomedInReDo_PCA_ReCluster_",setsname[j],".pdf",sep=""),height=5,width=6)
PCAPlot(dgetmp,1,2,do.return = TRUE,pt.size = 1,do.label=T)
PCAPlot(dgetmp,1,3,do.return = TRUE,pt.size = 1,do.label=T)
dev.off()

### Scree Plot for single PCA
numPCs=10
pdf(paste(dgefile,"dge_PCA_Variablel_variation.pdf",sep=""))
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
plot(dge@pca.obj[[1]]$sdev[1:150],type="b",ylab="Eigenvalue",xlab="PC",cex.lab=1.5) #check components representing greatest variation
abline(v=numPCs[i]+0.5,col="red")
text(numPCs[i]+0.5,0.8,col="red",paste(numPCs[i],"PCs"))
dev.off()
### density plot of Eigenvalue
eigenvalue=dge@pca.obj[[1]]$sdev[numPCs[i]]
eigenvalue # [1] 5.456186
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
# used res.0.2 for 4 clusters

### Build Cluster tree 
res="res.0.2";resi=i=1
dge <- SetAllIdent(dge, id = res[i])
dge <- BuildClusterTree(dge, do.reorder = T, reorder.numeric = T,pcs.use=1:numPCs[i])
dge@data.info[,res[i]]=dge@ident

###### Order the SPG clusters 
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

### Reordering cluster centroids by seriation
library(seriation)
n=ncluster=length(levels)
tmp=genecountsall[,levels]
da <- dist(t(as.matrix(tmp)), method = "euclidean")
pdf(file=paste0(dgename,"Centroid_norm_Seriation_Dissimilarity_",res[resi],".pdf"))
 dissplot(da, method="OLO")
 hmap(da) # # default method="OLO"
dev.off()
levelss=get_order(seriate(da,method="OLO"))
levels=levelss

### Reordered clusters for all cells
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

### save the re-ordered cluster IDs to dge file
ordered="res.0.2order" # this is the 4 SPG subtypes
dge=AddMetaData(dge,cells.ident.ordered,ordered)
dge=SetAllIdent(dge,ordered)
save(dge,file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster123cellsUMI1k.Robj"))

###### visualize 4 SPG subtypes in PCA & tSNE - Figure 4A
pdf(paste0(dgefile,"cluster_4SPGsubtypes.pdf"),width=5.5,height=5)
print(table(dge@data.info$clusters31_GeneFilter1,dge@ident))
PCAPlot(dge,pt.size = 1,do.label=TRUE)
TSNEPlot(dge,pt.size = 1,do.label=TRUE,label.size=8)  # Figure 4A
dev.off()
# save as Figure 4A

### Contribution of each batch to each cluster - Table S5A
  ## Absolute Number of Cells in Each Batch and Each Cluster
  ncellsbatchcluster = table(dge@data.info$orig.ident,dge@ident)
  ## % Cells Contributed to a Single Cluster from Each Batch
  percentcellsbatch = prop.table(ncellsbatchcluster,2)
  ## % Cells in Each Cluster from a Single Batch
  percentcellscluster = prop.table(ncellsbatchcluster,1)
  ## save as tables
  ncellscluster=rbind(ncellsbatchcluster,percentcellsbatch,percentcellscluster)
  write.table(ncellscluster,paste0(dgefile,"Cluster",setsname[j],res[i],"_ncellspercluster_batch.txt"),quote=F,row.names=T,col.names=T,sep="\t")
# save as Table S5A

### visualize 4 SPG subtypes for each individual batch in tSNE view
load(file =paste0(home,"data_DGE/MouseAdultST24genesUMI20cell15_31clusters_cluster123cellsUMI1k.Robj"))
dgefile=dgename="figJul2017_MouseAdultST24_ReCluster123filter_6_9/SPG1kUMIs_"
#dge <- SetAllIdent(dge, id = "res.0.2order")     # used this
table(dge@ident)
#  1   2   3   4 
#213 730 950 591 
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
myBrewerPalette=gg_color_hue(4)
# go to DimPlot to set xlim ylim
range(dge@tsne.rot[,1])
range(dge@tsne.rot[,2])
### tSNE: +coord_cartesian(ylim=c(-29,30),xlim=c(-46,44))
### remove legends
### plot all 2,484 SPG cells with >1k UMIs from all batches in grey color as background
### make dot size smaller
setGeneric("DimPlot", function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL, do.label = FALSE, label.size = 1, no.legend = FALSE) standardGeneric("DimPlot"))
#' @export
setMethod("DimPlot", "seurat",
          function(object,reduction.use="pca",dim.1=1,dim.2=2,cells.use=NULL,pt.size=1,do.return=FALSE,do.bare=FALSE,cols.use=NULL,group.by="ident",pt.shape=NULL, do.label = FALSE, label.size = 1, no.legend = FALSE) {
            cells.use=set.ifnull(cells.use,colnames(object@data))
            dim.code=translate.dim.code(reduction.use); dim.codes=paste(dim.code,c(dim.1,dim.2),sep="")
            data.plot=FetchData(object,dim.codes,cells.use)

            ident.use=as.factor(object@ident[cells.use])
            if (group.by != "ident") ident.use=as.factor(FetchData(object,group.by)[,1])
            data.plot$ident=ident.use
            x1=paste(dim.code,dim.1,sep=""); x2=paste(dim.code,dim.2,sep="")
            data.plot$x=data.plot[,x1]; data.plot$y=data.plot[,x2]
            data.plot$pt.size=pt.size
p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident),alpha=0.9),size=pt.size) + scale_colour_manual(values=myBrewerPalette) +coord_cartesian(ylim=c(-29,30),xlim=c(-46,44)) 
            if (!is.null(pt.shape)) {
              shape.val=FetchData(object,pt.shape)[cells.use,1]
              if (is.numeric(shape.val)) {
                shape.val=cut(shape.val,breaks = 5)
              }
              data.plot[,"pt.shape"]=shape.val
p=ggplot(data.plot,aes(x=x,y=y))+geom_point(aes(colour=factor(ident),shape=factor(pt.shape),size=pt.size,alpha=0.1))
            }
            if (!is.null(cols.use)) {
              p=p+scale_colour_manual(values=cols.use)
            }
            p2=p+xlab(x1)+ylab(x2) #+scale_size(range = c(pt.size, pt.size))
            p3=p2+gg.xax()+gg.yax()+no.legend.title+theme_bw() #+nogrid
            p3=p3+guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE)
            if (do.label) {
              data.plot %>% dplyr::group_by(ident) %>% summarize(x = median(x), y = median(y)) -> centers
              p3 <- p3 + geom_point(data = centers, aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, aes(label=ident), size = label.size)
            }
            if (no.legend) {
              p3 <- p3 + theme(legend.position = "none")
            }
            if (do.return) {
              if (do.bare) return(p)
              return(p3)
            }
            if (do.bare) print(p)
            else print(p3)
          }
)

sets=levels(dge@data.info$orig.ident)
data=data_bg=dge@tsne.rot
which(rownames(data)!=names(dge@ident))
data$ident=dge@ident
data %>% dplyr::group_by(ident) %>% summarize(tSNE_1 = median(tSNE_1), tSNE_2 = median(tSNE_2)) -> centers
label.size=6
plotset=NULL
for(i in 1:length(sets)){
set=sets[i]
# plot each cell type separately, adding all others as background
cells.use=names(dge@ident)[which(gsub("_.*","",names(dge@ident))==set)]
data.plot=data[cells.use,]
p=ggplot(data.plot,aes(x=tSNE_1,y=tSNE_2,color=ident))+
geom_point(data=data_bg,colour="grey",alpha=0.2,size=0.8)+
geom_point(alpha=0.8,size=0.8) + 
scale_colour_manual(values=myBrewerPalette) +
coord_cartesian(ylim=c(-29,30),xlim=c(-46,44)) +
gg.xax()+gg.yax()+no.legend.title+theme_bw()+
guides(size = FALSE,alpha=FALSE,colour=FALSE,fill=FALSE) #+theme(legend.title=element_blank())+gg.legend.pts(6)+gg.legend.text(12) # no legend by setting guides(xx=FALSE) #    # guides( colour = FALSE,fill=FALSE,
p3 <- p + geom_text(data=centers, aes(label=ident), colour="black", size = label.size)
plotset[[i]]=p3
}

### plot
jpeg(file=paste0(dgefile,"tSNE_4clusters_set_bg.jpeg"),res=300,height=4500,width=3000)
MultiPlotList(plotset,cols = 4)
dev.off()
pdf(paste(dgefile,"tSNE_4clusters_set_bg.pdf"),height=17,width=12)
MultiPlotList(plotset,cols = 4)
dev.off()

### end


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

  c1=c2=c3=c4=c5=NULL
  ncluster=length(unique(dge@ident))
## %Cell cycle genes per Cell
for(id in 1:ncluster){
	c1=rbind(c1,summary(exp1[which(dge@ident==id)]))
	c2=rbind(c2,summary(exp2[which(dge@ident==id)]))
	c3=rbind(c3,summary(exp3[which(dge@ident==id)]))
	c4=rbind(c4,summary(exp4[which(dge@ident==id)]))
	c5=rbind(c5,summary(exp5[which(dge@ident==id)]))
}
write.table(cbind(c1,c2,c3,c4,c5),paste0(dgefile,"Cluster_",setsname[j],"_cellcyclegenes.txt"),quote=F,row.names=F,col.names=F,sep="\t")


###### Finding differentially expressed genes for all SPG subtypes - Table S5B
library(Seurat)
library(FNN)
library(igraph)
library(dplyr)

### 09.25.2017 Update2: select markers based on fold change (log1.5) and p-value only
# change to min.pct=0 and min.diff.pct=0, and thresh.use=log(1.5)
dge=SetAllIdent(dge,id=res[1])
table(dge@ident)

markersall=FindAllMarkers(dge,only.pos=TRUE,min.pct=0,min.diff.pct=0,thresh.use=log(1.5),test.use="bimod",do.print = TRUE)
write.table(markersall,paste0(home,"data_DGE/MouseAdultST24_MarkersAll_UMI20cell15_31clusters_ReCluster",setsname[j],"_",res[i],"_pct0_diffpct0_thresh1.5fold_clustersPC1-14_9.25.2017.txt"),col.names=T,row.names=T,quote=F,sep="\t")

# Add mean expression level of each markers for each cluster for
#1) all cells in the cluster; and
#2) the cells positively expressing the specific marker in the cluster
table(markersallin$cluster[!grepl("V",markersallin$cluster)])
markersall=markersallin[!grepl("V",markersallin$cluster),]
markersall$cluster=as.numeric(markersall$cluster)
meanmarker=matrix(0,nrow=dim(markersall)[1],ncol=3)
meanmarker[,1]=markersall$gene
for(markerj in 1:dim(markersall)[1]){
marker=markersall[markerj,6]
cluster=markersall[markerj,5]
tmp=dge@data[marker,which(dge@ident==cluster)]
meanmarker[markerj,2]=expMean(tmp)
meanmarker[markerj,3]=expMean(tmp[tmp!=0])
print(markerj)
}
markersallmean=cbind(markersall,meanmarker)
which(markersallmean[,6] != markersallmean[,7])   # make sure nothing
write.table(markersallmean[,c(6,1:5,8,9)],paste0(home,"data_DGE/MouseAdultST24_MarkersMean_UMI20cell15_31clusters_ReCluster",setsname[j],"_",res[i],"_pct0.2_diffpct0.2_thresh2fold_clustersPC1-14_7.10.2017.txt"),col.names=T,row.names=F,quote=F,sep="\t")
# this is Table S5B

### markers compared against nearest neighbors
res="res.0.2order"
dge=SetAllIdent(dge,id=res[1])
table(dge@ident)

markers0102=FindMarkers(dge,1,2,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers0103=FindMarkers(dge,1,3,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers0203=FindMarkers(dge,2,3,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers0204=FindMarkers(dge,2,4,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")
markers0304=FindMarkers(dge,3,4,min.pct=0.2,min.diff.pct=0.2,thresh.use=log(2),test.use="bimod")

markers0102$cluster="1Vs2";markers0102$gene=rownames(markers0102)
markers0103$cluster="1Vs3";markers0103$gene=rownames(markers0103)
markers0203$cluster="2Vs3";markers0203$gene=rownames(markers0203)
markers0204$cluster="2Vs4";markers0204$gene=rownames(markers0204)
markers0304$cluster="3Vs4";markers0304$gene=rownames(markers0304)

markersneighbor=rbind(markers0102,markers0103,markers0203,markers0204,markers0304)
write.table(markersneighbor,paste0(home,"data_DGE/MarkersAll_Diff1_MouseAdultST24_UMI20cell15_31clusters_ReCluster123cellsUMI1k_res.0.2order_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=T,quote=F,sep="\t")
head(markersneighbor)
table(markersneighbor$cluster)
table(markersneighbor[markersneighbor$avg_diff>0,]$cluster)
table(markersneighbor[markersneighbor$avg_diff<0,]$cluster)
1Vs2 1Vs3 2Vs3 2Vs4 3Vs4 
  36  143   18  209   79 
  20   76    4   75   35 
  16   67   14  134   44

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
write.table(markersneighbormean[,c(6,1:5,8,9)],paste0(home,"data_DGE/MarkersAll_Diff1Mean_MouseAdultST24_UMI20cell15_31clusters_ReCluster123cellsUMI1k_res.0.2order_diff_pct0.2_diffpct0.2_thresh2fold_3.5.2018.txt"),col.names=T,row.names=F,quote=F,sep="\t")

###### Heatmap for markers across 4 SPG Subtypes - Figure 4B 
### Calculate cluster centroid
dge=SetAllIdent(dge,id="res.0.2order")
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

### Calcualte gene expression for each cluster centroid
# for each cluster, calculate average normalized expression of each gene
tmpdge=data.frame(t(as.matrix(dge@data[,cells.use])))
# make sure same order for cells.ident and dge before combining
which(names(cells.ident)!=colnames(dge@data)[cells.use])  # integer(0)
cells=data.frame(ident=cells.ident,t(as.matrix(dge@data[,cells.use])))
centroid=matrix(,dim(dge@data)[1],length(unique(cells$ident)))
rownames(centroid)=rownames(dge@data)
colnames(centroid)=unique(cells$ident)
for(i in unique(cells$ident)){
    centroid[,i]=apply(cells[which(cells$ident==i),-1],2,function(x) expMean(as.numeric(x))) # log(mean(exp(x) - 1) + 1)
    print(i)
}

### Genes Standardized Across Cell Types
centroid.std=(centroid-apply(centroid,1,mean))/apply(centroid,1,sd)
write.table(centroid.std,paste0(home,"data_DGE/SPG123UMIs1k4clusters_centroid_allgenes_std.txt"),row.names=T,col.names=T,quote=F,sep="\t")
centroid.std=read.table(paste0(home,"data_DGE/SPG123UMIs1k4clusters_centroid_allgenes_std.txt"),header=T,row.names=1)
genes=centroid.std

### plot heatmap for all markers for 4 SPG subtypes
colsep.use=cumsum(table(gsub("_.*","",levels))[levels])
col.lab=rep("",length(levels))
col.lab=gsub(".*_","",levels)

col.use=redblue100

ncluster=length(levels)
sidecol=matrix(0,2,length(levels))
sidecol[1,]=rep(rep(c("white","white"),each=12),3)[1:sum(ncluster)]
sidecol[2,]=myBrewerPalette[1:sum(ncluster)]
clab=cbind(sidecol[2,],sidecol[1,])
rlab=sidecol
rownames(rlab)=c("","Cell Type")
colnames(clab)=c("Cell Type","")

data.use=centroid.std[markersall$gene,]
row.lab=rownames(data.use)
jpeg(file=paste0(dgename,"centroid_std_markersall.jpeg"),res=300,height=1600,width=1000)
par(mar=c(4,4,1,1),mgp=c(2.5, 1, 0))
heatmap.3(data.use,dendrogram="none",Rowv=NA,Colv=NA,trace = "none",col=col.use,colsep = colsep.use,sepcolor="black",sepwidth=c(0.001,0.001),ColSideColors=clab,labCol=col.lab,labRow=row.lab,cexCol=1,cexRow=0.3,ColSideColorsSize = 3,RowSideColorsSize = 1.5,symm=F,symkey=F,symbreaks=F, scale="none",margins=c(7,3))
dev.off()
# save as Figure 4B

###### visualize expression of markers in tSNE space - Figure 4C
markers=c("Zbtb16","Gfra1","Neurog3","Nanos3","Lin28a","Stra8","Kit","Dmrtb1","Tcea3","Crabp1","Prdm9","Sycp3","Hormad1","Hormad2")

redblue100.alpha<-rgb(read.table(paste0(home,"data_DGE/redblue100.txt"),sep='\t',row.names=1,header=T),alpha=0.8)

object=dge
setname="SPG"
genes=markers

for(j in 1:length(genes)){
feature=features.plot=genes[j]

features.plot; dim.1 = 1; dim.2 = 2; cells.use = NULL; pt.size = 1;
pch.use = 16; reduction.use = "tsne";
use.imputed = FALSE; no.axes = TRUE; no.legend = FALSE

cells.use <- set.ifnull(cells.use, colnames(object@data))
dim.code <- translate.dim.code(reduction.use)
dim.codes <- paste(dim.code, c(dim.1, dim.2), sep = "")
data.plot <- object@tsne.rot

x1 <- paste(dim.code, dim.1, sep = "")
x2 <- paste(dim.code, dim.2, sep = "")

data.plot$x <- data.plot[, x1]
data.plot$y <- data.plot[, x2]
data.plot$pt.size <- pt.size
data.use <- data.frame(t(object@data[feature,]))
rownames(data.use)=feature

data.gene0 <- na.omit(data.frame(data.use[feature, ]))
data.plot$gene <- t(data.gene0)

st6<- data.plot
st6<-st6[order(st6[,6]),]
z<-st6[,6]

# edit gene name
if(features.plot=="Trdmt1"){features.plot="Dnmt2"}
if(features.plot=="Gm1564"){features.plot="Meioc"}

jpeg(paste0(dgename,setname,"_",feature,"_40-100redblue0.8.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[40:100][(z - min(z))/diff(range(z))*60 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

jpeg(paste0(dgename,setname,"_",feature,"_redblue0.8.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
z<-st6[,6]
zcolor <- redblue100.alpha[(z - min(z))/diff(range(z))*100 + 1]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()

## To maximize red color, use the same color for top 1% of expression for each marker
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

jpeg(paste0(dgename,setname,"_",feature,"_top1percent.jpeg"),height=850,width=800,res=300)
par(mar=c(0.5,0.5,2,0.5),mgp=c(1,0.5,0))
zcolor <- redblue100.alpha[c((z - min(z))/diff(range(z))*100 + 1, rep(100,top1n))]
plot(st6[,1],st6[,2], col=zcolor,pch=19,cex=0.3,cex.main=1.5,axes=F,main=features.plot,xlab="",ylab="")
dev.off()
}
# save as Figure 4C
